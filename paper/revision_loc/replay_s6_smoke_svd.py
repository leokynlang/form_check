# -*- coding: utf-8 -*-
"""replay_s6_smoke_svd.py — v5 重跑第二部分:S6 与工程台架 E 的 C/D SVD 重放。

机制:monkey-patch common._ols 不动(A/B 哨兵仍走原路径对照存档),只把
reset_test 的 HC3 协方差来源换成 SVD 因子——通过替换 common.reset_test 的
绑定实现为 rerun_cd_svd.reset_test_svd 的包装。wb_run/smoke_test 通过
`import common` 取 run_all_methods,故 patch 必须落在 common 模块属性上。

  - reset_test(variant, vcov='hc3') → SVD hat + SVD 协方差(paper 定位实验
    证明等效于修复 Gram 截断问题,且良态任务与原实现 3 位一致)
  - vcov='classical' 保持原实现(B 列与存档逐位可比,兼作哨兵)

输出:
  runs/cd_svd/s6_results_svd.csv     (wb_run 重放 results)
  runs/cd_svd/e_smoke_results_svd.csv(smoke_test 重放 results)
哨兵:重放完成后与原 runs/wb_real/results.csv、real_data/smoke_results.csv
的 p_fitted_f / p_fc 列逐行比对(应为 0 漂移,因 A/B 代码路径未动)。
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, r"D:\form_check")
sys.path.insert(0, r"D:\form_check\real_data")

import common                                    # noqa: E402
from rerun_cd_svd import reset_test_svd          # noqa: E402

HERE = r"D:\form_check"
OUT_DIR = os.path.join(HERE, "runs", "cd_svd")
ORIG_RESET = common.reset_test


def patched_reset(Xlog, ly, powers=(2, 3), variant="fitted", vcov="classical"):
    if vcov == "classical":
        return ORIG_RESET(Xlog, ly, powers=powers, variant=variant,
                          vcov="classical")
    return reset_test_svd(Xlog, ly, powers=powers, variant=variant)


def replay_wb():
    import wb_run
    wb_run.OUT = OUT_DIR                          # 结果落独立目录
    wb_run.common = common                        # 确保同一模块对象
    print("=== S6 重放(wb_run.main,缓存数据,零网络) ===", flush=True)
    wb_run.main()
    src = os.path.join(OUT_DIR, "results.csv")
    df = pd.read_csv(src, encoding="utf-8-sig")
    df.to_csv(os.path.join(OUT_DIR, "s6_results_svd.csv"),
              index=False, encoding="utf-8-sig")
    return df


def replay_smoke():
    import smoke_test
    print("=== 工程台架 E 重放(smoke_test.main,clean/ 缓存) ===", flush=True)
    out_path = os.path.join(OUT_DIR, "e_smoke_results_svd.csv")
    sys.argv = ["smoke_test.py", "--max-n", "800", "--out", out_path]
    smoke_test.main()
    return pd.read_csv(out_path)


def sentinel_check(new_csv, old_csv, keys, pcols=("p_fitted_f", "p_fc")):
    new = pd.read_csv(new_csv)
    old = pd.read_csv(old_csv, encoding="utf-8-sig")
    merged = new.merge(old, on=keys, suffixes=("_new", "_old"))
    report = {}
    for pc in pcols:
        cn, co = f"{pc}_new", f"{pc}_old"
        if cn not in merged or co not in merged:
            report[pc] = "cols missing"
            continue
        a = pd.to_numeric(merged[cn], errors="coerce")
        b = pd.to_numeric(merged[co], errors="coerce")
        m = a.notna() & b.notna()
        if m.sum() == 0:
            report[pc] = "no overlap"
            continue
        maxdiff = float(np.max(np.abs(a[m] - b[m])))
        report[pc] = {"n_compared": int(m.sum()), "max_abs_diff": maxdiff}
    return report


def main():
    common.reset_test = patched_reset
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "s6"):
        replay_wb()
    if which in ("all", "e"):
        replay_smoke()
    common.reset_test = ORIG_RESET
    # 哨兵:A/B 未动 → 与存档零漂移
    print("=== 哨兵比对 ===")
    r1 = sentinel_check(
        os.path.join(OUT_DIR, "s6_results_svd.csv"),
        os.path.join(HERE, "runs", "wb_real", "results.csv"),
        keys=["year", "task", "seed"])
    print("S6 :", json.dumps(r1, ensure_ascii=False))
    r2 = sentinel_check(
        os.path.join(OUT_DIR, "e_smoke_results_svd.csv"),
        os.path.join(HERE, "real_data", "smoke_results.csv"),
        keys=["dataset", "target", "features"])
    print("E  :", json.dumps(r2, ensure_ascii=False))


if __name__ == "__main__":
    main()
