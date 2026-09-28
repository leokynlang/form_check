"""smoke_test.py — 对 clean/ 全部候选幂律关系跑 form_check vs RESET 冒烟。

目的: 验证 39 个真实工程数据集能直接喂进 common.run_all_methods(A=fc 排列
p + B/C/D=RESET 三变体), 并给出第一眼判定差异。这不是功效/校准研究
(真实数据无已知真值), 只证明数据-管线兼容并留痕首轮判定。

用法:  python real_data/smoke_test.py [--max-n 800]
写入:  real_data/smoke_results.csv
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import common  # noqa: E402

CLEAN = os.path.join(HERE, "clean")


def prepare_sub(df, target, features, max_n, seed=42):
    sub = df[[target] + features].dropna()
    sub = sub[(sub[target] > 0)
              & (sub[features] > 0).all(axis=1)].reset_index(drop=True)
    if max_n and len(sub) > max_n:
        sub = sub.sample(n=max_n, random_state=seed).reset_index(drop=True)
    return sub


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=800,
                    help="单格抽样上限(排列检验耗时), 0=不抽样")
    ap.add_argument("--out", default=os.path.join(HERE, "smoke_results.csv"),
                    help="结果输出路径(默认原位覆盖 real_data/smoke_results.csv;"
                         "重放/复算请显式指向别处,避免覆盖已审计档案)")
    args = ap.parse_args()

    with open(os.path.join(CLEAN, "manifest.json"), encoding="utf-8") as f:
        manifest = json.load(f)

    rows = []
    for name, entry in sorted(manifest.items()):
        df = pd.read_csv(os.path.join(CLEAN, f"{name}.csv"))
        for reln in entry["candidate_relations"]:
            target, features = reln["target"], reln["features"]
            if not set([target] + features).issubset(df.columns):
                print(f"  ! {name}: 列缺失, 跳过 {target}~{features}")
                continue
            sub = prepare_sub(df, target, features, args.max_n)
            base = {
                "dataset": name, "target": target,
                "features": "+".join(features),
                "n_raw": int(len(df)), "n_used": int(len(sub)),
                "note": reln.get("note", ""),
            }
            if len(sub) < 10:
                base.update(status="insufficient_n")
                rows.append(base)
                print(f"  - {name} {target}~{'+'.join(features)}: "
                      f"n={len(sub)} 跳过")
                continue
            try:
                base.update(common.run_all_methods(
                    sub, target, features, seed=42))
                base["status"] = "ok"
            except Exception as exc:  # noqa: BLE001 — 留痕不抛
                base.update(status=f"runner_error: {exc!r}")
            rows.append(base)
            ptxt = "/".join(
                f"{base.get(f'p_{k}', float('nan')):.3f}"
                for k in ("fc", "fitted_f"))
            print(f"  ok {name}: {target}~{'+'.join(features)} "
                  f"n={base['n_used']} p_fc={ptxt}")

    res = pd.DataFrame(rows)
    out = args.out
    res.to_csv(out, index=False, encoding="utf-8-sig")

    ok = res[res["status"] == "ok"]
    print(f"\n=== 冒烟汇总: {len(ok)}/{len(res)} 格成功 ===")
    for k in ("fc", "fitted_f", "fitted_hc3", "exog_hc3"):
        col = f"rej_{k}"
        if col in ok:
            print(f"  {k}: 拒绝率 {ok[col].mean():.3f}, "
                  f"失败 {ok[f'fail_{k}'].notna().sum()} 格")
    print(f"done -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
