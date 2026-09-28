"""run_srsd.py — SRSD-Feynman 主基准 runner(任务2)。

任务池: srsd_data.registry_with_meta() 中 klass∈{H0,H1} 的公式(32 个,
13 H0 单项式 + 19 H1 非单项式),按官方采样范围/划分协议生成数据。
网格: n ∈ {50,200,1000} × rel ∈ {0,0.5%,1%,3%,5%} × seeds。
每种子每任务:
  - 训练集(n 个样本,含噪声) → 四方法 A-D;
  - test 集(2000 样本,官方域,同 rel 噪声)→ log 空间测试 RMSE(同分布外);
  - 外推集(域外 ×3,沿第一主变量,其余中位)→ 中位相对误差;
  - 记录 form_check 报警与三种 RESET 拒绝 → 报警-误差关联。
失败登记: 每行 fail_* 字段 + data_fail 列,绝不静默丢弃。
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from common import ALPHA, _ols, capture_env, run_all_methods
from srsd_data import make_dataset_fn, registry_with_meta

HERE = os.path.dirname(os.path.abspath(__file__))
N_GRID = [50, 200, 1000]
REL_GRID = [0.0, 0.005, 0.01, 0.03, 0.05]
N_TEST = 2000
FAR = 3.0                       # 外推 = 训练 log 域边界外推 ×3(半个十进位)
EXTRAP_REL_MIN = 0.005


def base_seed(s, n_idx, rel_idx):
    return 200_000 + s * 97 + n_idx * 13 + rel_idx * 5_000_000


def eval_errors(fn, eq_func, n, rel, train_df, feats, test_seed):
    """test RMSE(同分布)+ 外推中位相对误差。返回 dict。"""
    Xtr = np.column_stack([np.ones(len(train_df)),
                           np.log(train_df[feats].values)])
    b, _, _, _, _ = _ols(Xtr, np.log(train_df["target"].values))
    out = {"test_log_rmse": np.nan, "extrap_med_relerr": np.nan,
           "extrap_far_relerr": np.nan}
    # --- test 集:官方域重抽(独立流),同 rel 噪声 ---
    try:
        dft, _ = fn(seed=test_seed, n=N_TEST, rel=rel)
        Xt = np.column_stack([np.ones(len(dft)), np.log(dft[feats].values)])
        resid = np.log(dft["target"].values) - Xt @ b
        out["test_log_rmse"] = float(np.sqrt(np.mean(resid ** 2)))
    except Exception as e:  # noqa: BLE001 — 测试集生成失败写 data_fail,逐任务留痕
        out["test_log_rmse"] = np.nan
        out["data_fail"] = f"test_gen:{type(e).__name__}"
    # --- 外推:沿第一特征,log 域边界外扩 ln(FAR) ---
    try:
        Xlog_tr = np.log(train_df[feats].values)
        lo, hi = Xlog_tr.min(axis=0), Xlog_tr.max(axis=0)
        dlog = float(np.log(FAR))
        m = 12
        gx = np.concatenate([np.linspace(lo[0] - dlog, lo[0], m // 2),
                             np.linspace(hi[0], hi[0] + dlog, m - m // 2)])
        Xe = np.tile(np.median(Xlog_tr, axis=0), (len(gx), 1))
        Xe[:, 0] = gx
        # 真值:闭式 eq_func(原空间)→ 一次 log;relerr 用 log 差的稳定形式
        xs = np.exp(Xe)
        with np.errstate(all="ignore"):
            ye_true_log = np.log(np.array([float(eq_func(xs[i]))
                                           for i in range(len(xs))]))
            ye_pred = np.column_stack([np.ones(len(Xe)), Xe]) @ b
            relerr = np.abs(np.exp(ye_pred - ye_true_log) - 1.0)
        out["extrap_med_relerr"] = float(np.median(relerr))
        out["extrap_far_relerr"] = float(relerr[-1])
    except Exception as e:  # noqa: BLE001 — 外推评估失败写 data_fail,逐任务留痕
        out["data_fail"] = out.get("data_fail", "") + f"|extrap:{type(e).__name__}"
    return out


def one_cell(rec, n, rel, seeds):
    fn = make_dataset_fn(rec["cls"])
    inst = rec["cls"]()
    rows = []
    t_cell = time.perf_counter()
    for s in seeds:
        bs = base_seed(s, N_GRID.index(n), REL_GRID.index(rel))
        row = {"task": rec["eq_name"], "klass": rec["klass"], "n": n,
               "rel": rel, "seed": s, "base_seed": bs}
        try:
            df, feat_str = fn(seed=bs, n=n, rel=rel)
            feats = feat_str.split(",")
        except Exception as e:  # noqa: BLE001 — 训练集生成失败写 data_fail,该种子跳过
            row["data_fail"] = f"train_gen:{type(e).__name__}"
            rows.append(row)
            continue
        row["data_fail"] = None
        row.update(run_all_methods(df, "target", feats, seed=bs, n_boot=200))
        er = eval_errors(fn, inst.eq_func, n, rel, df, feats,
                         test_seed=bs + 777_777)
        row.update(er)
        rows.append(row)
    return pd.DataFrame(rows), time.perf_counter() - t_cell


def _flush(run_dir, buf):
    raw = pd.concat(buf, ignore_index=True)
    rp = os.path.join(run_dir, "raw_srsd.parquet")
    if os.path.exists(rp):
        old = pd.read_parquet(rp)
        raw = pd.concat([old, raw], ignore_index=True)
    raw.to_parquet(rp, index=False)
    raw.to_csv(os.path.join(run_dir, "raw_srsd.csv"), index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--rel", default="all")
    ap.add_argument("--n", default="all")
    ap.add_argument("--tasks", default="all")
    args = ap.parse_args()

    run_dir = os.path.join(HERE, "runs", args.tag)
    os.makedirs(run_dir, exist_ok=True)
    reg = [r for r in registry_with_meta() if r["klass"] in ("H0", "H1")]
    if args.tasks != "all":
        want = set(args.tasks.split(","))
        reg = [r for r in reg if r["eq_name"] in want]
    rel_grid = (REL_GRID if args.rel == "all"
                else [float(x) for x in args.rel.split(",")])
    n_grid = (N_GRID if args.n == "all"
              else [int(x) for x in args.n.split(",")])

    cfg = {"tag": args.tag, "seeds": args.seeds, "n_grid": N_GRID,
           "rel_grid": REL_GRID, "n_test": N_TEST, "far": FAR,
           "alpha": ALPHA,
           "seed_rule": "200_000 + s*97 + n_idx*13 + rel_idx*5_000_000",
           "task_pool": {r["eq_name"]: {"klass": r["klass"],
                                        "reason": r["reason"],
                                        "n_vars": r["n_vars"],
                                        "sampling": r["sampling_ranges"],
                                        "sympy": r["sympy_str"]}
                         for r in reg},
           "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1, default=str)
    with open(os.path.join(run_dir, "env.json"), "w", encoding="utf-8") as f:
        json.dump(capture_env(), f, ensure_ascii=False, indent=1)

    seeds = list(range(args.seeds))
    total = len(reg) * len(n_grid) * len(rel_grid)
    done = 0
    t0 = time.perf_counter()
    buf = []
    rp = os.path.join(run_dir, "raw_srsd.parquet")
    existing = pd.read_parquet(rp) if os.path.exists(rp) else None
    for rec in reg:
        for n in n_grid:
            for rel in rel_grid:
                done += 1
                if existing is not None:
                    m = ((existing.task == rec["eq_name"]) & (existing.n == n)
                         & (existing.rel == rel))
                    if int(m.sum()) >= args.seeds:
                        continue
                raw, dt = one_cell(rec, n, rel, seeds)
                buf.append(raw)
                print(f"[{done:3d}/{total}] {rec['eq_name']:22s} {rec['klass']:3s} "
                      f"n={n:4d} rel={rel:.3f} {dt:6.1f}s "
                      f"(elapsed {time.perf_counter()-t0:7.1f}s)", flush=True)
                if len(buf) >= 3:
                    _flush(run_dir, buf)
                    buf = []
    if buf:
        _flush(run_dir, buf)
    cfg["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    cfg["wall_seconds"] = time.perf_counter() - t0
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1, default=str)
    print(f"DONE {args.tag}: {total} cells -> {run_dir}")


if __name__ == "__main__":
    main()
