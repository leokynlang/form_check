"""run_experiment.py — form_check vs RESET 网格实验 runner。

用法:
  python run_experiment.py --tag pilot --seeds 50
  python run_experiment.py --tag full  --seeds 500

固定网格: n ∈ {50,200,1000} × rel ∈ {0,0.5%,1%,3%,5%} × 全部 19 场景。
base_seed = 100_000 + s*97 + n_idx*13 + rel_idx*5_000_000(网格位置入种子,
保证任何网格组合可独立复现)。逐种子落盘(parquet 每 250 行 flush),
中断可续(已存在的 (scenario,n,rel,seed) 跳过)。
H0/H1 分类在 scenarios.py::SCENARIO_META 预注册,本脚本不据结果调整。
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from common import (
    ALPHA,
    capture_env,
    extrapolation_error,
    fitting_quality,
    run_all_methods,
)

from scenarios import SCENARIO_META, SCENARIOS

HERE = os.path.dirname(os.path.abspath(__file__))
N_GRID = [50, 200, 1000]
REL_GRID = [0.0, 0.005, 0.01, 0.03, 0.05]
EXTRAP_REL_MIN = 0.005       # 零噪声场景无外推意义(残差退化)


def base_seed(s, n_idx, rel_idx):
    return 100_000 + s * 97 + n_idx * 13 + rel_idx * 5_000_000


def ctx_for(df, feats):
    lx = np.log(df[feats].values)
    ly = np.log(df["target"].values)
    Xd = np.column_stack([np.ones(len(ly)), lx])
    b = np.linalg.lstsq(Xd, ly, rcond=None)[0]
    pred = Xd @ b
    two = lx.shape[1] > 1
    ctx = {"m1": lx[:, 0].mean(), "m2": lx[:, 1].mean() if two else 0.0,
           "s1": lx[:, 0].std(), "med": np.median(lx, axis=0),
           "cross_sd": float(((lx[:, 0] - np.median(lx[:, 0]))
                              * (lx[:, 1] - np.median(lx[:, 1]))).std())
           if two else 1.0}
    return ctx, pred


def one_cell(scen, n, rel, seeds, run_dir):
    fn = SCENARIOS[scen]
    rows, extrap_rows = [], []
    t_cell = time.perf_counter()
    for s in seeds:
        bs = base_seed(s, N_GRID.index(n), REL_GRID.index(rel))
        df, feat_list = fn(seed=bs, n=n, rel=rel)
        feats = feat_list.split(",")
        row = {"scenario": scen, "klass": SCENARIO_META[scen][0],
               "family": SCENARIO_META[scen][1], "n": n, "rel": rel,
               "seed": s, "base_seed": bs}
        row.update(run_all_methods(df, "target", feats, seed=bs, n_boot=200))
        row.update(fitting_quality(np.log(df[feats].values),
                                   np.log(df["target"].values)))
        rows.append(row)
        if rel >= EXTRAP_REL_MIN:
            truth_fn = SCENARIO_META[scen][2]
            ctx, _ = ctx_for(df, feats)
            med, farpt = extrapolation_error(
                df, "target", feats,
                lambda X, _fn=truth_fn, _ctx=ctx: _fn(X, _ctx))
            er = {"scenario": scen, "n": n, "rel": rel, "seed": s,
                  "extrap_med_relerr": med, "extrap_far_relerr": farpt,
                  "rej_fc": row["rej_fc"], "p_fc": row["p_fc"],
                  "rej_fitted_f": row["rej_fitted_f"], "p_fitted_f": row["p_fitted_f"],
                  "rej_fitted_hc3": row["rej_fitted_hc3"], "p_fitted_hc3": row["p_fitted_hc3"],
                  "rej_exog_hc3": row["rej_exog_hc3"], "p_exog_hc3": row["p_exog_hc3"]}
            er["alarm_any"] = bool(row["rej_fc"] or row["rej_fitted_f"]
                                   or row["rej_fitted_hc3"] or row["rej_exog_hc3"])
            extrap_rows.append(er)
    return pd.DataFrame(rows), pd.DataFrame(extrap_rows), time.perf_counter() - t_cell


def load_existing(run_dir):
    p = os.path.join(run_dir, "raw.parquet")
    if os.path.exists(p):
        df = pd.read_parquet(p)
        return df, os.path.join(run_dir, "extrap.parquet")
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--scenarios", default="all",
                    help="逗号分隔场景名,默认全部")
    ap.add_argument("--rel", default="all", help="逗号分隔 rel 值,默认全部")
    ap.add_argument("--n", default="all", help="逗号分隔 n 值,默认全部")
    args = ap.parse_args()

    run_dir = os.path.join(HERE, "runs", args.tag)
    os.makedirs(run_dir, exist_ok=True)
    scen_names = (list(SCENARIOS) if args.scenarios == "all"
                  else args.scenarios.split(","))
    rel_grid = (REL_GRID if args.rel == "all"
                else [float(x) for x in args.rel.split(",")])
    n_grid = (N_GRID if args.n == "all"
              else [int(x) for x in args.n.split(",")])

    cfg = {"tag": args.tag, "seeds": args.seeds, "n_grid": N_GRID,
           "rel_grid": REL_GRID, "scenarios": scen_names, "alpha": ALPHA,
           "seed_rule": "100_000 + s*97 + n_idx*13 + rel_idx*5_000_000",
           "preregistration": (
               "H0/H1 分类由 scenarios.py::SCENARIO_META 在代码写入时固定,"
               "不依据检验结果调整;方法 A=dimreflex.form_check 本体+包装层"
               "排列经验 p 值;B/C/D=RESET fitted-F/fitted-HC3/exog-HC3;"
               "全部在 alpha=0.05 上比较。"),
           "scenario_meta": {k: {"klass": v[0], "family": v[1], "desc": v[3]}
                             for k, v in SCENARIO_META.items()},
           "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1)
    with open(os.path.join(run_dir, "env.json"), "w", encoding="utf-8") as f:
        json.dump(capture_env(), f, ensure_ascii=False, indent=1)

    seeds = list(range(args.seeds))
    total_cells = len(scen_names) * len(n_grid) * len(rel_grid)
    done = 0
    t0 = time.perf_counter()
    old_raw, _old_ex = load_existing(run_dir)
    buf, ebuf = [], []
    for scen in scen_names:
        for n in n_grid:
            for rel in rel_grid:
                # 断点续跑:该格已齐则跳过
                if old_raw is not None:
                    m = ((old_raw.scenario == scen) & (old_raw.n == n)
                         & (old_raw.rel == rel))
                    if int(m.sum()) >= args.seeds:
                        done += 1
                        continue
                raw, ex, dt = one_cell(scen, n, rel, seeds, run_dir)
                buf.append(raw)
                if len(ex):
                    ebuf.append(ex)
                done += 1
                print(f"[{done:3d}/{total_cells}] {scen:16s} n={n:4d} "
                      f"rel={rel:.3f} {dt:6.1f}s "
                      f"(elapsed {time.perf_counter()-t0:7.1f}s)", flush=True)
                if len(buf) >= 4:                      # 每 4 格 flush 一次
                    _flush(run_dir, buf, ebuf)
                    buf, ebuf = [], []
    if buf:
        _flush(run_dir, buf, ebuf)
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        cfg["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
        cfg["wall_seconds"] = time.perf_counter() - t0
        json.dump(cfg, f, ensure_ascii=False, indent=1)
    print(f"DONE {args.tag}: {total_cells} cells, "
          f"{time.perf_counter()-t0:.0f}s -> {run_dir}")


def _flush(run_dir, buf, ebuf):
    raw = pd.concat(buf, ignore_index=True)
    rp = os.path.join(run_dir, "raw.parquet")
    if os.path.exists(rp):
        old = pd.read_parquet(rp)
        raw = pd.concat([old, raw], ignore_index=True)
    raw.to_parquet(rp, index=False)
    raw.to_csv(os.path.join(run_dir, "raw.csv"), index=False)
    if ebuf:
        ex = pd.concat(ebuf, ignore_index=True)
        ep = os.path.join(run_dir, "extrap.parquet")
        if os.path.exists(ep):
            old = pd.read_parquet(ep)
            ex = pd.concat([old, ex], ignore_index=True)
        ex.to_parquet(ep, index=False)


if __name__ == "__main__":
    main()
