# -*- coding: utf-8 -*-
"""run_battery.py — 电池正式跑:合成 19 场景 + SRSD-Feynman H1 池。

配置冻结(与 P1/P2 同纪律):α=0.05, B=999 wild, n=200, 500 种子,
rel∈{0.5%,1%,3%,5%};SRSD 侧沿用官方域(homo)+hetero 两噪声型,
不跑 range 缩放(P2 已证结论对采样范围稳健,省 2/3 算力)。
种子纪律:base=900_000_000 + cell_idx*10_000_000 + seed*97,格间不重叠。

断点续跑:结果按 cell 追加到 runs/<tag>/raw_battery.parquet,
重跑同 tag 自动跳过已完成的 cell(cell 键=场景/方程×rel×noise)。
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from scenarios import SCENARIOS, SCENARIO_META  # noqa: E402
from srsd_data import registry_with_meta, make_dataset_fn  # noqa: E402
from battery import battery_run, SYNTHETIC_TRUTH, ALPHA, B_WILD  # noqa: E402

FROZEN = {"n_train": 200, "n_seeds": 500,
          "rel_grid": [0.005, 0.01, 0.03, 0.05], "alpha": ALPHA,
          "b_wild": B_WILD}


def synth_cells():
    """合成格:全部 19 场景 × rel。合成侧无 hetero 维度——噪声型
    (HEAVYTAIL/VAR_GROWING/AUTOCORR/OUTLIER*)已编为独立场景;hetero
    参数只属于 SRSD 采样接口。"""
    cells = []
    for name in SCENARIO_META:
        for rel in FROZEN["rel_grid"]:
            cells.append(("synth", name, rel, False))
    return cells


def srsd_cells():
    """SRSD 格:P2 选出的 H1 入检方程 × rel × {homo,hetero}。

    复用 P2 的分层选择结果(selection_table.csv),不重新选——
    选择先于检验结果,冻结复用。
    """
    sel_path = os.path.join(HERE, "runs", "feynman_families",
                            "selection_table.csv")
    if not os.path.exists(sel_path):
        print(f"警告: {sel_path} 不存在,SRSD 侧跳过(可先跑 P2 或单独跑合成)")
        return []
    sel = pd.read_csv(sel_path)
    eqs = sel[(sel.eligible) & (sel.klass == "H1")]["eq_name"].tolist()
    cells = []
    for eq in eqs:
        for rel in FROZEN["rel_grid"]:
            for het in (False, True):
                cells.append(("srsd", eq, rel, het))
    return cells


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="battery_full")
    ap.add_argument("--seeds", type=int, default=FROZEN["n_seeds"])
    ap.add_argument("--only", choices=("synth", "srsd"), default=None,
                    help="只跑一侧(先点火合成,SOS 侧继续补)")
    args = ap.parse_args()

    run_dir = os.path.join(HERE, "runs", args.tag)
    os.makedirs(run_dir, exist_ok=True)
    raw_path = os.path.join(run_dir, "raw_battery.parquet")

    cells = synth_cells()
    if args.only in (None, "srsd"):
        cells += srsd_cells()
    if args.only == "srsd":
        cells = srsd_cells()
    elif args.only == "synth":
        cells = synth_cells()

    done_cells = set()
    if os.path.exists(raw_path):
        old = pd.read_parquet(raw_path)
        done_cells = set(old["cell_id"].unique())
        print(f"续跑:已有 {len(done_cells)} 格完成,共 {len(cells)} 格")

    cfg = {"frozen": FROZEN, "n_cells": len(cells),
           "seed_rule": "900_000_000 + cell_idx*10_000_000 + seed*97",
           "truth_registry_synth": {k: (v or "none") for k, v
                                    in SYNTHETIC_TRUTH.items()},
           "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1)

    seeds = list(range(args.seeds))
    buf = []
    t0 = time.perf_counter()
    n_run = 0
    for ci, (side, name, rel, het) in enumerate(cells):
        cell_id = f"{side}:{name}:rel={rel}:het={int(het)}"
        if cell_id in done_cells:
            continue
        if side == "synth":
            fn = SCENARIOS[name]
            truth = SYNTHETIC_TRUTH.get(name) or "none"
        else:
            rec = None
            for r in registry_with_meta():
                if r["eq_name"] == name:
                    rec = r
                    break
            if rec is None:
                print(f"跳过未知方程 {name}")
                continue
            fn = make_dataset_fn(rec["cls"])
            truth = "none"          # H1 方程无电池族(那是 P2 语义,不是电池语义)
        rows = []
        for s in seeds:
            bs = 900_000_000 + ci * 10_000_000 + s * 97
            try:
                if side == "synth":
                    # 合成场景签名 fn(seed,n,rel)——噪声型编在场景名里,
                    # 没有 hetero 参数( hetero 只属于 SRSD 采样接口)。
                    df, feat_str = fn(seed=bs, n=FROZEN["n_train"], rel=rel)
                else:
                    df, feat_str = fn(seed=bs, n=FROZEN["n_train"], rel=rel,
                                      hetero=het)
                feats = feat_str.split(",")
            except ValueError as e:
                rows.append({"cell_id": cell_id, "side": side, "name": name,
                             "rel": rel, "hetero": het, "truth_family": truth,
                             "seed": s, "data_fail": str(e)})
                continue
            rec_row = {"cell_id": cell_id, "side": side, "name": name,
                       "rel": rel, "hetero": het, "truth_family": truth,
                       "seed": s, "data_fail": None}
            rec_row.update(battery_run(df, "target", feats, seed=bs))
            rows.append(rec_row)
        buf.append(pd.DataFrame(rows))
        n_run += 1
        if len(buf) >= 4:
            _append(raw_path, buf)
            buf = []
        print(f"[{ci+1}/{len(cells)}] {cell_id} "
              f"({time.perf_counter()-t0:.0f}s)", flush=True)
    if buf:
        _append(raw_path, buf)
    cfg["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    cfg["cells_run"] = n_run
    cfg["wall_seconds"] = time.perf_counter() - t0
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1)
    print(f"DONE {args.tag}: {n_run} 格 -> {run_dir}")


def _append(raw_path, buf):
    new = pd.concat(buf, ignore_index=True)
    if os.path.exists(raw_path):
        old = pd.read_parquet(raw_path)
        new = pd.concat([old, new], ignore_index=True)
    new.to_parquet(raw_path, index=False)


if __name__ == "__main__":
    main()
