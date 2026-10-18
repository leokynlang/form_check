"""Read-only, machine-readable inventory of the R2 experimental record."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
KEYS = ["task", "n", "rel", "seed"]
INPUTS = [
    "runs/srsd_full/raw_srsd.parquet",
    "runs/srsd_full/quality_recompute.parquet",
    "runs/cd_svd/s2_cd_svd.csv",
    "runs/gate_cmp/results.csv",
]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="runs/p01_inventory")
    args = parser.parse_args()
    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    frames, records = {}, {}
    for name in INPUTS:
        path = ROOT / name
        df = (pd.read_parquet(path) if path.suffix == ".parquet"
              else pd.read_csv(path, low_memory=False))
        frames[name] = df
        key = KEYS if "rel" in df else ["bench", "task", "n", "seed"]
        duplicates = int(df.duplicated(key).sum())
        if duplicates:
            raise ValueError(f"Duplicate keys in {name}: {duplicates}")
        records[name] = {
            "sha256": digest(path), "bytes": path.stat().st_size,
            "rows": len(df), "tasks": int(df.task.nunique()),
            "columns": df.columns.tolist(), "keys": key,
            "duplicate_keys": duplicates,
            "n_grid": sorted(int(x) for x in df.n.unique()),
            "noise_grid": sorted(df.rel.unique().tolist()) if "rel" in df else [0.03],
        }
    raw, quality, svd, s5 = [frames[name] for name in INPUTS]
    joins = {}
    for name, other in [("quality", quality), ("svd", svd)]:
        merged = raw[KEYS].merge(other[KEYS], on=KEYS, how="outer",
                                 validate="one_to_one", indicator=True)
        joins[name] = merged._merge.value_counts().to_dict()
        if not merged._merge.eq("both").all():
            raise ValueError(f"Incomplete {name} join")
    s5_srsd = s5[s5.bench == "srsd"]
    old_seeds = set(raw.base_seed.astype(int))
    s5_seeds = set(s5_srsd.seed.astype(int))
    result = {
        "scope": "Inventory only; no regenerated data or performance estimates.",
        "inputs": records, "s2_one_to_one_joins": joins,
        "s2_s5_dataset_seed_overlap": len(old_seeds & s5_seeds),
        "s2_positive_noise_rows": int(raw.rel.gt(0).sum()),
        "s2_task_classes": raw[["task", "klass"]].drop_duplicates().to_dict("records"),
        "svd_false_sentinel_rows": int(svd.sentinel_ok.eq(False).sum()),
        "svd_false_sentinel_tasks": sorted(svd.loc[svd.sentinel_ok.eq(False), "task"].unique()),
        "warnings": [
            "S5 cannot be joined to S2 on dataset identity: different actual seeds.",
            "S2 cached C/D original columns are historical Gram-path results.",
            "Old sample-dependent endpoint is not a fixed extrapolation domain.",
            "Already studied equations are not untouched confirmation tasks.",
        ],
    }
    target = out / "inventory.json"
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"output": str(target), "s2_rows": len(raw),
                      "s5_rows": len(s5), "joins_complete": True,
                      "s2_s5_seed_overlap": len(old_seeds & s5_seeds)}))


if __name__ == "__main__":
    main()
