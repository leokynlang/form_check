# -*- coding: utf-8 -*-
"""battery_probe2.py — 修法探针:wild bootstrap(Rademacher)零分布。

探针 1 的诊断:v(y)=y 与 v(y)=1/y 两族在乘性噪声下残差异方差,
排列零分布可交换性破裂 → rational 假升级 30%。修法:零分布改为
pred + resid*sign(Rademacher)——对异方差稳健的经典修法
(wild bootstrap, Liu 1988 / Davidson&Flachaire 2008)。
同种子重跑四场景×四族,与排列版逐格对照。

预测(运行前写):
W1. H0_BASELINE power 族 wild 版仍校准(0.02-0.12,wild 在同方差下
    通常略保守);
W2. RATIONAL viable(rational) 从 0.66/0.70 → ≥0.90(修主病灶);
W3. LOG viable(logarithmic) 保持 ≥0.90,且 rational 陪跑率下降;
W4. THRESHOLD 升级触发保持 ≥0.95(wild 不应牺牲结构功效);
W5. 代价检查:EXP_MODEL viable(exp) 保持 ≥0.90。
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from scenarios import SCENARIOS  # noqa: E402

B_PERM = 999
ALPHA = 0.05
N = 200
SEEDS = 200
RELS = (0.01, 0.03)

FAMILIES = {
    "power":       (lambda x: np.log(x), lambda y: np.log(y)),
    "exponential": (lambda x: x,           lambda y: np.log(y)),
    "logarithmic": (lambda x: np.log(x), lambda y: y),
    "rational":    (lambda x: 1.0 / x,   lambda y: 1.0 / y),
}
FAM_IDX = {f: i for i, f in enumerate(FAMILIES)}

TRUTH_FAMILY = {
    "H0_BASELINE": "power",
    "EXP_MODEL": "exponential",
    "LOG_MODEL": "logarithmic",
    "RATIONAL": "rational",
    "THRESHOLD": None,
}
SCENARIO_LIST = list(TRUTH_FAMILY)


def spec_test_wild(x, y, uf, vf, seed):
    """wild 版:零分布 = pred + resid*Rademacher(异方差稳健)。"""
    u_raw = uf(x)
    v = vf(y)
    u = (u_raw - u_raw.mean()) / u_raw.std()
    base = np.column_stack([np.ones_like(u), u])
    design = np.column_stack([base, u ** 2])
    coef, _, _, _ = np.linalg.lstsq(design, v, rcond=None)
    t_obs = float(coef[-1])
    b, _, _, _ = np.linalg.lstsq(base, v, rcond=None)
    pred = base @ b
    resid = v - pred
    if float(resid @ resid) <= 1e-12 * float(v @ v):
        return {"p": np.nan, "t": np.nan, "fail": "degenerate_zero_residual"}
    rng = np.random.default_rng(seed)
    pinv_row = np.linalg.pinv(design)[-1, :]
    c2_base = float(pinv_row @ pred)
    signs = rng.choice((-1.0, 1.0), size=(B_PERM, len(v)))
    null = c2_base + (resid[None, :] * signs) @ pinv_row
    mu, sd = float(null.mean()), float(null.std())
    if not np.isfinite(sd) or sd <= 0:
        return {"p": np.nan, "t": np.nan, "fail": "degenerate_null"}
    p = (1.0 + float(np.sum(np.abs(null - mu) >= abs(t_obs - mu)))) / (B_PERM + 1.0)
    return {"p": float(p), "t": float((t_obs - mu) / sd), "fail": None}


def main():
    rows = []
    t0 = time.perf_counter()
    total = len(SCENARIO_LIST) * len(RELS) * SEEDS
    done = 0
    for scen in SCENARIO_LIST:
        fn = SCENARIOS[scen]
        for rel_i, rel in enumerate(RELS):
            for s in range(SEEDS):
                done += 1
                df, feat_str = fn(seed=700_000 + s * 97 + rel_i * 13,
                                  n=N, rel=rel)
                x = df[feat_str.split(",")[0]].values
                y = df["target"].values
                rec = {"scenario": scen, "rel": rel, "seed": s}
                viable = set()
                for fam, (uf, vf) in FAMILIES.items():
                    bs = 700_000 + s * 97 + rel_i * 13 + FAM_IDX[fam] * 1_000_003
                    r = spec_test_wild(x, y, uf, vf, bs)
                    rec[f"rej_{fam}"] = bool(r["p"] < ALPHA)
                    if np.isfinite(r["p"]) and r["p"] >= ALPHA:
                        viable.add(fam)
                rec["n_viable"] = len(viable)
                truth = TRUTH_FAMILY[scen]
                rec["containment"] = (truth in viable) if truth else None
                rec["escalation_fired"] = (len(viable) == 0)
                rows.append(rec)
                if done % 500 == 0:
                    print(f"[{done}/{total}] ({time.perf_counter()-t0:.0f}s)",
                          flush=True)

    raw = pd.DataFrame(rows)
    out_dir = os.path.join(HERE, "runs", "battery_probe_wild")
    os.makedirs(out_dir, exist_ok=True)
    raw.to_parquet(os.path.join(out_dir, "raw.parquet"), index=False)

    print("\n=== wild 版拒绝率(200 种子)===")
    summary = {}
    for scen in SCENARIO_LIST:
        for rel in RELS:
            g = raw[(raw.scenario == scen) & (raw.rel == rel)]
            rates = {f: float(g[f"rej_{f}"].mean()) for f in FAMILIES}
            cont = (float(g["containment"].dropna().astype(float).mean())
                    if g["containment"].notna().any() else None)
            esc = float(g.escalation_fired.mean())
            med_v = float(g.n_viable.median())
            summary[f"{scen}@rel={rel}"] = {"rej": rates, "containment": cont,
                                            "escalation": esc,
                                            "median_n_viable": med_v}
            print(f"{scen:12s} rel={rel:.2f} "
                  + " ".join(f"{f}:{rates[f]:.3f}" for f in FAMILIES)
                  + f" | med_v={med_v:.0f} cont="
                  + (f"{cont:.3f}" if cont is not None else "n/a")
                  + f" esc={esc:.3f}")

    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump({"summary": summary,
                   "config": {"null": "wild_rademacher", "B_PERM": B_PERM,
                              "ALPHA": ALPHA, "N": N, "SEEDS": SEEDS,
                              "RELS": list(RELS)}},
                  f, ensure_ascii=False, indent=1)
    print(f"\nDONE {time.perf_counter()-t0:.0f}s -> {out_dir}")


if __name__ == "__main__":
    main()
