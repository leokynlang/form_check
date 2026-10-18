"""run_pysr_gate.py — 端到端：PySR 产出幂律候选 → 门卫检验。

补论文 §6 的缺口："no symbolic-regression method is run anywhere
（…）not the pass rate of method-emitted wrong forms"。

设计（可辩护）：
- PySR 在 **log 空间**搜索，算子集限制为 {'+', '-'}。
  于是 log(y) = c0 + Σ aᵢ·log(xᵢ) 是它唯一能表达的形式，
  这恰是幂律 y = e^{c0}·Π xᵢ^{aᵢ} 的对数形式。
  即：PySR 的假设类 = 门卫的零假设类，候选天然落在可评价空间内。
- 方法决定的是变量选择与指数/系数，这正是受限 SR 管线的自由度。

协议对齐现有 S2/S7：
- 只用 S2 中幂律可近似良好的 H1 任务（避免 R² 被结构性压死）。
- 门卫 A（form_check 置换）、B/C/D（RESET，C/D 走 SVD 路径）。
- 记录 R² 与外推误差，按论文的 R²≥0.99 与 far>10% 两条门槛分层。
"""
from __future__ import annotations

import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from common import form_check_with_p, reset_test  # noqa: E402
import srsd_data  # noqa: E402

OUT = os.path.join(HERE, "runs", "pysr_gate")
os.makedirs(OUT, exist_ok=True)

# 任务分两组：
#  DANGEROUS 组 —— S2 实测「拟合好但外推危险」的来源（rmse<0.1 且 far>10%
#    的样本只集中在 ii.38.14 / bonus.15 / i.30.5；i.30.5 已 withdraw）
#  SAFE 组 —— 幂律可近似良好且外推安全，作为对照
DANGEROUS_TASKS = ["FeynmanIICh38Eq14", "FeynmanBonus15"]
SAFE_TASKS = [
    "FeynmanIIICh4Eq32", "FeynmanIIICh4Eq33", "FeynmanIICh11Eq27",
    "FeynmanIICh11Eq28", "FeynmanICh10Eq7", "FeynmanIICh13Eq34",
    "FeynmanIICh13Eq23", "FeynmanICh48Eq2", "FeynmanICh41Eq16",
]
GOOD_TASKS = DANGEROUS_TASKS + SAFE_TASKS


def parse_linear_powerlaw(expr, feature_names):
    """expr: log(y) 关于 log(xᵢ) 的表达式。返回 (prefactor, exponents) 或 None。

    判定：对每个符号 xᵢ 的偏导必须是常数（不含符号），否则不是线性组合。
    """
    import sympy
    syms = sorted(expr.free_symbols, key=lambda s: str(s))
    mapping = {}
    for s in syms:
        m = re.fullmatch(r"x(\d+)", str(s))
        if m is None or int(m.group(1)) >= len(feature_names):
            return None
        mapping[s] = feature_names[int(m.group(1))]

    exponents = {}
    for s in syms:
        try:
            d = sympy.simplify(sympy.diff(expr, s))
        except Exception:
            return None
        if d.free_symbols:
            return None                      # 非线性项 → 非幂律
        v = float(d)
        if not np.isfinite(v):
            return None
        exponents[mapping[s]] = v
    for f in feature_names:                  # 未出现的变量指数为 0
        exponents.setdefault(f, 0.0)
    try:
        const = float(expr.subs({s: 0.0 for s in syms}))
    except Exception:
        return None
    if not np.isfinite(const):
        return None
    return float(np.exp(const)), exponents


def predict_powerlaw(prefactor, exponents, X, feature_names):
    out = np.full(X.shape[0], float(prefactor))
    for j, f in enumerate(feature_names):
        out = out * X[:, j] ** exponents[f]
    return out


def far_probe(eq_cls, X, feature_names, prefactor, exponents, m=12, mult=3.0):
    """一个数量级外的外推误差（相对真值中位数）。"""
    inst = eq_cls()
    Xf = np.tile(X.mean(axis=0), (m, 1))
    for j in range(Xf.shape[1]):
        Xf[:, j] = np.linspace(X[:, j].max(), X[:, j].max() * mult, m)
    with np.errstate(all="ignore"):
        yt = np.array([float(inst.eq_func(Xf[i])) for i in range(m)])
    yp = predict_powerlaw(prefactor, exponents, Xf, feature_names)
    ok = np.isfinite(yt) & np.isfinite(yp) & (yt > 0) & (yp > 0)
    if ok.sum() < 3:
        return np.nan
    return float(np.median(np.abs(yp[ok] - yt[ok]) / yt[ok]))


def run_one(task, eq_cls, n, seed, rel, niterations=25, eq_name=""):
    from pysr import PySRRegressor

    fn = srsd_data.make_dataset_fn(eq_cls)
    df, feats = fn(seed, n, rel)
    names = feats.split(",")
    X = df[names].to_numpy(float)
    y = df["target"].to_numpy(float)
    lx, ly = np.log(X), np.log(y)

    t0 = time.perf_counter()
    model = PySRRegressor(
        niterations=niterations, populations=12,
        binary_operators=["+", "-"], unary_operators=[],
        maxsize=14, deterministic=True, parallelism="serial",
        random_state=seed, temp_equation_file=False,
        verbosity=0, progress=False,
    )
    model.fit(lx, ly)
    dt = time.perf_counter() - t0

    expr = model.sympy()
    pred_log = model.predict(lx)
    rmse = float(np.sqrt(np.mean((ly - pred_log) ** 2)))
    ss_res = float(np.sum((ly - pred_log) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2_log = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    row = {"task": task, "eq_name": eq_name, "n": n, "seed": seed, "rel": rel,
           "pysr_seconds": dt, "pysr_expr": str(expr),
           "rmse_log": rmse, "r2_log": r2_log}

    parsed = parse_linear_powerlaw(expr, names)
    if parsed is None:
        row["parsed"] = False
        row["fail"] = "not_linear_in_logspace"
        return row
    prefactor, exponents = parsed
    row.update({"parsed": True, "prefactor": prefactor,
                "exponents": json.dumps(exponents, sort_keys=True)})

    # ---- 四个门卫，全部作用在 PySR 给出的候选上 ----
    fc = form_check_with_p(df, "target", names, seed)
    row["p_fc"] = fc.get("p_fc", np.nan)
    row["fail_fc"] = fc.get("fail_fc")

    b = reset_test(lx, ly, variant="fitted", vcov="classical")
    c = reset_test(lx, ly, variant="fitted", vcov="hc3")
    d = reset_test(lx, ly, variant="exog", vcov="hc3")
    row.update({"p_fitted_f": b.get("p", np.nan),
                "p_fitted_hc3": c.get("p", np.nan),
                "p_exog_hc3": d.get("p", np.nan),
                "fail_reset": "|".join(str(x.get("fail")) for x in (b, c, d)
                                       if x.get("fail"))})

    row["extrap_far_relerr"] = far_probe(eq_cls, X, names, prefactor,
                                         exponents)
    return row


def main(tasks=None, ns=(50,), seeds=(0,), rel=0.01, niterations=25):
    reg = srsd_data.registry_with_meta()
    by_key = {r["key"].lower(): r for r in reg}

    wanted = tasks or GOOD_TASKS
    found = []
    for t in wanted:
        r = by_key.get(t.lower())
        if r is not None:
            found.append(r)
        else:
            print(f'  [跳过] 未在 registry 中找到 {t}')
    print(f'任务数: {len(found)}')

    rows = []
    total = len(found) * len(ns) * len(seeds)
    k = 0
    for r in found:
        for n in ns:
            for seed in seeds:
                k += 1
                print(f'[{k}/{total}] {r["key"]} n={n} seed={seed} ...',
                      flush=True)
                try:
                    row = run_one(r["key"], r["cls"], n, seed, rel,
                                  niterations, eq_name=r.get("eq_name", ""))
                except Exception as exc:
                    row = {"task": r["key"], "eq_name": r.get("eq_name", ""),
                           "n": n, "seed": seed, "rel": rel,
                           "fail": f"{type(exc).__name__}: {exc}"[:120]}
                rows.append(row)
                if row.get("parsed"):
                    print(f'    rmse={row["rmse_log"]:.4f} '
                          f'R2={row["r2_log"]:.4f} '
                          f'p=({row.get("p_fc")},{row.get("p_fitted_f")},'
                          f'{row.get("p_fitted_hc3")},{row.get("p_exog_hc3")}) '
                          f'far={row.get("extrap_far_relerr")}', flush=True)
                else:
                    print(f'    !! {row.get("fail")}', flush=True)

    out = pd.DataFrame(rows)
    path = os.path.join(OUT, "pysr_gate_runs.csv")
    out.to_csv(path, index=False)
    print(f'\n已写入 {path}  ({len(out)} 行)')
    return out


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0,
                    help="只跑前 N 个任务（0=全部）")
    ap.add_argument("--n", type=int, default=50)
    ap.add_argument("--seeds", type=int, default=1,
                    help="种子数，从 0 开始")
    args = ap.parse_args()

    chosen = GOOD_TASKS[:args.limit] if args.limit else None
    main(tasks=chosen, ns=(args.n,), seeds=tuple(range(args.seeds)))