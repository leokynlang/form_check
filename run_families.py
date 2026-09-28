# -*- coding: utf-8 -*-
"""run_families.py — P2:Feynman 方程按函数族分层覆盖(逐方程报告)。

与 P1(run_srsd.py)的关系:P1 回答"32 个正值任务上的整体 FPR/功效";
P2 回答"哪类物理方程形式会被哪类检验漏检"。差异:
- 按 sympy 表达式结构分族(精确幂律/多项相加/指对数/三角/有理/多变量交互);
- 分层选择 20-30 个代表方程(正值可 log 者入检,非正值标 not-applicable);
- 检验配置沿用 P1 冻结值:powers=(2,3)、B=999 排列、α=0.05、rel∈{0.5%,1%,3%,5%};
- 新增异方差情形(SD 随 log y 线性 0.3→1.7×rel)与采样范围敏感性(官方域 ×0.5 / ×2);
- 逐方程结果表 + 按族汇总,绝不只报宏平均。

族判据只依赖表达式结构,选择先于任何检验结果。
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from srsd_data import registry_with_meta, make_dataset_fn    # noqa: E402
from common import run_all_methods, ALPHA                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# P1 冻结配置(P2 不得重调):
FROZEN = {"powers": (2, 3), "b_perm": 999, "alpha": 0.05,
          "rel_grid": [0.005, 0.01, 0.03, 0.05],
          "n_train": 200, "n_seeds": 30, "range_scales": [1.0, 0.5, 2.0]}


# ---------------------------------------------------------------- 族分类
def _monomial_exponents(f):
    """廉价单项式判据:顶层 Mul 分解,每个含变量的因子必须是 Symbol 或
    Symbol**常数。返回指数 dict 或 None(非单项式)。无需 simplify,O(因子数)。"""
    import sympy
    f = sympy.sympify(f)
    if isinstance(f, sympy.Add):
        return None
    factors = f.as_ordered_factors() if isinstance(f, sympy.Mul) else [f]
    exps = {}
    for fac in factors:
        if not fac.free_symbols:
            continue                       # 纯常数因子(exp(2)、π 等)
        if isinstance(fac, sympy.Symbol):
            exps[fac.name] = 1
        elif (isinstance(fac, sympy.Pow) and isinstance(fac.base, sympy.Symbol)
              and fac.exp.is_number):
            exps[fac.base.name] = fac.exp
        else:
            return None                    # exp(...)/log(...)/Add(1+x)**-1 等
    return exps if exps else None


def classify_family(sympy_eq):
    """按表达式结构分族。判据透明,写死,先于实验。"""
    import sympy
    f = sympy.sympify(sympy_eq)
    syms = list(f.free_symbols)
    if _monomial_exponents(f) is not None:
        return "exact_power_law"
    # 2) 结构扫描
    has = {"exp": False, "log": False, "trig": False, "add": False,
           "div": False, "sqrt": False, "tanh": False}
    for node in sympy.preorder_traversal(f):
        name = type(node).__name__
        if name == "exp":
            has["exp"] = True
        elif name == "log":
            has["log"] = True
        elif name in ("sin", "cos", "tan", "asin", "acos", "atan"):
            has["trig"] = True
        elif name == "tanh":
            has["tanh"] = True
        elif name == "Add" and node is not f:
            has["add"] = True
        elif name == "Pow" and (node.exp == -1 or
                                (hasattr(node.exp, "is_Rational") and
                                 node.exp.is_Rational and node.exp < 0)):
            has["div"] = True
        elif name == "Pow" and hasattr(node.exp, "is_Rational") and \
                node.exp.is_Rational and node.exp != 1 and node.exp != 2:
            has["sqrt"] = True
    if has["trig"] or has["tanh"]:
        return "trigonometric"
    if has["log"]:
        return "logarithmic"
    if has["exp"]:
        return "exponential"
    if has["div"] and has["add"]:
        return "rational"
    if has["add"] and len(syms) >= 3:
        return "multivar_additive"
    if has["add"]:
        return "additive"
    if has["sqrt"]:
        return "fractional_power"
    return "other"


def build_family_table():
    """全部 120 方程:名称/族/变量数/是否入检(P1 正值判据)/理由。"""
    rows = []
    for r in registry_with_meta():
        inst = r["cls"]()
        fam = classify_family(inst.sympy_eq)          # 用对象,别用截断的 srepr
        rows.append({"eq_name": r["eq_name"], "family": fam,
                     "klass": r["klass"], "n_vars": r["n_vars"],
                     "eligible": r["klass"] in ("H0", "H1"),
                     "reason": r["reason"], "cls": r["cls"]})
    return pd.DataFrame(rows)


def select_representatives(fam_table, per_family_cap=6, total_cap=30):
    """分层选择:每族至多 6 个。族内排序 = eligible 优先 → 变量数少 → 名称,
    避免 EXCLUDE 成员挤占同族可检方程的名额。总数超限则按族轮转裁剪。"""
    picked = []
    for fam, g in fam_table.groupby("family"):
        g2 = g.assign(_e=(~g["eligible"]).astype(int)).sort_values(
            ["_e", "n_vars", "eq_name"])
        cap = per_family_cap if g2.eligible.any() else 2
        picked.append(g2.head(cap))
    sel = pd.concat(picked)
    if len(sel) > total_cap:
        keep = []
        counts = {}
        for _, r in sel.iterrows():          # 族间轮转,保平衡
            counts[r["family"]] = counts.get(r["family"], 0) + 1
            if counts[r["family"]] <= max(2, total_cap // sel.family.nunique()):
                keep.append(r.name)
        sel = sel.loc[keep]
    return sel


def eval_cell(fn, rel, seed_base, seeds, hetero, range_scale):
    """一个 (方程,rel,noise,scale) 格 × seeds → 每方法拒绝/失败/时间。"""
    rows = []
    for s in seeds:
        bs = seed_base + s * 97
        row = {"seed": s}
        try:
            df, feat_str = fn(seed=bs, n=FROZEN["n_train"], rel=rel,
                              hetero=hetero, range_scale=range_scale)
            feats = feat_str.split(",")
        except ValueError as e:
            row["data_fail"] = str(e)
            rows.append(row)
            continue
        row["data_fail"] = None
        res = run_all_methods(df, "target", feats, seed=bs,
                              n_boot=FROZEN["b_perm"])
        for k in ("rej_fc", "rej_fitted_f", "rej_fitted_hc3", "rej_exog_hc3",
                  "p_fc", "p_fitted_f", "p_fitted_hc3", "p_exog_hc3",
                  "fail_fc", "fail_fitted_f", "fail_fitted_hc3", "fail_exog_hc3"):
            if k in res:
                row[k] = res[k]
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="feynman_families")
    ap.add_argument("--seeds", type=int, default=FROZEN["n_seeds"])
    args = ap.parse_args()

    run_dir = os.path.join(HERE, "runs", args.tag)
    os.makedirs(run_dir, exist_ok=True)
    fam_table = build_family_table()
    fam_table.drop(columns=["cls"]).to_csv(
        os.path.join(run_dir, "family_table_all120.csv"), index=False)
    sel = select_representatives(fam_table)
    sel.drop(columns=["cls"]).to_csv(
        os.path.join(run_dir, "selection_table.csv"), index=False)

    seeds = list(range(args.seeds))
    grid = [(r["eq_name"], rel, het, sc)
            for _, r in sel.iterrows()
            if r["eligible"] and r["klass"] == "H1"     # 有曲率可检的对象
            for rel in FROZEN["rel_grid"]
            for het in (False, True)
            for sc in FROZEN["range_scales"]]
    # H0 族(单项式)只跑官方域(它们的角色是各族的 FPR 校准,不测 scale 交互)
    grid_h0 = [(r["eq_name"], rel, het, 1.0)
               for _, r in sel.iterrows()
               if r["klass"] == "H0"
               for rel in FROZEN["rel_grid"]
               for het in (False, True)]

    cfg = {"frozen": {k: list(v) if isinstance(v, tuple) else v
                      for k, v in FROZEN.items()},
           "alpha": ALPHA,
           "note": "配置沿用 P1 冻结值;族分类只用于分层选择,不改变 H0/H1 判据",
           "started": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1)

    rows = []
    t0 = time.perf_counter()
    total = len(grid) + len(grid_h0)
    done = 0
    buf = []
    for eq, rel, het, sc in grid + grid_h0:
        done += 1
        rec = sel[sel.eq_name == eq].iloc[0]
        fn = make_dataset_fn(rec["cls"])
        cell = eval_cell(fn, rel, seed_base=300_000 + done * 1_000_000,
                         seeds=seeds, hetero=het, range_scale=sc)
        cell.insert(0, "range_scale", sc)
        cell.insert(0, "hetero", het)
        cell.insert(0, "rel", rel)
        cell.insert(0, "eq_name", eq)
        cell.insert(0, "family", rec["family"])
        cell.insert(0, "klass", rec["klass"])
        rows.append(cell)
        buf.append(cell)
        na = int(cell.data_fail.notna().sum()) if "data_fail" in cell else 0
        print(f"[{done:3d}/{total}] {eq:22s} {rec['family']:16s} "
              f"rel={rel:.3f} het={int(het)} sc={sc:.1f} "
              f"na={na}/{len(cell)} ({time.perf_counter()-t0:6.0f}s)", flush=True)
        if len(buf) >= 8:
            _append(run_dir, buf)
            buf = []
    if buf:
        _append(run_dir, buf)
    cfg["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    cfg["wall_seconds"] = time.perf_counter() - t0
    with open(os.path.join(run_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, ensure_ascii=False, indent=1)
    print(f"DONE {args.tag}: {total} cells -> {run_dir}")


def _append(run_dir, buf):
    new = pd.concat(buf, ignore_index=True)
    p = os.path.join(run_dir, "raw_families.parquet")
    if os.path.exists(p):
        old = pd.read_parquet(p)
        new = pd.concat([old, new], ignore_index=True)
    new.to_parquet(p, index=False)
    new.to_csv(os.path.join(run_dir, "raw_families.csv"), index=False)


if __name__ == "__main__":
    main()
