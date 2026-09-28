# -*- coding: utf-8 -*-
"""battery_probe.py — 探针:多函数族"形式合格性电池"是否可行(核心代码零改动)。

[动机] core form_check 只检验一种零形式(log-log 直线=幂律),报警时只能说
"幂律可能不适",说不出"那什么才合适"。P2 实锤各检测器有族级盲区
(fc 漏 exponential 0.11;fitted-F 漏 rational 0.13)。本探针验证:
同一套"线性化坐标 + 二次曲率列 + 残差排列零分布"机制,对每个可线性化
函数族各跑一遍,能否 (a) 保持幂律 H0 校准 (b) 正确路由到真函数族。

[语义转变(与 P1/P2 的关键区别)] 电池里每族的检验是独立的"形式合格性":
family F 的检验 H0 = "真值条件均值在 F 的线性化坐标里是线性的"。
可行集 viable = {F: p_F >= alpha}(非拒绝),输出候选列表而非单一赢家
——与项目宪法"输出证据与候选列表"一致。每族各自 alpha=0.05,
不做 max-T 剧场(非拒绝语义下多重性不膨胀第一类错误;
max-T 留给"选最小 p 当路由"的用法,那是 promotion 阶段的决策)。

[预注册预测(写于运行前,用于对答案)]
P1. H0_BASELINE(幂律真值): rej(power)∈[0.02,0.10](校准);
    其余三族 rej≥0.9(幂律数据在 exp/log/rat 坐标下确有曲率);
    containment(真族∈viable)≥0.95;中位 |viable| ≤ 1.2。
P2. EXP_MODEL: rej(power)=1.00;viable(exponential)≥0.90;
    containment≥0.95。
P3. LOG_MODEL: viable(logarithmic)≥0.85——这是最有压力的格:
    log 族在 raw 空间拟合,乘性噪声在 raw 空间是异方差的,
    排列零分布的可交换性假设受压;此格测的就是这个。
P4. RATIONAL: viable(rational)≥0.90(精确线性化 1/y=0.1+0.2/x);
    rej(power)≈1.0(饱和曲线在 log-log 强凹)。
P5. THRESHOLD(分段,尖锐变点): 四族全拒 → 升级路径触发率≥0.90
    (无单一线性化形式能吸收折角 → "无可行族→升级 piecewise/auto")。

[边界(诚实留白)] P2 报告的 exponential 盲区(Feynman 真实方程,
T≫采样量程时 log-log 本来就近直线)是域内不可约——任何域内检验都低功效,
电池不解决它;那格属于量纲层(x/T 量纲一化)与外推告警的分工。
多特征扩展、additive-noise H0 下 raw 空间族的 FPR 校准,留正式版。

[与 fc 的关系] power 族的排列 p 值 = form_check 本体统计量的 p 值化
(同一 c2,同一零分布构造)。首格做一致性 sanity(|Δt|<30%,n_boot 不同)。

配置:沿用 P1/P2 冻结纪律 α=0.05, B=999;n=200, seeds=200,
rel∈{0.01,0.03};种子纪律 base=700000+seed*97+rel_idx*13+fam_idx*1000003。
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

# 族注册:线性化坐标 u(x) 与响应变换 v(y)。判据透明,先于实验。
FAMILIES = {
    "power":       (lambda x: np.log(x), lambda y: np.log(y)),
    "exponential": (lambda x: x,           lambda y: np.log(y)),
    "logarithmic": (lambda x: np.log(x), lambda y: y),
    "rational":    (lambda x: 1.0 / x,   lambda y: 1.0 / y),
}
FAM_IDX = {f: i for i, f in enumerate(FAMILIES)}

# 预注册:场景 -> 真函数族(THRESHOLD=无单族,应触发升级)
TRUTH_FAMILY = {
    "H0_BASELINE": "power",
    "EXP_MODEL": "exponential",
    "LOG_MODEL": "logarithmic",
    "RATIONAL": "rational",
    "THRESHOLD": None,
}
SCENARIO_LIST = list(TRUTH_FAMILY)


def spec_test(x, y, uf, vf, seed):
    """单族形式合格性检验:线性化坐标下 二次曲率系数的残差排列 p 值。

    与 form_check 本体同构:base=[1,u],design=[1,u,u²](u 标准化定标),
    统计量=design 最后一列系数;零分布=线性拟合残差排列重建伪数据集。
    c2 是 v 的线性泛函 → pinv 行技巧向量化 B=999(与 common.py 同法)。
    """
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
    idx = np.argsort(rng.random((B_PERM, len(v))), axis=1)
    null = c2_base + resid[idx] @ pinv_row
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
                ps = {}
                for fam, (uf, vf) in FAMILIES.items():
                    bs = 700_000 + s * 97 + rel_i * 13 + FAM_IDX[fam] * 1_000_003
                    r = spec_test(x, y, uf, vf, bs)
                    rec[f"p_{fam}"] = r["p"]
                    rec[f"rej_{fam}"] = bool(r["p"] < ALPHA)
                    rec[f"fail_{fam}"] = r["fail"]
                    ps[fam] = r["p"]
                viable = {f for f, p in ps.items() if np.isfinite(p) and p >= ALPHA}
                rec["viable"] = ",".join(sorted(viable))
                rec["n_viable"] = len(viable)
                truth = TRUTH_FAMILY[scen]
                rec["containment"] = (truth in viable) if truth else None
                rec["escalation_fired"] = (len(viable) == 0)
                rows.append(rec)
                if done % 500 == 0:
                    print(f"[{done}/{total}] {scen} rel={rel} "
                          f"({time.perf_counter()-t0:.0f}s)", flush=True)

    raw = pd.DataFrame(rows)
    out_dir = os.path.join(HERE, "runs", "battery_probe")
    os.makedirs(out_dir, exist_ok=True)
    raw.to_parquet(os.path.join(out_dir, "raw.parquet"), index=False)

    print("\n=== 拒绝率表(rej rate @α=0.05,200 种子)===")
    print("场景 | rel | " + " | ".join(f"rej_{f}" for f in FAMILIES)
          + " | 中位n_viable | containment | escalation")
    summary = {}
    for scen in SCENARIO_LIST:
        for rel in RELS:
            g = raw[(raw.scenario == scen) & (raw.rel == rel)]
            rates = {f: float(g[f"rej_{f}"].mean()) for f in FAMILIES}
            med_v = float(g.n_viable.median())
            cont = (float(g.continuity.astype(float).mean())
                    if False else
                    (float(g["containment"].dropna().astype(float).mean())
                     if g["containment"].notna().any() else None))
            esc = float(g.escalation_fired.mean())
            summary[f"{scen}@rel={rel}"] = {
                "rej_rates": rates, "median_n_viable": med_v,
                "containment": cont, "escalation": esc}
            print(f"{scen:12s} | {rel:.2f} | "
                  + " | ".join(f"{rates[f]:.2f}" for f in FAMILIES)
                  + f" | {med_v:.1f} | "
                  + (f"{cont:.2f}" if cont is not None else "n/a")
                  + f" | {esc:.2f}")

    # viable 组合分布(H0 上"只有 power 可行"的纯度)
    print("\n=== viable 组合分布(scen@rel) ===")
    combos = {}
    for scen in SCENARIO_LIST:
        for rel in RELS:
            g = raw[(raw.scenario == scen) & (raw.rel == rel)]
            top = g.viable.value_counts(normalize=True)
            combos[f"{scen}@{rel}"] = {k: round(float(v), 3)
                                       for k, v in top.head(4).items()}
            print(f"{scen}@{rel}: " + str(combos[f"{scen}@{rel}"]))

    # fc 一致性 sanity(首格):power 族 t vs dimreflex form_check t
    try:
        sys.path.insert(0, os.path.dirname(HERE))       # 桌面根(包父目录)
        from dimreflex.discover import form_check as fc_native
        df, feat_str = SCENARIOS["H0_BASELINE"](seed=700_000, n=N, rel=0.01)
        nat = fc_native(df, target="target", features=["x1"],
                        n_boot=200, random_state=700_000)
        mine = spec_test(df["x1"].values, df["target"].values,
                         *FAMILIES["power"], seed=700_000 + FAM_IDX["power"] * 1_000_003)
        ratio = abs(mine["t"] - nat["t_stat"]) / max(abs(nat["t_stat"]), 1e-9)
        print(f"\nfc sanity: native t={nat['t_stat']:.3f} "
              f"probe t={mine['t']:.3f} rel_diff={ratio:.1%} "
              f"(n_boot 200 vs 999,预期同号且同量级)")
        summary["fc_sanity"] = {"native_t": nat["t_stat"], "probe_t": mine["t"]}
    except Exception as e:                               # noqa: BLE001
        print(f"\nfc sanity 跳过: {e}")

    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump({"summary": summary, "combos": combos,
                   "config": {"B_PERM": B_PERM, "ALPHA": ALPHA, "N": N,
                              "SEEDS": SEEDS, "RELS": list(RELS)}},
                  f, ensure_ascii=False, indent=1)
    print(f"\nDONE {time.perf_counter()-t0:.0f}s -> {out_dir}")


if __name__ == "__main__":
    main()
