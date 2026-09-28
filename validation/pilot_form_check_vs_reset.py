"""pilot_form_check_vs_reset.py — 主张1成败gate的1小时pilot

[PAPER_PLANS 预注册 gate]:
  噪声 1%、n>=200 时, form_check 的 TPR@FPR<=5% >= 0.8,
  真幂律误报率 <= 0.05, 且功效/误报不劣于 RESET(教科书实现)。
  全面劣于 RESET -> 删检测器 claim, 论文降级为现象学 benchmark。

[设计] 复用 validate_friction 的实验台(同物理、同 pi 群):
  H1: truth=prandtl_karman (形式错, 应报警) x 50 种子
  H0: truth=power_law      (形式对, 不应误报) x 50 种子
  噪声 1% 乘性, n=200, Re 范围 [4e3,1e7]。
  两检测器在完全相同的 log-log (lx, ly) 数据上跑:
    form_check: dimreflex 偏残差曲率 + 残差bootstrap零分布, |t|>3
    RESET: Ramsey(1969) 教科书实现 — 线性拟合后把 yhat^2, yhat^3
           加入回归, 对两个新增系数做 F 检验, p<0.05 报警
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import os as _os
_here = _os.path.dirname(_os.path.abspath(__file__))  # .../form_check/validation
_root = _os.path.dirname(_here)
sys.path.insert(1, _here)     # validate_friction 同目录(位 1,让位给 _bootstrap)
sys.path.insert(0, _root)
import _bootstrap  # noqa: E402 — vendored dimreflex 强制优先
import dimreflex as dr
from validate_friction import make_experiment

NOISE = 0.01
N = 200
N_SEEDS = 50


def pi_columns(df):
    """框架自己的 pi 群(与 validate_friction.shape_probes 同源):
    pi0 = dp/(rho v^2), pi1 = rho v D / mu(=Re)。"""
    pi0 = (df["dp"] / (df["rho"] * df["v"] ** 2)).values
    pi1 = (df["rho"] * df["v"] * df["D"] / df["mu"]).values
    return pd.DataFrame({"pi0": pi0, "pi1": pi1})


def reset_test(lx, ly, alpha=0.05):
    """Ramsey RESET 教科书实现: ly ~ 1 + lx 受限, 加 yhat^2/yhat^3 重估,
    F 检验新增系数联合为零。返回 (alarm, p_value, F)。"""
    n = len(ly)
    X_r = np.column_stack([np.ones(n), lx])
    beta_r = np.linalg.lstsq(X_r, ly, rcond=None)[0]
    yhat = X_r @ beta_r
    X_u = np.column_stack([X_r, yhat ** 2, yhat ** 3])
    rss_r = float(np.sum((ly - yhat) ** 2))
    beta_u = np.linalg.lstsq(X_u, ly, rcond=None)[0]
    rss_u = float(np.sum((ly - X_u @ beta_u) ** 2))
    q, k = 2, X_u.shape[1]                 # 新增约束数 q=2, 非约束参数数 k
    F = ((rss_r - rss_u) / q) / (rss_u / (n - k))
    p = float(stats.f.sf(F, q, n - k))
    return bool(p < alpha), p, float(F)


def run_one(seed, truth):
    df = make_experiment(seed=90000 + seed, noise=NOISE, n=N,
                         re_hi=1e7, truth=truth)
    pi = pi_columns(df)
    fc = dr.form_check(pi, target="pi0", features=["pi1"], n_boot=200,
                       random_state=seed)
    lx = np.log(pi["pi1"].values)
    ly = np.log(pi["pi0"].values)
    alarm_r, p_r, F_r = reset_test(lx, ly)
    return {"seed": seed, "truth": truth,
            "fc_alarm": bool(fc["alarm"]), "fc_t": fc["t_stat"],
            "reset_alarm": alarm_r, "reset_p": p_r, "reset_F": F_r}


def main():
    rows = []
    for truth in ("prandtl_karman", "power_law"):
        for s in range(N_SEEDS):
            rows.append(run_one(s, truth))
            if (s + 1) % 10 == 0:
                print(f"  {truth}: {s + 1}/{N_SEEDS}", flush=True)
    res = pd.DataFrame(rows)

    pk = res[res.truth == "prandtl_karman"]
    pl = res[res.truth == "power_law"]
    fc_tpr = float(pk.fc_alarm.mean())
    fc_fpr = float(pl.fc_alarm.mean())
    rs_tpr = float(pk.reset_alarm.mean())
    rs_fpr = float(pl.reset_alarm.mean())

    print()
    print(f"===== pilot 结果 (noise=1%, n={N}, {N_SEEDS} seeds/truth) =====")
    print(f"form_check  TPR={fc_tpr:.2f}  FPR(真幂律)={fc_fpr:.2f}  "
          f"mean|t| PK={pk.fc_t.abs().mean():.1f}  PL={pl.fc_t.abs().mean():.1f}")
    print(f"RESET       TPR={rs_tpr:.2f}  FPR(真幂律)={rs_fpr:.2f}  "
          f"median p PK={pk.reset_p.median():.2e}  PL={pl.reset_p.median():.3f}")

    gate_tpr = fc_tpr >= 0.8
    gate_fpr = fc_fpr <= 0.05
    not_worse = (fc_tpr >= rs_tpr - 1e-9) or (fc_fpr <= rs_fpr - 1e-9)
    print()
    print("===== 预注册 gate 判定 =====")
    print(f"  TPR>=0.8     : {fc_tpr:.2f} -> {'PASS' if gate_tpr else 'FAIL'}")
    print(f"  误报<=0.05   : {fc_fpr:.2f} -> {'PASS' if gate_fpr else 'FAIL'}")
    print(f"  不劣于RESET  : TPR {fc_tpr:.2f} vs {rs_tpr:.2f}, "
          f"FPR {fc_fpr:.2f} vs {rs_fpr:.2f} -> "
          f"{'PASS' if not_worse else 'FAIL'}")
    overall = gate_tpr and gate_fpr and not_worse
    print("  总判定: %s" % ("GATE PASS (可冲 ML4PS 检测器主张)" if overall
                           else "GATE FAIL (按预注册降级/补实验)"))

    out_csv = Path(__file__).resolve().parent / "pilot_form_check_vs_reset_results.csv"
    res.to_csv(out_csv, index=False)
    print(f"明细已存 {out_csv}")
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()
