"""validate_friction.py — 负对照:光滑圆管摩擦系数(Prandtl-Karman 对数律)

[背景] 框架的默认假设是幂律。正对照(泵相似律)证明"真值是幂律时能恢复";
本脚本证明相反的方向:真值**不是**幂律时,框架是否诚实地暴露形式失配。
只会说"拟合得好"的框架没有可信度;能说"这个形式不合适"才是卖点。

[真值] 1/sqrt(f) = 2*log10(Re*sqrt(f)) - 0.8
  Prandtl-Karman 光滑管律(Nikuradse 实验,教科书圣杯)。
  关键性质:隐式(f 两边出现)、非幂律(f ~ 1/log^2(Re),比任何幂都衰减得慢),
  因此局部指数随 Re 持续漂移——任何有限幂律都覆盖不了的指纹。

[数据] 合成实验台(seed=20260901,可复现):光滑圆管 D=0.1m, L/D=100,
  水为工质, T∈[20,60]C(给 rho/mu 真实散布), Re∈[4e3,1e7] 60 点,
  dp 加 0.5% 乘性噪声(取对数后=加性高斯,恰落在框架 log 空间拟合的
  舒适区——公平对抗:只在"形式"这一维度设卡)。

[A] 预注册判据(先写判据后跑;初稿 2026-09-01,同日噪声带标定后修订):
  C1 形式失配有证据(观测 = 有噪管道自己的 pi 群):
     滑窗指数跨度 > 本设计 H0 噪声带 q95, 且 |曲率 t| > 3
     [标定] 长量程(4e3..1e7) H0:跨度 q95=0.017,曲率 sigma=0.00032
     [观测] 跨度 0.127(q95 的 7.4 倍),曲率 t=+21   -> 预期报警
  C2 幂律路径不静默:method='power_law' 须对形式失配输出警告。
     [2026-09-01 首版] 实测静默(R²=0.991/Level 3/零警告) = 缺口实锤,
     当时预期 FAIL;同日 form_check 进入 discover.py(偏残差曲率 +
     残差 bootstrap 零分布)后,本判据按预注册翻转为:必须留痕。
     (审计轨迹:预注册预期随框架版本翻转,翻转本身记录在案)
  C3 候选排序:auto 路径 rational_log rank 1 且 recommended,
     cv_r_squared >= 0.99(标定值 0.9999);power_law 排其后
  C4 诚实失败:symbolic(PySR 未装)success=False 且 reason 说明,
     discover 不抛异常
  C5 特异性(防狼来了):Blasius 区(4e3..1e5,真值 0.3164*Re^-0.25 确为
     幂律)用**该设计自己的** H0 噪声带判定,不得报警
     [标定] 短量程 H0:跨度 q95=0.043,曲率 sigma=0.00187
     [教训,已内化] 阈值是设计条件的函数,不可跨设计移植:
     曲率 sigma 反比于 x 量程跨度的平方(5.9x = (3.4/1.4)^2,实测)。
     框架启示:form_check 必须"对本数据网格+本噪声水平"逐例标定,
     写死常数阈值(如 anomaly 的 条件数>30)同属此脆弱类别。
  C6 外推考试(记录,不判):两候选外推 Re=1e8,对比真值
     (幂律外推的死亡方式是本脚本的展示品,不是判据)

用法: python validate_friction.py   # exit 0 = 全部与预注册一致, 1 = 偏离
"""

import contextlib
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq

SEED = 20260901
NOISE = 0.005            # 0.5% 乘性测量噪声
N = 60
D_PIPE = 0.1             # 管径 m
L_OVER_D = 100.0


# ---------------------------------------------------------------- 真值与实验台
def prandtl_karman_f(re):
    """给定雷诺数,数值求解光滑管摩擦系数 f。

    f 在方程两边 → 隐式方程,移不出 f = 表达式。
    brentq 在给定区间内找 g(f)=0 的根:二分法家族,每猜一次中点,
    按"高了/低了"把区间砍半,给了区间必收敛。
    """
    def g(f):
        return 1.0 / np.sqrt(f) - 2.0 * np.log10(re * np.sqrt(f)) + 0.8
    return brentq(g, 1e-6, 0.1)


def _fluid(seed, n):
    """水, 20-60C:给 rho/mu 真实散布(否则是常量列,误触框架已知短板,
    负对照就测不到想测的东西了——对抗数据的原则:只在目标维度设卡)。"""
    rng = np.random.default_rng(seed)
    T = rng.uniform(20.0, 60.0, n)
    rho = 997.0 * (1.0 - 2.1e-4 * (T - 20.0))      # 线性近似,够用
    mu = 1.002e-3 * 10.0 ** ((20.0 - T) / 12.0)    # 指数近似,量级正确
    return rho, mu


def make_experiment(seed=SEED, noise=NOISE, n=N, re_hi=1e7,
                    truth="prandtl_karman"):
    """生成一次"实验",返回 DataFrame(列名 = dimreflex 变量名)。

    default_rng(seed) 固定随机种子 = 考卷印死:任何人重跑得到完全相同的
    数据,判据钉在固定数据上,不是每次重抽。
    truth: "prandtl_karman"(对数律真值) | "power_law"(H0 对照与 Blasius 区)
    """
    rng = np.random.default_rng(seed + 1)          # +1: 噪声流与 _fluid 独立
    rho, mu = _fluid(seed, n)
    Re = np.logspace(np.log10(4e3), np.log10(re_hi), n)
    v = Re * mu / (rho * D_PIPE)                   # 由 Re 反解流速
    if truth == "prandtl_karman":
        # brentq 只吃标量 → 列表推导逐个求根(不是所有函数都能向量化)
        f = np.array([prandtl_karman_f(r) for r in Re])
    else:                                          # H0: 纯幂律(Blasius 系数)
        f = 0.3164 * Re ** -0.25
    dp = f * (rho * v ** 2 / 2.0) * L_OVER_D       # Darcy-Weisbach 压降
    # 乘性噪声:仪器标定是"读数的±0.5%",误差跟着读数走
    dp = dp * (1.0 + noise * rng.standard_normal(n))
    return pd.DataFrame({"dp": dp, "v": v,
                         "D": np.full(n, D_PIPE),  # np.full:造常数列
                         "rho": rho, "mu": mu})


# ---------------------------------------------------------------- 形状体检
def shape_probes(df, win=12):
    """对框架自己的 pi 群做两项"形状体检"——判据 C1/C5 的观测值。

    体检对象 = PI_0 = dp/(rho*v^2) 与 PI_1 = rho*v*D/mu(=Re),
    正是 discover() 输出的那两个 pi 群(差常数倍不影响形状)。
    判据建立在框架的表示上,框架无法辩解"你量错了东西"。

    log-log 空间里幂律 = 直线。两把独立的尺子量"弯不弯":
      二次曲率:拟合 y = c2*x^2 + b*x + a,纯幂律的 c2 = 0(全局视角)
      滑窗指数:每 win 点拟合局部斜率,真幂律的斜率处处相同(局部视角)
    """
    pi0 = df["dp"] / (df["rho"] * df["v"] ** 2)
    pi1 = df["rho"] * df["v"] * df["D"] / df["mu"]
    x = np.log10(pi1.values)
    y = np.log10(pi0.values)
    # polyfit 返回系数从高次到低次排列(经典坑:第一次几乎人人栽在这)
    c2, _, _ = np.polyfit(x, y, 2)
    local = np.array([np.polyfit(x[i:i + win], y[i:i + win], 1)[0]
                      for i in range(len(x) - win + 1)])
    return {"curvature": float(c2),
            "slope_span": float(local.max() - local.min()),
            "local_slopes": local}


def noise_band(re_hi, truth_exp="power", n_sims=200):
    """本设计的 H0 噪声带:真值换成纯幂律、其余一切相同,模拟 n_sims 次。

    判据阈值不是拍脑袋的常数:在"纯幂律+同噪声"的平行宇宙里看体检
    统计量能虚报到多少,取分位数。零分布逐设计标定(C5 的教训)。
    """
    spans, curvs = [], []
    for s in range(n_sims):
        p = shape_probes(make_experiment(seed=5000 + s, re_hi=re_hi,
                                         truth=truth_exp))
        spans.append(p["slope_span"])
        curvs.append(p["curvature"])
    spans, curvs = np.array(spans), np.array(curvs)
    return {"span_q95": float(np.quantile(spans, 0.95)),
            "curv_mean": float(curvs.mean()),
            "curv_std": float(curvs.std())}


def judge(probe, band):
    """观测值对照本设计的噪声带。返回 (是否报警, 曲率 t 值)。"""
    t_curv = ((probe["curvature"] - band["curv_mean"])
              / band["curv_std"])
    alarm = probe["slope_span"] > band["span_q95"] or abs(t_curv) > 3.0
    return bool(alarm), float(t_curv)


# ---------------------------------------------------------------- 框架调用
def run_discover(data, method):
    """调 dimreflex(2026-09-05 拆分:vendored 包在本文件夹 ../dimreflex/)。"""
    import os as _os
    _root = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
    sys.path.insert(1, _root)   # 位 1:让位给上面已插入的本目录位 0
    import _bootstrap  # noqa: E402 — vendored dimreflex 强制优先
    import dimreflex as dr
    s = dr.System("smooth_pipe")
    s.add_variable("dp", "Pa", "target")
    s.add_variable("v", "m/s", "state")
    s.add_variable("D", "m", "geometry")
    s.add_variable("rho", "kg/m**3", "material_property")
    s.add_variable("mu", "Pa*s", "material_property")
    return dr.discover(data, s, target="dp", method=method)


FORM_WORDS = ("形式", "曲率", "幂律偏差", "弯曲", "form", "curvat", "shape")


def check_c2(method="power_law"):
    """C2: 幂律路径对形式失配必须留痕(form_check, v0.5.0 起)。

    首版(2026-09-01 上午)本函数检查"是否静默"并预期 FAIL;
    form_check 进入框架后,按预注册翻转为"必须报警"。"""
    res = run_discover(make_experiment(), method)
    fit = res.get("fit") or {}
    warnings = (res.get("data_report") or {}).get("warnings") or []
    form_warnings = [w for w in warnings
                     if any(k in w for k in FORM_WORDS)]
    flagged = bool(form_warnings)
    fc = res.get("form_check") or {}
    desc = (f"formula={fit.get('formula')!r}  R2={fit.get('r_squared'):.4f}"
            f"  form_check.alarm={fc.get('alarm')}"
            f"  form_check.t={fc.get('t_stat', 0):+.1f}"
            f"  形式相关警告数={len(form_warnings)}")
    return flagged, desc


def check_c3_c4():
    """C3+C4: auto 路径的候选排序与诚实失败。"""
    res = run_discover(make_experiment(), "auto")
    auto = res.get("autonomous") or {}
    assessments = {a["family"]: a for a in auto.get("assessments", [])}

    rat = assessments.get("rational_log", {})
    pw = assessments.get("power_law", {})
    c3_ok = bool(
        rat.get("success") and rat.get("recommended")
        and rat.get("rank") == 1
        and (rat.get("validation") or {}).get("cv_r_squared", 0) >= 0.99
        and pw.get("rank", 99) > (rat.get("rank") or 0))
    c3_desc = (f"rational_log rank={rat.get('rank')} recommended="
               f"{rat.get('recommended')} cv_r2="
               f"{(rat.get('validation') or {}).get('cv_r_squared'):.5f}; "
               f"power_law rank={pw.get('rank')} cv_r2="
               f"{(pw.get('validation') or {}).get('cv_r_squared'):.5f}")

    sym = assessments.get("symbolic", {})
    c4_ok = bool(sym and not sym.get("success")
                 and "PySR" in (sym.get("reason") or ""))
    c4_desc = f"symbolic success={sym.get('success')} reason={sym.get('reason')!r}"
    return c3_ok, c3_desc, c4_ok, c4_desc, res


def extrapolation_exam(auto_res):
    """C6: 用框架自己拟合的参数外推到 Re=1e8,和真值对答案(只记录)。"""
    assessments = {a["family"]: a
                   for a in (auto_res.get("autonomous") or {}).get("assessments", [])}
    re_star = 1e8
    f_true = prandtl_karman_f(re_star)
    pi0_true = f_true * L_OVER_D / 2.0          # PI_0 = dp/(rho v^2) = f/2*L/D

    pw = assessments["power_law"]["model"]["params"]
    b = pw["exponents"]["PI_1_v^1_D^1_rho^1_mu^-1"]
    pi0_pw = pw["prefactor"] * re_star ** b

    rl = assessments["rational_log"]["model"]["params"]
    c0, c1 = rl["num_coefs"]                    # 分子: c0 + c1*L
    d0 = rl["den_coefs"][0]                     # 分母: 1 + d0*L
    L = np.log10(re_star) - rl["log_center"]    # 居中的 log10(PI_1)
    pi0_rl = 10.0 ** ((c0 + c1 * L) / (1.0 + d0 * L))

    print(f"\n=== C6 外推考试(记录,不判): Re = {re_star:.0e}, "
          f"真值 PI_0 = {pi0_true:.4f} ===")
    print(f"  power_law   PI_0 = {pi0_pw:.4f}   误差 "
          f"{(pi0_pw / pi0_true - 1) * 100:+.1f}%")
    print(f"  rational_log PI_0 = {pi0_rl:.4f}   误差 "
          f"{(pi0_rl / pi0_true - 1) * 100:+.1f}%")
    print("  (指数 b 的 std_error 极小 —— 统计上'精确'与形式上'正确'是两回事)")


# ---------------------------------------------------------------- 主判据
def main():
    print("=== C1 形式失配证据(PK 真值, 全量程 4e3..1e7) ===")
    data = make_experiment()
    probe = shape_probes(data)
    band_long = noise_band(1e7)
    c1_alarm, t_curv = judge(probe, band_long)
    print(f"观测: 跨度 {probe['slope_span']:.3f}, 曲率 {probe['curvature']:+.4f} "
          f"(t={t_curv:+.1f})")
    print(f"H0 带: 跨度 q95 {band_long['span_q95']:.3f}, "
          f"曲率 3σ {3 * band_long['curv_std']:.4f}")
    print(f"滑窗指数首→末: {probe['local_slopes'][0]:+.3f} → "
          f"{probe['local_slopes'][-1]:+.3f}  (真幂律应恒定)")
    print(f"C1 {'PASS: 报警' if c1_alarm else 'FAIL: 未报警'}")

    print("\n=== C2 幂律路径不静默(v0.5.0 起:必须对形式失配留痕) ===")
    flagged, desc = check_c2()
    print(f"  {desc}")
    c2_msg = ("PASS: 框架已对形式失配留痕(缺口已闭合)" if flagged
              else "FAIL: 框架对形式失配零警告(静默缺口复发?)")
    print(f"C2 {c2_msg}")

    print("\n=== C3/C4 auto 路径: 候选排序与诚实失败 ===")
    c3_ok, c3_desc, c4_ok, c4_desc, auto_res = check_c3_c4()
    print(f"  {c3_desc}\n  C3 {'PASS' if c3_ok else 'FAIL'}")
    print(f"  {c4_desc}\n  C4 {'PASS' if c4_ok else 'FAIL'}")

    print("\n=== C5 特异性(Blasius 区, 真值确为幂律, 不得报警) ===")
    data_b = make_experiment(re_hi=1e5, truth="power_law")
    probe_b = shape_probes(data_b)
    band_short = noise_band(1e5)
    alarm_b, t_b = judge(probe_b, band_short)
    print(f"观测: 跨度 {probe_b['slope_span']:.3f}, 曲率 "
          f"{probe_b['curvature']:+.4f} (t={t_b:+.1f})")
    print(f"H0 带(本设计): 跨度 q95 {band_short['span_q95']:.3f}, "
          f"曲率 3σ {3 * band_short['curv_std']:.4f}")
    print(f"C5 {'FAIL: 误报' if alarm_b else 'PASS: 不报警'}")

    # C5b(记录,不判): 同一短窗口装回 PK 真值 —— 检查应当报警,
    # 因为 PK 即使在 Blasius 窗口内也真实地不是幂律(指数仍漂移)。
    # 教训存档: 首跑曾在此处忘传 truth="power_law", 检查对 PK 短窗
    # 报警(t=+15), exit 1 判"偏离预注册" —— 审计机制审计了审计者。
    probe_pk = shape_probes(make_experiment(re_hi=1e5))
    alarm_pk, t_pk = judge(probe_pk, band_short)
    print(f"C5b(记录): PK 短窗 跨度 {probe_pk['slope_span']:.3f}, "
          f"t={t_pk:+.1f} -> {'报警(真信号, 敏感性余量大)' if alarm_pk else '未报警'}")

    extrapolation_exam(auto_res)

    verdict = c1_alarm and flagged and c3_ok and c4_ok and not alarm_b
    print("\n=== 总判 ===")
    print("与预注册一致(C1 报警, C2 留痕, C3/C4 过, C5 不误报):",
          "YES" if verdict else "NO -> 审计偏离, 查预注册 [A]")
    return 0 if verdict else 1


if __name__ == "__main__":
    # 控制台编码重配置失败不阻断验证(GBK 控制台仅影响显示,不影响判据)
    with contextlib.suppress(Exception):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.exit(main())
