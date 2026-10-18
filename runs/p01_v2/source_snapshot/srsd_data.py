"""srsd_data.py — SRSD-Feynman 数据接入(官方公式类,零改动加载)。

机制:
- 上游仓库浅克隆到 %TEMP%/srsd-benchmark(官方 mit 协议);
- datasets/feynman.py 的 120 个公式类 import 时依赖 torchdistill(仅文件工具)
  与 graphviz(仅画图)→ 以 stub 模块注入绕开,公式代码本身一行不改;
- H0/H1 分类不改上游:由本模块读每类的 sympy 表达式 + 官方采样范围,
  按"严格单项式(所有指数有理化为非零)/非单项式"预注册划分,判据与
  任何检验结果无关(写死在分类函数里,先分类后实验)。
- 壁纸式筛选:输入输出严格为正 + log-log 适配(var>0, y>0, 变量数<=3)。
"""
import os
import sys
import tempfile
import types

import numpy as np
import pandas as pd

_UPSTREAM = os.path.join(tempfile.gettempdir(), "srsd-benchmark")


def _ensure_upstream():
    if not os.path.isdir(_UPSTREAM):
        import subprocess
        subprocess.run(
            ["git", "clone", "--depth", "1", "--quiet",
             "https://github.com/omron-sinicx/srsd-benchmark",
             _UPSTREAM], check=True, timeout=300)


def _install_stubs():
    """绕开 torchdistill / graphviz(公式类只用它们的边角工具)。"""
    if "torchdistill" in sys.modules:
        return
    td = types.ModuleType("torchdistill")
    cf = types.ModuleType("torchdistill.common")
    fu = types.ModuleType("torchdistill.common.file_util")

    def make_parent_dirs(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    fu.make_parent_dirs = make_parent_dirs
    yu = types.ModuleType("torchdistill.common.yaml_util")
    import yaml
    yu.load_yaml_file = lambda p: yaml.safe_load(
        open(p, encoding="utf-8"))  # noqa: SIM115 — 短生命周期桩模块,进程即退
    cf.file_util = fu
    cf.yaml_util = yu
    td.common = cf
    sys.modules["torchdistill"] = td
    sys.modules["torchdistill.common"] = cf
    sys.modules["torchdistill.common.file_util"] = fu
    sys.modules["torchdistill.common.yaml_util"] = yu

    gv = types.ModuleType("graphviz")
    gv.Digraph = gv.Source = object
    sys.modules["graphviz"] = gv


def load_feynman():
    """返回 FEYNMAN_EQUATION_CLASS_DICT(name->class)。"""
    _ensure_upstream()
    _install_stubs()
    if _UPSTREAM not in sys.path:
        sys.path.insert(0, _UPSTREAM)
    for stale in [k for k in list(sys.modules) if k.startswith("datasets")]:
        del sys.modules[stale]
    import datasets.feynman as feynman_mod
    return feynman_mod.FEYNMAN_EQUATION_CLASS_DICT


# ---------------------------------------------------------------- 分类与筛选
def _is_monomial(sympy_eq):
    """严格单项式判定:表达式结构是有限个符号的常数次幂之积。
    实现:逐符号做 n(log f / log x) 微商,能全为常数 → 单项式。"""
    import sympy
    f = sympy.sympify(sympy_eq)
    syms = sorted(f.free_symbols, key=lambda s: s.name)
    if not syms:
        return False                      # 纯常数,不是回归任务
    for s in syms:
        try:
            d = sympy.simplify(sympy.diff(sympy.log(f), s) * s)
        except Exception:  # noqa: BLE001 — sympy 对病态式会抛各种异常,统一按"非单项式"处理
            return False
        if d.free_symbols:
            return False                  # 该变量的指数含其他符号 → 非单项式
    return True


def _draw(so, rng, n):
    """按官方采样类的语义抽 n 个点(与上游 sampling.py 逐类对齐)。"""
    lo, hi = float(so.min_value), float(so.max_value)
    name = type(so).__name__
    pos_only = not getattr(so, "uses_negative", True)
    neg_only = not getattr(so, "uses_positive", True)
    if name == "DefaultSampling":                  # log10 域均匀
        v = 10.0 ** rng.uniform(np.log10(max(abs(lo), 1e-300)),
                                np.log10(max(abs(hi), 1e-300)), size=n)
    elif name == "SimpleSampling":                 # 线性域均匀
        v = rng.uniform(lo, hi, size=n)
    elif name == "IntegerSampling":                # 整数均匀
        v = rng.integers(int(lo), int(hi), size=n).astype(float)
    else:                                          # 未知类:保守按线性
        v = rng.uniform(lo, hi, size=n)
    if neg_only:
        v = -np.abs(v)
    elif not pos_only:                             # 双侧:50/50 混合(上游语义)
        v = np.where(rng.random(n) < 0.5, v, -v)
    return v


def classify(eq_cls):
    """预注册分类:返回 (klass in {H0,H1,EXCLUDE}, reason)。
    判据只依赖公式结构与官方采样范围,与任何检验结果无关。"""
    inst = eq_cls()
    f = inst.sympy_eq
    rng = np.random.default_rng(20260903)
    ok_pos, ok_finite = True, True
    for _ in range(64):
        xs = [_draw(so, rng, 1)[0] for so in inst.sampling_objs]
        try:
            y = float(inst.eq_func(np.array(xs, dtype=float)))
        except Exception:  # noqa: BLE001 — 病态采样点按不可用处理,重抽/剔除
            ok_finite = False
            break
        if not np.isfinite(y):
            ok_finite = False
            break
        if y <= 0 or any(v <= 0 for v in xs):
            ok_pos = False
            break
    if not ok_finite:
        return "EXCLUDE", "nonfinite_under_official_sampling"
    if not ok_pos:
        return "EXCLUDE", "nonpositive_inputs_or_output"
    if inst.get_var_count() > 3:
        return "EXCLUDE", "too_many_variables_for_harness"
    if _is_monomial(f):
        return "H0", "strict_monomial"
    return "H1", "non_monomial_structure"


def registry_with_meta():
    """返回 list of dict: 每个公式的名称/类/klass/理由/变量数/采样范围/sympy。"""
    import sympy
    out = []
    for name, cls in load_feynman().items():
        inst = cls()
        klass, reason = classify(cls)
        ranges = [(float(so.min_value), float(so.max_value),
                   type(so).__name__,
                   bool(getattr(so, "uses_positive", True)),
                   bool(getattr(so, "uses_negative", True)))
                  for so in inst.sampling_objs]
        out.append({"key": name, "eq_name": inst._eq_name, "cls": cls,
                    "klass": klass, "reason": reason,
                    "n_vars": inst.get_var_count(),
                    "sampling_ranges": ranges,
                    "sympy_str": sympy.srepr(inst.sympy_eq)[:200]})
    return out


def make_dataset_fn(eq_cls, n_vars_limit=3):
    """构造 run_experiment 兼容的场景函数 fn(seed,n,rel)->(df, feats)。
    严格按官方采样范围(log10 域均匀),乘性相对高斯噪声(rel),
    输出列名 target,x1..xk(k<=3)。y<=0 或非有限 → 抛 ValueError(上游
    数据性质,登记失败,不静默)。

    hetero=True 时改用方差随 log 预测值线性增大的异方差噪声
    (SD 从 0.3·rel 单调升到 1.7·rel,与合成台架 H1_VAR_GROWING 同构);
    range_scale 对每个变量的官方 log10 域做乘性缩放(1.0=官方,0.5=半程,
    2.0=两倍),注意 2.0 倍域可能使原本正值的函数出负值——那是真实性质
    (敏感性结论的一部分),由调用方登记为 not-applicable,不做平移。"""
    inst = eq_cls()

    def fn(seed, n, rel, hetero=False, range_scale=1.0):
        rng = np.random.default_rng(seed * 7919 + 104729)
        cols = []
        for so in inst.sampling_objs:
            lo, hi = float(so.min_value), float(so.max_value)
            name = type(so).__name__
            neg_only = not getattr(so, "uses_positive", True)
            pos_only = not getattr(so, "uses_negative", True)
            if name == "DefaultSampling":
                v = 10.0 ** rng.uniform(np.log10(max(abs(lo) * range_scale, 1e-300)),
                                        np.log10(max(abs(hi) * range_scale, 1e-300)),
                                        size=n)
            else:                                   # 线性/整数域
                v = rng.uniform(lo * range_scale, hi * range_scale, size=n)
            if neg_only:
                v = -np.abs(v)
            elif not pos_only:
                v = np.where(rng.random(n) < 0.5, v, -v)
            cols.append(v)
        X = np.column_stack(cols)
        with np.errstate(all="ignore"):
            y = np.array([float(inst.eq_func(X[i])) for i in range(n)])
        if (not np.isfinite(y).all()) or (not np.isfinite(X).all()):
            raise ValueError("nonfinite under official sampling")
        if (y <= 0).any() or (X <= 0).any():
            raise ValueError("nonpositive under official sampling")
        ly_c = np.log(y)
        if hetero:
            span = ly_c.max() - ly_c.min()
            sd_vec = rel * (0.3 + 1.4 * (ly_c - ly_c.min())
                            / (span if span > 1e-12 else 1.0))
            ly = ly_c + sd_vec * rng.standard_normal(n)
        else:
            ly = ly_c + rel * rng.standard_normal(n)
        data = {"target": np.exp(ly)}
        for j in range(X.shape[1]):
            data[f"x{j+1}"] = X[:, j]
        return pd.DataFrame(data), ",".join(f"x{j+1}" for j in range(X.shape[1]))

    return fn
