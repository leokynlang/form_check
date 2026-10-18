"""Prospective positive-domain task registry; no fitted outcomes are consulted.

All registry entries are JSON data. Source equations preserve their upstream
expressions and coefficients. Designed response laws are synthetic challenges,
not additional physical discoveries. A latent truth is the noiseless log-location
target; ``sample`` adds independent Gaussian noise in log target space.
"""
from __future__ import annotations

import json
import math
from functools import cache, lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special

import srsd_data

ROOT = Path(__file__).resolve().parents[2]

_OLD_GROUPS = {
    "feynman-i.10.7": "lorentz_radical",
    "feynman-i.48.2": "lorentz_radical",
    "feynman-ii.13.23": "lorentz_radical",
    "feynman-ii.13.34": "lorentz_radical",
    "feynman-i.26.2": "oscillatory_interference",
    "feynman-i.30.5": "oscillatory_interference",
    "feynman-i.37.4": "oscillatory_interference",
    "feynman-bonus.15": "oscillatory_interference",
    "feynman-bonus.20": "oscillatory_interference",
    "feynman-i.27.6": "rational_denominator",
    "feynman-i.39.11": "rational_denominator",
    "feynman-i.43.43": "rational_denominator",
    "feynman-ii.11.27": "rational_denominator",
    "feynman-ii.11.28": "rational_denominator",
    "feynman-ii.38.14": "rational_denominator",
    "feynman-i.41.16": "bose_planck",
    "feynman-iii.4.32": "bose_planck",
    "feynman-iii.4.33": "bose_planck",
    "feynman-bonus.5": "sum_radical",
}


@lru_cache(maxsize=1)
def _source_instances():
    return {cls()._eq_name: cls() for cls in srsd_data.load_feynman().values()}


def _bounds_at(bounds, scale, quantiles):
    lo, hi = map(float, bounds)
    if scale == "log":
        return [float(np.exp(np.log(lo) + q * (np.log(hi) - np.log(lo))))
                for q in quantiles]
    values = [lo + q * (hi - lo) for q in quantiles]
    if scale == "integer":
        return [math.ceil(values[0]), math.floor(values[1])]
    return [float(v) for v in values]


def _scales(instance):
    return [{"DefaultSampling": "log", "SimpleSampling": "linear",
             "IntegerSampling": "integer"}[type(obj).__name__]
            for obj in instance.sampling_objs]


def _source_task(name, klass, split, group, outer, train, *, source_kind,
                 validity, regime="fixed_interior", original=None):
    instance = _source_instances()[name]
    scales = _scales(instance)
    return {
        "task": name, "klass": klass, "split": split,
        "structural_group": group, "source_kind": source_kind,
        "source_equations": [name], "expression": str(instance.sympy_eq),
        "n_features": len(outer), "outer_bounds": outer,
        "train_bounds": train, "coordinate_scales": scales,
        "validity_argument": validity,
        "domain_policy": "Declared closed positive boxes fixed before sampling; no rejection sampling; latent truth is noiseless log-location.",
        "regime_kind": regime, "family": "upstream_source",
        "params": {"source_id": name},
        "integer_support": "inclusive integer endpoints; uniform integer draws; grids enumerate registered integer coordinates",
        "original_sampling": original,
    }


def _old_tasks():
    config = json.loads((ROOT / "runs/srsd_full/config.json").read_text(encoding="utf-8"))
    tasks = []
    for name, meta in sorted(config["task_pool"].items()):
        instance = _source_instances()[name]
        scales = _scales(instance)
        original = [[float(obj.min_value), float(obj.max_value)]
                    for obj in instance.sampling_objs]
        if len(scales) != meta["n_vars"]:
            raise ValueError(f"Source variable count changed: {name}")
        outer = [_bounds_at(b, s, (.05, .95)) for b, s in zip(original, scales, strict=False)]
        train = [_bounds_at(b, s, (.25, .75)) for b, s in zip(original, scales, strict=False)]
        validity = "All coordinates are strictly positive on the declared compact box; products, positive constants and positive powers are finite."
        if name == "feynman-i.37.4":
            outer[2], train[2] = [.1, 1.4], [.425, 1.075]
            validity += " The angle is at most1.4<pi/2, so cos(angle)>0 and the interference sum is strictly positive."
        if name in {"feynman-i.26.2", "feynman-i.30.5"}:
            validity += " Angular coordinates are strictly inside(0,pi/2); sine numerators and denominators are strictly positive."
        if name in {"feynman-i.39.11", "feynman-i.43.43"}:
            outer[0][0] = max(1.05, outer[0][0])
            validity += " The shifted-pole coordinate is at least1.05, so x0-1>=0.05."
        if _OLD_GROUPS.get(name) == "lorentz_radical":
            velocity_index = 1
            validity += f" Speed is at most{outer[velocity_index][1]:.17g}<299792458;1-(v/c)^2 has a strict positive margin."
        if name in {"feynman-ii.11.27", "feynman-ii.11.28"}:
            margin = 1 - outer[0][1] * outer[1][1] / 3
            if margin <= 0:
                raise ValueError(f"Invalid dielectric source box: {name}")
            validity += f" The denominator1-x0*x1/3 is at least{margin:.17g}>0."
        if _OLD_GROUPS.get(name) == "bose_planck":
            validity += " The exponential argument is strictly positive; expm1(argument)>0. Stable expm1 evaluation preserves the exact source law."
        if name == "feynman-bonus.15":
            beta = outer[0][1] / 299792458
            validity += f" beta<={beta:.17g}<1; radical positive and denominator1+beta*cos(theta)>=1-beta>0."
        if name == "feynman-bonus.20":
            validity += " Writing r=x0/x1>0, the bracket is r+1/r-sin(theta)^2>=1."
        if name == "feynman-bonus.5":
            validity += " Both masses are strictly positive, so their sum and the square-root denominator are positive."
        group = "H0_null_family" if meta["klass"] == "H0" else _OLD_GROUPS[name]
        rec = _source_task(name, meta["klass"], "development", group, outer, train,
                           source_kind="reused_source_restricted_domain", validity=validity,
                           original=original)
        if meta["klass"] == "H0":
            rec["design_group"] = f"legacy_development_{name}"
        tasks.append(rec)
    return tasks


def _new_source_tasks():
    specs = [
        ("feynman-i.6.20", "gaussian", [[.1, 2.], [.8, 1.5]],
         "sigma>=0.8; Gaussian exponential strictly positive, exponent>=-3.125; normalizing denominator is positive."),
        ("feynman-i.6.20a", "gaussian", [[.1, 2.]],
         "The Gaussian exponential is strictly positive, exponent>=-2, and its normalization is a positive constant."),
        ("feynman-i.6.20b", "gaussian", [[.2, 2.], [.1, .5], [.8, 1.5]],
         "sigma>=0.8; |x-mu|<=1.9, so Gaussian exponent>=-2.8203125; output strictly positive."),
        ("feynman-i.44.4", "positive_log", [[1e25, 1e27], [10., 1000.], [3e-4, 1e-3], [1e-5, 1e-4]],
         "All dimensional coordinates are positive; V2/V1>=3, so log(V2/V1)>=log(3)>0; the original Boltzmann prefactor is unchanged."),
        ("feynman-ii.35.18", "hyperbolic_response", [[1e23, 1e25], [1e-25, 1e-23], [.001, .1], [10., 1000.]],
         "All coordinates positive;0<argument<=0.007242971; exp(argument)+exp(-argument)>=2, with no overflow or pole."),
        ("feynman-ii.35.21", "hyperbolic_response", [[1e23, 1e25], [1e-25, 1e-23], [.001, .1], [10., 1000.]],
         "All coordinates positive;0<argument<=0.007242971; tanh(argument)>0 and all original dimensional prefactors are retained."),
    ]
    tasks = []
    for name, group, outer, validity in specs:
        scales = _scales(_source_instances()[name])
        train = [_bounds_at(b, s, (.25, .75)) for b, s in zip(outer, scales, strict=False)]
        tasks.append(_source_task(name, "H1", "confirmation", group, outer, train,
                                  source_kind="new_source_positive_restriction", validity=validity))
    return tasks


def _designed_tasks():
    # Four windows are population-level design regimes, never independent
    # mechanism units. The fixed outer box is the same for all four windows.
    families = {
        "erf": ("gaussian", "confirmation", "erf(x0)", [.05, 2.],
                "erf(x)>0 for x>0, since it is the integral of a strictly positive Gaussian density."),
        "atan": ("inverse_trig_phase", "confirmation", "atan(x0)", [.05, 3.],
                 "atan(x)>0 and finite for every x>0."),
        "asinh": ("inverse_hyperbolic", "reserve", "asinh(x0)", [.05, 3.],
                  "asinh(x)=log(x+sqrt(1+x^2))>0 for x>0; no branch ambiguity."),
        "lambertw": ("inverse_exponential", "reserve", "LambertW(x0)", [.05, 3.],
                     "The principal Lambert W branch is unique, real, positive and finite for x>0."),
        "ellipk": ("elliptic_integral", "confirmation", "elliptic_k(x0)", [.02, .7],
                   "0<m<=0.7;1-m*sin(theta)^2>=0.3 throughout the defining integral, hence finite positive truth."),
        "besseli0": ("bessel_integral", "confirmation", "besseli(0,x0)", [.05, 3.],
                     "I0(x)=pi^-1*integral(exp(x*cos(theta)),theta=0..pi); the integrand is positive and bounded for the compact box."),
        "threshold": ("threshold_response", "stress", "x0*(1+Max(x0-1,0)**2)", [.1, 3.],
                      "x>0 and1+max(x-1,0)^2>=1; continuous positive piecewise response. Below-threshold training has no information about the external regime."),
        "resonance": ("damped_resonance", "confirmation", "1/((1-x0**2)**2+x1**2*x0**2)", [.1, 1.6],
                      "x>0 and gamma>=0.1; denominator is a sum of squares with gamma^2*x^2>=0.0001; no resonance pole."),
        "gompertz": ("nested_exponential", "confirmation", "exp(-exp(-x0))", [.05, 3.],
                    "For x>0,0<exp(-exp(-x))<1; both exponential levels are finite on this compact box."),
        "damped_oscillation": ("oscillatory_interference", "development", "2+exp(-x0)*cos(x1*x0)", [.05, 3.],
                              "x>0;|exp(-x)*cos(omega*x)|<1, hence truth is strictly greater than1. Grouped conservatively with old angular/interference formulas."),
        "softplus": ("positive_log", "confirmation", "log(1+exp(x0))", [.05, 3.],
                     "log(1+exp(x))>0 for all finite x; evaluated using logaddexp for numerical stability."),
        "sech": ("hyperbolic_response", "confirmation", "1/cosh(x0)", [.05, 2.],
                 "cosh(x)>=1, so sech(x)>0; arguments are compact and cannot overflow."),
    }
    tasks = []
    windows = [( .15, .45), (.30, .60), (.45, .75), (.20, .80)]
    for family, (group, split, formula, outer_x, validity) in families.items():
        for variant, window in enumerate(windows, 1):
            outer = [list(outer_x)]
            scales = ["linear" if family == "ellipk" else "log"]
            train = [_bounds_at(outer_x, scales[0], window)]
            regime = "shape_window"
            if family == "threshold":
                train[0] = [[.2, .5], [.4, .8], [.8, 1.3], [1.1, 2.]][variant - 1]
                regime = "external_only_shift" if variant <= 2 else "training_threshold_or_post_threshold"
            if family == "resonance":
                outer.append([.1, 1.])
                scales.append("log")
                train.append([.3, .6])
            elif family == "damped_oscillation":
                outer.append([.2, 2.])
                scales.append("log")
                train.append([.5, 1.2])
            outer.append([.1, 10.])
            scales.append("log")
            train.append([.3, 3.])
            amp = len(outer) - 1
            tasks.append({
                "task": f"designed-{family}-r{variant}", "klass": "H1", "split": split,
                "structural_group": group, "source_kind": "synthetic_positive_mechanism",
                "source_equations": [], "expression": f"x{amp}*({formula})",
                "n_features": len(outer), "outer_bounds": outer, "train_bounds": train,
                "coordinate_scales": scales, "validity_argument": validity + " The amplitude A is strictly positive.",
                "domain_policy": "Fixed dimensionless synthetic challenge box; common amplitude rule A_outer=[0.1,10],A_train=[0.3,3] fixed for every family before outcomes; no rejection sampling.",
                "regime_kind": regime, "family": family,
                "params": {"regime": variant, "amplitude_index": amp, "threshold": 1., "threshold_power": 2.},
                "integer_support": None,
            })
    return tasks


def _null_tasks():
    tasks = []
    for split, family, offset in [("calibration", "log_uniform", 0), ("confirmation", "log_beta", 1)]:
        for dim in range(1, 5):
            for variant in range(1, 4):
                # Odd numerator /17 for calibration, even numerator /19 for
                # confirmation. Exponents never coincide with old integer or
                # half-integer nulls, or with the other role's signatures.
                denominator = 17 if offset == 0 else 19
                numerators = [(2 * (variant + 3 * j) + 1 + offset) * (-1 if j % 2 else 1)
                              for j in range(dim)]
                exponents = [v / denominator for v in numerators]
                coefficient = (1. + variant / 7.) if offset == 0 else (1. + variant / 11.)
                expression = str(coefficient) + "*" + "*".join(
                    f"x{j}**({num}/{denominator})" for j, num in enumerate(numerators))
                tasks.append({
                    "task": f"null-{family}-p{dim}-r{variant}", "klass": "H0", "split": split,
                    "structural_group": "H0_null_family",
                    "design_group": f"{split}_{family}_p{dim}",
                    "source_kind": "synthetic_null_design", "source_equations": [],
                    "expression": expression, "n_features": dim,
                    "outer_bounds": [[.1, 10.] for _ in range(dim)],
                    "train_bounds": [[.3, 3.] for _ in range(dim)],
                    "coordinate_scales": ["log"] * dim,
                    "validity_argument": "Positive coefficient and strictly positive compact input boxes ensure finite positive monomials; exact log-affine conditional mean of log target.",
                    "domain_policy": "H0 is one shared functional family; only design distributions and exponent signatures are held out. Calibration uses uniform log coordinates; confirmation uses Beta(2,2) log coordinates.",
                    "regime_kind": "null_design", "family": "monomial_null",
                    "params": {"coefficient": coefficient, "exponents": exponents,
                               "exponent_numerators": numerators, "exponent_denominator": denominator,
                               "design_distribution": family}, "integer_support": None,
                })
    return tasks


def build_registry():
    """Return an outcome-independent, deterministic, JSON-only task manifest."""
    tasks = _old_tasks() + _new_source_tasks() + _designed_tasks() + _null_tasks()
    if len(tasks) != 110 or len({r["task"] for r in tasks}) != 110:
        raise AssertionError("Unexpected registry cardinality or duplicate IDs")
    json.dumps(tasks, allow_nan=False)
    for task in tasks:
        if len(task["train_bounds"]) != task["n_features"]:
            raise AssertionError("Bad feature count")
        for outer, train in zip(task["outer_bounds"], task["train_bounds"], strict=False):
            if not 0 < outer[0] <= train[0] < train[1] <= outer[1]:
                raise AssertionError(f"Invalid positive containment: {task['task']}")
    return tasks


@cache
def _source_function(name):
    import sympy as sp

    expr = _source_instances()[name].sympy_eq
    # Replace precisely exp(t)-1, preserving all prefactors and arguments.
    def is_expm1(node):
        return (isinstance(node, sp.Add) and len(node.args) == 2
                and sp.Integer(-1) in node.args
                and any(arg.func == sp.exp for arg in node.args))

    def replace_expm1(node):
        exponential = next(arg for arg in node.args if arg.func == sp.exp)
        return sp.Function("expm1")(exponential.args[0])

    expr = expr.replace(is_expm1, replace_expm1)
    symbols = sorted(expr.free_symbols, key=lambda symbol: int(str(symbol)[1:]))
    return sp.lambdify(symbols, expr, modules=[{"expm1": np.expm1}, "numpy"], cse=True)


def evaluate(task, X):
    """Evaluate the exact registered noiseless truth at n-by-p coordinates."""
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] != task["n_features"]:
        raise ValueError("X must be an n-by-registered-feature-count array")
    if not np.isfinite(X).all() or np.any(X <= 0):
        raise ValueError("Evaluation coordinates must be finite and strictly positive")
    if task["family"] == "upstream_source":
        result = _source_function(task["params"]["source_id"])(*[X[:, j] for j in range(X.shape[1])])
    elif task["family"] == "monomial_null":
        result = np.exp(np.log(task["params"]["coefficient"]) + np.log(X) @ np.asarray(task["params"]["exponents"]))
    else:
        x, family = X[:, 0], task["family"]
        if family == "erf":
            shape = special.erf(x)
        elif family == "atan":
            shape = np.arctan(x)
        elif family == "asinh":
            shape = np.arcsinh(x)
        elif family == "lambertw":
            shape = special.lambertw(x, k=0).real
        elif family == "ellipk":
            shape = special.ellipk(x)
        elif family == "besseli0":
            shape = special.i0(x)
        elif family == "threshold":
            shape = x * (1 + np.maximum(x - 1., 0.) ** 2)
        elif family == "resonance":
            shape = 1 / ((1 - x ** 2) ** 2 + X[:, 1] ** 2 * x ** 2)
        elif family == "gompertz":
            shape = np.exp(-np.exp(-x))
        elif family == "damped_oscillation":
            shape = 2 + np.exp(-x) * np.cos(X[:, 1] * x)
        elif family == "softplus":
            shape = np.logaddexp(0., x)
        elif family == "sech":
            shape = 1 / np.cosh(x)
        else:
            raise ValueError(f"Unknown registered family: {family}")
        result = X[:, task["params"]["amplitude_index"]] * shape
    result = np.broadcast_to(np.asarray(result, dtype=float), (len(X),)).copy()
    if not np.isfinite(result).all() or np.any(result <= 0):
        raise ValueError(f"Nonpositive or nonfinite registered truth: {task['task']}")
    return result


def sample(task, seed, n, rel):
    """Draw one fresh design without rejection; noise is Gaussian in log target."""
    rng = np.random.default_rng(seed)
    columns = []
    beta_design = task.get("params", {}).get("design_distribution") == "log_beta"
    for bounds, scale in zip(task["train_bounds"], task["coordinate_scales"], strict=False):
        lo, hi = bounds
        if scale == "integer":
            values = rng.integers(int(lo), int(hi) + 1, size=n).astype(float)
        else:
            coordinates = rng.beta(2., 2., size=n) if beta_design else rng.random(n)
            values = (np.exp(np.log(lo) + coordinates * (np.log(hi) - np.log(lo)))
                      if scale == "log" else lo + coordinates * (hi - lo))
        columns.append(values)
    X = np.column_stack(columns)
    truth = evaluate(task, X)
    target = np.exp(np.log(truth) + float(rel) * rng.standard_normal(n))
    if not np.isfinite(target).all() or np.any(target <= 0):
        raise ValueError("Observation noise generated invalid positive targets")
    features = [f"x{j + 1}" for j in range(task["n_features"])]
    data = pd.DataFrame({"target": target, **{name: X[:, j] for j, name in enumerate(features)}})
    return data, features
