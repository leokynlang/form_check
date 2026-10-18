"""Structure-only inventory of all upstream SRSD equations for fixed-domain work.

No samples, candidate fits, gate statistics, or diagnostic outcomes are computed.
Historical inclusion is read from the existing immutable registry. Suggested
groups deliberately err toward over-grouping; the protocol must freeze the final
group assignments before generating P0-2 data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import itertools
import json
from collections import defaultdict
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _normalize(expr):
    """Remove a global constant multiplier, retaining all internal constants."""
    expr = sp.factor_terms(expr)
    variables = tuple(sorted(expr.free_symbols, key=str))
    if variables:
        _, expr = expr.as_independent(*variables, as_Add=False)
    return expr


def _tree(expr, ignore_numbers=False, ignore_variables=False, exponent=False):
    if isinstance(expr, sp.Symbol):
        return ("V",) if ignore_variables else ("S", str(expr))
    if not expr.free_symbols:
        return ("N", "C" if ignore_numbers and not exponent else sp.srepr(expr))
    if isinstance(expr, sp.Pow):
        return ("Pow", _tree(expr.base, ignore_numbers, ignore_variables),
                _tree(expr.exp, ignore_numbers, ignore_variables, exponent=True))
    children = [_tree(arg, ignore_numbers, ignore_variables) for arg in expr.args]
    if isinstance(expr, (sp.Add, sp.Mul)):
        children.sort(key=repr)
    return (expr.func.__name__, *children)


def _rename_tree(tree, rename):
    if tree[0] == "S":
        return ("V", rename[tree[1]])
    if tree[0] in {"V", "N"}:
        return tree
    children = [_rename_tree(child, rename) for child in tree[1:]]
    if tree[0] in {"Add", "Mul"}:
        children.sort(key=repr)
    return (tree[0], *children)


def canonical_structure(expr, ignore_numbers=False):
    """Canonical under every variable renaming and a nonzero scalar prefactor.

    All 120 upstream tasks have <=8 variables, so exhaustive permutations are
    practical and avoid a heuristic graph-isomorphism claim. This is syntactic
    canonicalization after factor_terms, not general algebraic equivalence.
    """
    expr = _normalize(expr)
    names = sorted(str(s) for s in expr.free_symbols)
    tree = _tree(expr, ignore_numbers=ignore_numbers)
    return min(repr(_rename_tree(tree, dict(zip(names, perm, strict=True))))
               for perm in itertools.permutations(range(len(names))))


def _monomial_structure(expr):
    """Pure power-law structure, without running the historical sample filter."""
    if not expr.free_symbols:
        return True
    if isinstance(expr, sp.Symbol):
        return True
    if isinstance(expr, sp.Mul):
        return all(_monomial_structure(a) for a in expr.args)
    if isinstance(expr, sp.Pow):
        return not expr.exp.free_symbols and _monomial_structure(expr.base)
    return False


def nonlinear_core(expr):
    """Conservative near-twin key, ignoring separable monomial multipliers.

    Constants are masked except exponents. Variables become generic leaves.
    E.g. m/sqrt(1-v^2/c^2) and p*v/sqrt(1-v^2/c^2) share a core. An identical
    key is a leakage warning, not a claim that tasks are independent or equal.
    """
    expr = _normalize(expr)
    if _monomial_structure(expr):
        return "MONOMIAL"
    if isinstance(expr, sp.Mul):
        args = [a for a in expr.args if not _monomial_structure(a)]
        if len(args) == 1:
            return nonlinear_core(args[0])
        return repr(("MulCore", *sorted(nonlinear_core(a) for a in args)))
    return repr(_tree(expr, ignore_numbers=True, ignore_variables=True))


def features(expr):
    nodes = list(sp.preorder_traversal(expr))
    numerator, denominator = sp.fraction(sp.together(expr))
    unary = sorted({n.func.__name__ for n in nodes if isinstance(n, sp.Function)})
    return {
        "strict_monomial_structure": _monomial_structure(expr),
        "contains_sum": any(isinstance(n, sp.Add) for n in nodes),
        "contains_product": any(isinstance(n, sp.Mul) for n in nodes),
        "contains_ratio": bool(denominator.free_symbols),
        "contains_noninteger_power": any(
            isinstance(n, sp.Pow) and not n.exp.free_symbols and n.exp.is_integer is not True
            for n in nodes),
        "unary_functions": unary,
        "numerator": str(numerator), "denominator": str(denominator),
    }


def connected_groups(records):
    """Union exact/constant-masked/core keys and same source equation variants."""
    parent = {r["eq_name"]: r["eq_name"] for r in records}

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    edges = []
    for key in ("canonical_structure", "constant_masked_structure", "nonlinear_core", "source_family"):
        by_key = defaultdict(list)
        for record in records:
            by_key[record[key]].append(record["eq_name"])
        for members in by_key.values():
            for name in members[1:]:
                first = members[0]
                if find(first) != find(name):
                    parent[find(name)] = find(first)
                edges.append({"a": first, "b": name, "reason": key})
    groups = defaultdict(list)
    for record in records:
        groups[find(record["eq_name"])].append(record)
    output = []
    for index, members in enumerate(sorted(groups.values(), key=lambda x: min(r["eq_name"] for r in x))):
        gid = f"structural-{index:03d}"
        old = sorted(r["eq_name"] for r in members if r["historical32"])
        for record in members:
            record["suggested_group"] = gid
            record["historically_tainted_group"] = bool(old)
        output.append({"group": gid, "members": sorted(r["eq_name"] for r in members),
                       "historical32_members": old, "historically_tainted": bool(old)})
    return output, edges


def build_inventory():
    import srsd_data

    registry_path = ROOT / "runs/srsd_full/task_registry.csv"
    with registry_path.open(encoding="utf-8", newline="") as stream:
        old = {r["eq_name"]: r for r in csv.DictReader(stream)}
    config_path = ROOT / "runs/srsd_full/config.json"
    task_pool = json.loads(config_path.read_text(encoding="utf-8"))["task_pool"]
    records = []
    classes = srsd_data.load_feynman()
    for class_name, cls in classes.items():
        equation = cls()
        name, expr = equation._eq_name, equation.sympy_eq
        historical = old.get(name, {})
        frozen_historical = task_pool.get(name)
        source_family = name.rstrip("abtx")
        sampling = [{"class_name": type(obj).__name__, "attributes": vars(obj)}
                    for obj in equation.sampling_objs]
        record = {
            "eq_name": name, "class_name": class_name,
            "n_vars": equation.get_var_count(), "variables": [str(x) for x in equation.x],
            "expression": str(expr), "sympy_srepr": sp.srepr(expr),
            "upstream_docstring": inspect.getdoc(cls),
            "upstream_source_line": inspect.getsourcelines(cls)[1],
            "sampling_objects": sampling,
            "historical32": name in task_pool,
            "historical_class": frozen_historical["klass"] if frozen_historical else "EXCLUDE",
            "historical_reason": historical.get("reason", "NOT_REGISTERED"),
            "legacy_registry_class": historical.get("klass", "NOT_REGISTERED"),
            "legacy_registry_disagrees": bool(frozen_historical and historical.get("klass") != frozen_historical["klass"]),
            "source_family": source_family,
            "canonical_structure": canonical_structure(expr),
            "constant_masked_structure": canonical_structure(expr, ignore_numbers=True),
            "nonlinear_core": nonlinear_core(expr),
            **features(expr),
        }
        records.append(record)
    groups, edges = connected_groups(records)
    upstream_path = Path(inspect.getfile(next(iter(classes.values()))))
    result = {
        "scope": "Structure/source inventory only; no samples or diagnostic outcomes read.",
        "upstream_source": str(upstream_path), "upstream_sha256": sha256(upstream_path),
        "historical_registry": str(registry_path), "historical_registry_sha256": sha256(registry_path),
        "authoritative_historical_config": str(config_path), "authoritative_historical_config_sha256": sha256(config_path),
        "grouping_rule": "Connected components of exact, constant-masked, nonlinear-core and source-family keys; all historical32-tainted components are development only unless conservatively merged further.",
        "counts": {
            "all": len(records), "historical32": sum(r["historical32"] for r in records),
            "historical_excluded": sum(r["historical_class"] == "EXCLUDE" for r in records),
            "historical_H0": sum(r["historical_class"] == "H0" for r in records),
            "historical_H1": sum(r["historical_class"] == "H1" for r in records),
            "legacy_registry_disagreements": sum(r["legacy_registry_disagrees"] for r in records),
            "suggested_groups": len(groups),
            "historically_tainted_groups": sum(g["historically_tainted"] for g in groups),
            "untainted_groups": sum(not g["historically_tainted"] for g in groups),
            "new_at_most_3_variables": sum(not r["historical32"] and r["n_vars"] <= 3 for r in records),
        },
        "groups": groups, "group_edges": edges, "equations": records,
        "limitations": [
            "Identical keys are conservative leakage flags, not independent scientific mechanisms.",
            "Different keys do not guarantee dissimilarity; final manual grouping may merge components.",
            "Positive branch restriction changes original upstream sampling; it must be explicit in the protocol.",
            "All monomial structures are one component because they become the same log-linear hypothesis family.",
        ],
    }
    return result


def source_domain_candidates():
    """Analytic positive boxes for parent review, never final assignments.

    Bounds use only the formula and the original positive sampling branch.
    They do not use fitted-model errors, screening results or gate outcomes.
    """
    pi = float(sp.pi)
    return [
        {"eq_name": "feynman-i.6.20a", "outer_box": [[0.1, 3.0]], "proof": "exp(-theta^2/2)/sqrt(2*pi)>0 for finite theta; bounded exponent avoids underflow.", "merge_note": "All Gaussian variants together."},
        {"eq_name": "feynman-i.6.20", "outer_box": [[0.1, 2.0], [1.0, 10.0]], "proof": "sigma>=1, |theta/sigma|<=2, positive Gaussian density without underflow.", "merge_note": "All Gaussian variants together."},
        {"eq_name": "feynman-i.6.20b", "outer_box": [[0.1, 2.0], [0.1, 2.0], [1.0, 10.0]], "proof": "sigma>=1, |(theta-theta1)/sigma|<=1.9, positive Gaussian density.", "merge_note": "All Gaussian variants together."},
        {"eq_name": "feynman-i.15.3t", "outer_box": [[5e-5, 1e-4], [1e5, 1e7], [1.0, 100.0]], "proof": "v<c and t-v*x/c^2 >=5e-5-1.12e-8>0.", "merge_note": "Conservatively merge with all historical Lorentz/Doppler structures; development only."},
        {"eq_name": "feynman-i.15.3x", "outer_box": [[60.0, 100.0], [1e6, 1e7], [1e-6, 5e-6]], "proof": "v<c and x-v*t>=60-50>0.", "merge_note": "Conservatively merge with all historical Lorentz/Doppler structures; development only."},
        {"eq_name": "feynman-i.16.6", "outer_box": [[1e6, 1e8], [1e6, 1e8]], "proof": "Both velocities positive and below c; denominator 1+v*w/c^2>1 and composed speed<c.", "merge_note": "Relativistic family; conservatively merge historical Lorentz/Doppler group."},
        {"eq_name": "feynman-i.18.12", "outer_box": [[0.1, 10.0], [0.1, 10.0], [pi/6, 5*pi/6]], "proof": "All factors positive and sin(theta)>=1/2.", "merge_note": "Merge positive sine/cosine phase shifts; historical trigonometric group may taint."},
        {"eq_name": "feynman-i.30.3", "outer_box": [[10.0, 1000.0], [10, 20], [0.005, 0.02]], "proof": "Integer N in [10,20], phi/2 in [.0025,.01], N*phi/2 in [.025,.2]; neither sine vanishes.", "sampling_note": "Retain integer sampling for N; do not silently turn slit count into continuous feature.", "merge_note": "Trigonometric ratio/diffraction family; possible old sine-ratio overlap."},
        {"eq_name": "feynman-i.34.10", "outer_box": [[1e9, 1e11], [1e5, 1e7]], "proof": "nu>0 and 1-v/c>=1-1e7/c>0.", "merge_note": "Merge Doppler/Lorentz family for conservative historical exposure."},
        {"eq_name": "feynman-i.34.14", "outer_box": [[1e6, 1e8], [1e9, 1e11]], "proof": "v<c, nu>0; both (1+v/c) and sqrt(1-v^2/c^2) positive.", "merge_note": "Merge Doppler/Lorentz family for conservative historical exposure."},
        {"eq_name": "feynman-ii.6.11", "outer_box": [[1e-22, 1e-20], [pi/12, 5*pi/12], [1e-10, 1e-8]], "proof": "Positive charge/radius and cos(theta)>=cos(5*pi/12)>0.", "merge_note": "Positive sine/cosine phase-shift family."},
        {"eq_name": "feynman-ii.6.15b", "outer_box": [[1e-22, 1e-20], [pi/12, 5*pi/12], [1e-10, 1e-8]], "proof": "Positive charge/radius, theta in (0,pi/2), so sin(theta)*cos(theta)>0.", "merge_note": "sin(theta)cos(theta)=sin(2theta)/2; merge sine/cosine family."},
        {"eq_name": "feynman-ii.15.4", "outer_box": [[1e-25, 1e-23], [0.001, 0.1], [2*pi/3, 5*pi/6]], "proof": "Positive dipole/field and cos(theta)<0 imply -p*E*cos(theta)>0.", "merge_note": "Exact formula twin II.15.5; positive phase-shift sine/cosine family."},
        {"eq_name": "feynman-ii.15.5", "outer_box": [[1e-22, 1e-20], [10.0, 1000.0], [2*pi/3, 5*pi/6]], "proof": "Positive factors and cos(theta)<0 imply positive product with leading minus.", "merge_note": "Exact formula twin II.15.4; positive phase-shift sine/cosine family."},
        {"eq_name": "feynman-ii.21.32", "outer_box": [[0.001, 0.1], [1.0, 100.0], [1e6, 1e8]], "proof": "All factors positive and v<c keeps Doppler denominator positive.", "merge_note": "Same nonlinear core as I.34.10; merge Doppler/Lorentz historical family."},
        {"eq_name": "feynman-ii.24.17", "outer_box": [[3e10, 1e11], [0.05, 0.1]], "proof": "omega^2/c^2-pi^2/d^2>=1.11265e-17*(3e10)^2-pi^2/.05^2>6000.", "merge_note": "Difference-of-squares square-root; review alongside bonus.7."},
        {"eq_name": "feynman-ii.37.1", "outer_box": [[1e-25, 1e-23], [1e4, 1e6], [0.001, 0.1]], "proof": "All factors positive, 1+x1>1.", "merge_note": "Affine-shift numerator versus historical shifted denominator differs mainly by sign in log space; conservative merge."},
        {"eq_name": "feynman-iii.8.54", "outer_box": [[1e-18, 1e-17], [1e-18, 1e-17]], "proof": "Phase is in approximately [.00948,.9483], inside (0,pi), so sin^2>0.", "merge_note": "sin^2(theta)=(1-cos(2theta))/2; merge III.15.12, diffraction and other phase-shift trig if broad grouping."},
        {"eq_name": "feynman-iii.15.12", "outer_box": [[1e-18, 1e-16], [0.1, 10.0], [1e-10, 1e-8]], "proof": "Analytically phase is in [1e-11,1e-7], so 1-cos(phase)>0.", "numerical_caveat": "Original numpy 1-cos rounds to zero at small phase. Requires frozen algebraically equivalent 2*sin(phase/2)^2 evaluator before reuse; do not silently sample-reject.", "merge_note": "Exact trig identity twin of squared sine up to phase mapping."},
        {"eq_name": "feynman-iii.17.37", "outer_box": [[1e-18, 1e-16], [1e-18, 1e-16], [pi/12, 5*pi/12]], "proof": "1+a*cos(theta)>0.", "numerical_caveat": "a<=1e-16 makes nonlinearity at/below double precision relative to the constant. Not a meaningful new H1 at original ranges in float64; keep inventory only unless explicitly changing physical design.", "merge_note": "Affine cosine family; do not call computationally collapsed curve a substantive H1."},
        {"eq_name": "feynman-bonus.8", "outer_box": [[1e-24, 1e-22], [pi/12, 11*pi/12]], "proof": "lambda>0 and 1+positive*(1-cos(theta))>=1.", "numerical_caveat": "Within original lambda range nonlinear denominator correction is <=2.45e-9; physically valid but close to monomial.", "merge_note": "Compton/trig-rational family; not guaranteed independent of trigonometric pool."},
        {"eq_name": "feynman-bonus.10", "outer_box": [[pi/12, pi/3], [1e6, 1e8]], "proof": "cos(theta)>=.5 and beta<=.334 imply numerator>.166; denominator>=1-.334>0.", "merge_note": "Relativistic aberration; merge historical Doppler/Lorentz family conservatively."},
        {"eq_name": "feynman-bonus.18", "outer_box": [[10.0, 1000.0], [1e8, 1e10], [1.0, 100.0]], "proof": "Both c^2*x0/x1^2 and x2^2 strictly positive.", "merge_note": "Positive sum-of-monomials family; compare with old bonus.5 and new synthetic sum/mixture group before split."},
    ]


def write_inventory(out):
    result = build_inventory()
    out.mkdir(parents=True, exist_ok=True)
    (out / "inventory.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    flat = []
    for r in result["equations"]:
        flat.append({k: r[k] for k in ["eq_name", "class_name", "n_vars", "historical32", "historical_class", "historical_reason", "suggested_group", "historically_tainted_group", "expression", "nonlinear_core"]})
    with (out / "inventory.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    lines = ["# P0-2 upstream structural inventory", "", result["scope"], "", "```json", json.dumps(result["counts"], indent=2), "```", "", "## Suggested conservative groups", ""]
    for g in result["groups"]:
        label = "development (historical exposure)" if g["historically_tainted"] else "new source structure"
        lines.append(f"- {g['group']} — {label}: {', '.join(g['members'])}")
    lines += ["", "## New, previously excluded, <=3-variable candidates", "", "These are candidates for domain review, not automatically valid tasks.", ""]
    for r in result["equations"]:
        if not r["historical32"] and r["n_vars"] <= 3:
            lines.append(f"- {r['eq_name']} ({r['suggested_group']}): `{r['expression']}`; original exclusion: {r['historical_reason']}.")
    (out / "INVENTORY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    candidates = source_domain_candidates()
    by_name = {r["eq_name"]: r for r in result["equations"]}
    for candidate in candidates:
        source = by_name[candidate["eq_name"]]
        candidate["expression"] = source["expression"]
        candidate["suggested_group"] = source["suggested_group"]
        candidate["historically_tainted_group"] = source["historically_tainted_group"]
    (out / "source_domain_review.json").write_text(json.dumps({
        "status": "Review proposals only, not a frozen protocol or final split.",
        "selection_basis": "Source formula, physical variable descriptions, positive official branches and analytic bounds only; no fits/outcomes.",
        "candidates": candidates,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="runs/p02_inventory")
    args = parser.parse_args()
    result = write_inventory(ROOT / args.out)
    print(json.dumps(result["counts"]))


if __name__ == "__main__":
    main()
