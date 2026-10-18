"""Structure-only leakage assertions for the pre-sampling P0-2 registry.

Manual family assignments are primary. Canonical syntax and nonlinear cores add
conservative checks; neither establishes general algebraic independence. Null
tasks deliberately share the log-linear functional family, so their resampling
unit is the frozen design group rather than the universal monomial group.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

import sympy as sp

from .inventory import canonical_structure, nonlinear_core

ROOT = Path(__file__).resolve().parents[2]


def group_unit(task):
    """Return the independent unit used in design/family bootstrap summaries."""
    klass = task.get("klass")
    if klass not in {"H0", "H1"}:
        raise ValueError(f"Invalid task class {klass!r} for {task.get('task')}")
    field, prefix = ("design_group", "design:") if klass == "H0" else ("structural_group", "structure:")
    value = task.get(field)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Task {task.get('task')} needs a nonempty {field}")
    return prefix + value


def _parse_expression(task):
    """Parse local registry expressions, including stable numerical spellings.

    Inputs are authored, trusted local protocol strings. No external document
    instruction or diagnostic result is read here. Scalar frozen parameters are
    substituted before structural comparisons; changing their window does not
    create an independent functional family.
    """
    expression = task.get("expression")
    if expression is None:
        raise ValueError(f"Missing expression for {task.get('task')}")
    local = {
        "log1p": lambda value: sp.log(1 + value),
        "expm1": lambda value: sp.exp(value) - 1,
        "hypot": lambda x, y: sp.sqrt(x*x + y*y),
        "abs": sp.Abs, "Abs": sp.Abs,
        "arcsinh": sp.asinh, "asinh": sp.asinh,
        "arccosh": sp.acosh, "acosh": sp.acosh,
        "arctanh": sp.atanh, "atanh": sp.atanh,
        "lambertw": sp.LambertW, "LambertW": sp.LambertW,
        "pi": sp.pi,
    }
    try:
        parsed = sp.sympify(expression, locals=local)
        if not isinstance(parsed, sp.Expr):
            raise TypeError("expression did not parse as a scalar SymPy expression")
        substitutions = {}
        parameters = task.get("params", {})
        if isinstance(parameters, dict):
            for symbol in parsed.free_symbols:
                value = parameters.get(str(symbol))
                if isinstance(value, (float, int, str)) and not isinstance(value, bool):
                    candidate = sp.sympify(value)
                    if candidate.is_number:
                        substitutions[symbol] = candidate
        return parsed.subs(substitutions)
    except (TypeError, ValueError, SyntaxError, sp.SympifyError) as exc:
        raise ValueError(f"Cannot parse frozen expression for {task.get('task')}: {exc}") from exc


@lru_cache(maxsize=1)
def _historical_inventory():
    path = ROOT / "runs/p02_inventory/inventory.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    from .inventory import build_inventory
    return build_inventory()


def _aliases(task):
    values = [task.get("task", "")]
    source = task.get("source_equations", [])
    if isinstance(source, str):
        source = [source]
    if source is None:
        source = []
    if not isinstance(source, (list, tuple)):
        raise ValueError(f"Invalid source_equations for {task.get('task')}")
    values.extend(source)
    return {str(value).strip().lower() for value in values if value}


def _null_signature(expr):
    """Power-law slopes invariant to variable renaming and scalar prefactors."""
    # Freeze decimal exponents exactly before simplification. Floating exponent
    # subtraction can otherwise leave spurious variable powers near 1e-16.
    expr = expr.xreplace({value: sp.Rational(str(value)) for value in expr.atoms(sp.Float)})
    exponents = []
    for variable in expr.free_symbols:
        exponent = sp.simplify(variable * sp.diff(expr, variable) / expr)
        if exponent.free_symbols:
            raise ValueError(f"H0 expression is not a monomial: {expr}")
        exponents.append(sp.srepr(exponent))
    if not exponents:
        raise ValueError("A constant expression is not a null regression task")
    return tuple(sorted(exponents))


def _cross_role_collisions(records, key):
    by_key = defaultdict(list)
    for record in records:
        by_key[record[key]].append(record)
    collisions = []
    for value, members in by_key.items():
        roles = sorted({m["split"] for m in members})
        if len(roles) > 1:
            collisions.append({"key": str(value), "splits": roles,
                               "tasks": sorted(m["task"] for m in members)})
    return collisions


def assert_group_isolation(tasks, *, h0_structure_exception=None):
    """Validate registry isolation before data generation; return audit evidence.

    ``h0_structure_exception`` must be an explicit nonempty justification to
    permit identical frozen null slope signatures across roles. Sharing the
    general monomial family is always permitted and explicitly reported.
    """
    tasks = list(tasks)
    if not tasks:
        raise ValueError("Cannot certify an empty registry")
    names = [t.get("task") for t in tasks]
    if any(not isinstance(n, str) or not n.strip() for n in names):
        raise ValueError("Every task needs a nonempty string identity")
    if len(set(names)) != len(names):
        raise ValueError("Duplicate task identities in frozen registry")
    if h0_structure_exception is not None and (
        not isinstance(h0_structure_exception, str) or not h0_structure_exception.strip()
    ):
        raise ValueError("Null structure exception needs an explicit justification")

    inventory = _historical_inventory()
    alias_map = {}
    old_canonical, old_core = set(), set()
    for source in inventory["equations"]:
        alias_map[source["eq_name"].lower()] = source
        alias_map[source["class_name"].lower()] = source
        if source["historical32"] and source["historical_class"] == "H1":
            old_canonical.add(source["canonical_structure"])
            old_core.add(source["nonlinear_core"])

    normalized, old_seen, tainted_seen = [], [], []
    for task in tasks:
        split = task.get("split")
        if not isinstance(split, str) or not split.strip():
            raise ValueError(f"Missing split for {task['task']}")
        unit = group_unit(task)
        old_aliases, tainted_aliases = [], []
        for alias in _aliases(task):
            source = alias_map.get(alias)
            if source and source["historical32"]:
                old_aliases.append(source["eq_name"])
            if source and source["historically_tainted_group"]:
                tainted_aliases.append(source["eq_name"])
        manually_tainted = any(task.get(key) is True for key in (
            "historically_tainted", "historical_exposure", "development_only"))
        if (old_aliases or tainted_aliases or manually_tainted) and split != "development":
            raise ValueError(f"Historical exposure requires development split: {task['task']}")
        old_seen.extend(old_aliases)
        tainted_seen.extend(tainted_aliases)
        expr = _parse_expression(task)
        record = {"task": task["task"], "split": split, "klass": task["klass"], "unit": unit}
        if task["klass"] == "H0":
            record["null_signature"] = _null_signature(expr)
        else:
            record["canonical"] = canonical_structure(expr)
            record["nonlinear_core"] = nonlinear_core(expr)
            if (record["canonical"] in old_canonical or record["nonlinear_core"] in old_core) and split != "development":
                raise ValueError(f"Historical formula/core clone requires development split: {task['task']}")
        normalized.append(record)

    unit_collisions = _cross_role_collisions(normalized, "unit")
    if unit_collisions:
        raise ValueError(f"A structural/design group crosses splits: {unit_collisions}")
    h1 = [r for r in normalized if r["klass"] == "H1"]
    syntax_collisions = _cross_role_collisions(h1, "canonical")
    core_collisions = _cross_role_collisions(h1, "nonlinear_core")
    if syntax_collisions or core_collisions:
        raise ValueError(f"H1 syntactic/core twins cross splits: {syntax_collisions + core_collisions}")
    null_collisions = _cross_role_collisions([r for r in normalized if r["klass"] == "H0"], "null_signature")
    if null_collisions and h0_structure_exception is None:
        raise ValueError(f"Exact null slope clones cross roles without explicit exception: {null_collisions}")
    by_split = defaultdict(lambda: {"H0": 0, "H1": 0, "units": set()})
    for record in normalized:
        by_split[record["split"]][record["klass"]] += 1
        by_split[record["split"]]["units"].add(record["unit"])
    return {
        "passed": True, "n_tasks": len(tasks),
        "class_counts": dict(Counter(t["klass"] for t in tasks)),
        "splits": {name: {"H0": counts["H0"], "H1": counts["H1"],
                          "n_units": len(counts["units"]), "units": sorted(counts["units"])}
                   for name, counts in sorted(by_split.items())},
        "historical32_identities_included": sorted(set(old_seen)),
        "historically_tainted_source_aliases": sorted(set(tainted_seen)),
        "h1_canonical_cross_split_collisions": syntax_collisions,
        "h1_nonlinear_core_cross_split_collisions": core_collisions,
        "h0_exact_slope_cross_split_collisions": null_collisions,
        "h0_function_family_exception": "H0 tasks intentionally share the monomial/log-linear hypothesis family; design_group is their independent unit.",
        "h0_exact_structure_exception": h0_structure_exception,
        "limitations": [
            "Canonical syntax/nonlinear cores are conservative checks, not complete algebraic equivalence tests.",
            "Manual structural_group assignments are primary and may merge additional source or mechanism near-twins.",
            "Parameter-window, coordinate-scale or prefactor changes do not create an independent structural family.",
        ],
    }
