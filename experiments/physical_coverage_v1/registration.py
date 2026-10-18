"""Build a120-equation physical-domain census without samples or fitted results.

Allowed inputs are pinned source/structure metadata and historical task identity
manifests only. This module deliberately has no sampling, fitting, gate or result
analysis API. Its output is a metadata review package, not an execution freeze.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import mpmath as mp
import sympy as sp

from experiments.fixed_domain.inventory import canonical_structure

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
INVENTORY = ROOT / "runs/p02_inventory/inventory.json"
P02_TASKS = ROOT / "runs/p02_fixed_domain/task_registry.json"
PROTOCOL = PACKAGE / "protocol.json"

# Deliberately conservative mechanisms, declared from formula/source inspection.
# Generic shared operators alone are not a claim of independence. In particular,
# powers/reciprocals and positive sums are not promoted as new mechanisms merely
# by changing the output power, number of terms or monomial feature products.
FAMILY_MERGES = {
    "gaussian_activation": [
        "feynman-i.6.20", "feynman-i.6.20a", "feynman-i.6.20b", "feynman-i.40.1",
    ],
    "relativistic_or_offset_radical": [
        "feynman-i.10.7", "feynman-i.15.10", "feynman-i.15.3t", "feynman-i.15.3x",
        "feynman-i.16.6", "feynman-i.34.10", "feynman-i.34.14", "feynman-i.48.2",
        "feynman-ii.13.23", "feynman-ii.13.34", "feynman-ii.21.32", "feynman-ii.24.17",
        "feynman-bonus.6", "feynman-bonus.10", "feynman-bonus.15",
    ],
    "positive_sum_or_norm": [
        "feynman-i.11.19", "feynman-i.13.4", "feynman-i.24.6",
        "feynman-ii.6.15a", "feynman-ii.36.38", "feynman-iii.10.19",
        "feynman-bonus.5", "feynman-bonus.17", "feynman-bonus.18",
    ],
    "affine_or_additive_rational": [
        "feynman-i.18.4", "feynman-i.27.6", "feynman-i.39.11", "feynman-i.43.43",
        "feynman-ii.10.9", "feynman-ii.11.27", "feynman-ii.11.28",
        "feynman-ii.37.1", "feynman-ii.38.14",
    ],
    "bose_exponential_difference": [
        "feynman-i.41.16", "feynman-iii.4.32", "feynman-iii.4.33", "feynman-iii.14.14",
    ],
    "hyperbolic_response": ["feynman-ii.35.18", "feynman-ii.35.21"],
    "positive_log": ["feynman-i.44.4"],
    "geometric_separation": ["feynman-i.8.14", "feynman-i.9.18"],
    "frequency_difference_resonance": ["feynman-i.32.17", "feynman-ii.11.3"],
    "positive_source_difference": ["feynman-i.13.12", "feynman-ii.2.42"],
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _write_json(path, payload):
    Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def _positive_support(sampling):
    """Condition original permitted sign branches; never invent a new scale."""
    attributes, kind = sampling["attributes"], sampling["class_name"]
    lo, hi = attributes["min_value"], attributes["max_value"]
    positive = attributes.get("uses_positive", True)
    negative = attributes.get("uses_negative", True)
    if kind == "IntegerSampling":
        raw = set(range(int(lo), int(hi)))  # upstream randint is high-exclusive
        support = set()
        if positive:
            support.update(raw)
        if negative:
            support.update(-value for value in raw)
        support = sorted(value for value in support if value > 0)
        if len(support) < 2:
            raise ValueError("positive_integer_branch_has_fewer_than2_values")
        a = support[max(0, math.ceil(.05 * len(support)) - 1)]
        b = support[max(0, math.ceil(.95 * len(support)) - 1)]
        return [int(a), int(b)], "integer", {"original_positive_integer_support": support}
    if kind == "DefaultSampling":
        if not positive:
            raise ValueError("no_original_positive_input_support")
        if lo <= 0 or hi <= lo:
            raise ValueError("positive_log_support_not_established")
        bounds = [math.exp(math.log(lo) + q * (math.log(hi) - math.log(lo)))
                  for q in (.05, .95)]
        return bounds, "log", {"original_positive_support": [lo, hi]}
    if kind != "SimpleSampling":
        raise ValueError(f"unsupported_source_sampling_class:{kind}")
    intervals = []
    if positive and hi > 0:
        intervals.append([max(float(lo), 0.), float(hi)])
    if negative and lo < 0:
        intervals.append([max(float(-hi), 0.), float(-lo)])
    if not intervals:
        raise ValueError("no_original_positive_continuous_branch")
    # All source intervals crossing zero are symmetric or have a single
    # permitted branch. Do not approximate an unsupported disconnected mixture.
    merged = [min(a for a, _ in intervals), max(b for _, b in intervals)]
    if any(interval != merged for interval in intervals):
        raise ValueError("positive_branch_mixture_requires_exact_cdf_review")
    a, b = merged
    return [a + .05 * (b - a), a + .95 * (b - a)], "linear", {
        "original_positive_support": merged,
    }


def _bounds(value):
    return mp.mpf(value._mpi_[0]), mp.mpf(value._mpi_[1])


def _as_interval_number(number):
    if number == sp.pi:
        return mp.iv.pi
    if isinstance(number, (sp.Integer, sp.Rational, sp.Float)):
        rational = sp.Rational(number)
        return mp.iv.mpf(str(rational.p)) / mp.iv.mpf(str(rational.q))
    raise ValueError(f"unsupported_numeric_atom:{number}")


def interval_value(expression, variables):
    """Outward interval evaluation of a whole box; no point evaluation."""
    if isinstance(expression, sp.Symbol):
        return variables[str(expression)]
    if not expression.free_symbols and (expression == sp.pi or isinstance(expression, (sp.Integer, sp.Rational, sp.Float))):
        return _as_interval_number(expression)
    if isinstance(expression, sp.Add):
        result = mp.iv.mpf(0)
        for argument in expression.args:
            result += interval_value(argument, variables)
        return result
    if isinstance(expression, sp.Mul):
        result = mp.iv.mpf(1)
        for argument in expression.args:
            result *= interval_value(argument, variables)
        return result
    if isinstance(expression, sp.Pow):
        base = interval_value(expression.base, variables)
        lower, upper = _bounds(base)
        exponent = expression.exp
        if exponent.free_symbols:
            raise ValueError("variable_exponent_requires_proof")
        if exponent.is_integer:
            if int(exponent) < 0 and lower <= 0 <= upper:
                raise ValueError("denominator_interval_contains_zero")
            return base ** int(exponent)
        if lower < 0:
            raise ValueError("fractional_power_base_not_proven_nonnegative")
        if lower == 0 and exponent < 0:
            raise ValueError("negative_power_base_can_be_zero")
        if exponent == sp.Rational(1, 2):
            return mp.iv.sqrt(base)
        return base ** interval_value(exponent, variables)
    if expression.func in {sp.exp, sp.log, sp.sin, sp.cos, sp.tanh}:
        argument = interval_value(expression.args[0], variables)
        if expression.func == sp.log:
            if _bounds(argument)[0] <= 0:
                raise ValueError("log_argument_not_proven_positive")
            return mp.iv.ln(argument)
        if expression.func == sp.tanh:
            exponential = mp.iv.exp(2 * argument)
            return (exponential - 1) / (exponential + 1)
        function = {sp.exp: mp.iv.exp, sp.sin: mp.iv.sin, sp.cos: mp.iv.cos}[expression.func]
        return function(argument)
    raise ValueError(f"interval_operation_not_supported:{expression.func}")


def _safe_identity(record, expression):
    """Two general positive identities; never change a task's source law."""
    source_expression = sp.sympify(record["sympy_srepr"])
    variables = sorted(source_expression.free_symbols, key=lambda x: int(str(x)[1:]))
    rewrites = []
    if record["eq_name"] == "feynman-bonus.20":
        a, b, angle = variables
        original = a / b + b / a - sp.sin(angle) ** 2
        identity = (a - b) ** 2 / (a * b) + 1 + sp.cos(angle) ** 2
        if sp.trigsimp(sp.cancel(original - identity)) != 0:
            raise AssertionError("Invalid registered positive-ratio identity")
        expression = expression.xreplace({original: identity})
        rewrites.append("r+1/r-sin(theta)^2=(a-b)^2/(a*b)+1+cos(theta)^2>=1")
    if record["eq_name"] == "feynman-bonus.13":
        _, _, a, b, angle = variables
        original = a ** 2 - 2 * a * b * sp.cos(angle) + b ** 2
        identity = (a - b) ** 2 + 2 * a * b * (1 - sp.cos(angle))
        if sp.expand(original - identity) != 0:
            raise AssertionError("Invalid geometric-separation identity")
        expression = expression.xreplace({original: identity})
        rewrites.append("r^2+d^2-2*r*d*cos(theta)=(r-d)^2+2*r*d*(1-cos(theta))")
    return expression, rewrites


def _source_constraints(record, variables):
    doc = record.get("upstream_docstring") or ""
    tail = doc.split("- Constraints:", 1)[1] if "- Constraints:" in doc else ""
    constraints = [line.strip()[2:] for line in tail.splitlines() if line.strip().startswith("- ")]
    checks, obligations, warnings = [], [], []
    local = {name: sp.Symbol(name, real=True) for name in variables}
    local.update({"sin": sp.sin, "cos": sp.cos, "pi": sp.pi})
    for raw in constraints:
        indices = [int(value) for value in re.findall(r"x\[(\d+)\]", raw)]
        if indices and max(indices) >= record["n_vars"]:
            reason = "source_constraint_index_out_of_declared_dimension"
            warnings.append({"constraint": raw, "reason": reason})
            obligations.append(f"Resolve source constraint indexing without inventing a law: {raw}")
            checks.append({"constraint": raw, "certified": False, "reason": reason})
            continue
        rewritten = re.sub(r"x\[(\d+)\]", r"x\1", raw).replace("np.", "")
        parts = re.split(r"\s*(!=|>=|<=|>|<)\s*", rewritten, maxsplit=1)
        if len(parts) != 3:
            obligations.append(f"Translate and prove source constraint: {raw}")
            checks.append({"constraint": raw, "certified": False, "reason": "unparsed_constraint"})
            continue
        left, operation, right = parts
        try:
            difference = sp.sympify(left, locals=local) - sp.sympify(right, locals=local)
            difference, identities = _safe_identity(record, difference)
            lo, hi = _bounds(interval_value(difference, variables))
            certified = {"!=": lo > 0 or hi < 0, ">": lo > 0, ">=": lo >= 0,
                         "<": hi < 0, "<=": hi <= 0}[operation]
            checks.append({"constraint": raw, "certified": bool(certified),
                           "difference_interval": [mp.nstr(lo, 35), mp.nstr(hi, 35)],
                           "symbolically_checked_identities": identities})
            if not certified:
                obligations.append(f"Prove this source constraint on the entire unchanged box: {raw}")
        except (ValueError, ZeroDivisionError, TypeError, KeyError) as exc:
            checks.append({"constraint": raw, "certified": False, "reason": str(exc)})
            obligations.append(f"Resolve/prove source constraint over the whole unchanged box: {raw}")
    # Named source velocities are physical assumptions in addition to symbolic
    # denominator validity; source coordinate labels and units remain untouched.
    for line in doc.splitlines():
        match = re.match(r"\s*- x\[(\d+)\]:\s+(v|u|w)\s+\(", line)
        if match:
            index = int(match.group(1))
            if index < record["n_vars"]:
                lo, hi = _bounds(variables[f"x{index}"])
                certified = lo > 0 and hi < mp.mpf("299792458")
                checks.append({"constraint": f"source physical velocity {match.group(2)}:0<v<c",
                               "certified": bool(certified)})
                if not certified:
                    obligations.append(f"Prove source velocity x{index} remains subluminal on the full box")
    return checks, obligations, warnings


def _source_variable_type_review(record):
    checks = []
    for description in record["variable_descriptions"]:
        match = re.match(r"x\[(\d+)\]:.*?\(([^)]*)\)", description)
        if not match or "integer" not in match.group(2).lower():
            continue
        index = int(match.group(1))
        if index >= record["n_vars"]:
            checks.append({"description": description, "status": "source_type_review_pending",
                           "reason": "variable_description_index_outside_source_dimension"})
            continue
        sampling = record["sampling_objects"][index]["class_name"]
        relaxation = "real" in match.group(2).lower()
        status = ("source_integer_sampling_consistent" if sampling == "IntegerSampling"
                  else "source_documented_integer_to_real_relaxation" if relaxation
                  else "source_type_review_pending")
        checks.append({"variable_index": index, "description": description,
                       "original_sampling_class": sampling, "status": status,
                       "relaxation_evidence": description if relaxation else None})
    return checks


def _training_boxes(outer, scales, widths):
    result, obligations = {}, []
    for width in widths:
        boxes = []
        for j, ((lo, hi), scale) in enumerate(zip(outer, scales, strict=True)):
            a, b = (math.log(lo), math.log(hi)) if scale == "log" else (lo, hi)
            center = (a + b) / 2
            lower, upper = center - width * (b - a) / 2, center + width * (b - a) / 2
            if scale == "log":
                bounds = [math.exp(lower), math.exp(upper)]
            elif scale == "integer":
                bounds = [math.ceil(lower), math.floor(upper)]
                if bounds[1] - bounds[0] < 1:
                    obligations.append(f"width={width}:integer_x{j}_has_fewer_than2_training_values")
                if bounds[0] <= lo and bounds[1] >= hi:
                    obligations.append(f"width={width}:integer_x{j}_has_no_external_support")
            else:
                bounds = [lower, upper]
            boxes.append(bounds)
        result[str(width)] = boxes
    return result, obligations


def _monomial(expression):
    if not expression.free_symbols or expression.is_Symbol:
        return True
    if isinstance(expression, sp.Mul):
        return all(_monomial(value) for value in expression.args)
    return (isinstance(expression, sp.Pow) and not expression.exp.free_symbols
            and _monomial(expression.base))


def _power_core(expression):
    expression = sp.factor_terms(expression)
    if _monomial(expression):
        return sp.Integer(1)
    if isinstance(expression, sp.Mul):
        arguments = [argument for argument in expression.args if not _monomial(argument)]
        if len(arguments) == 1:
            return _power_core(arguments[0])
        return sp.Mul(*[_power_core(argument) for argument in arguments])
    if isinstance(expression, sp.Pow) and not expression.exp.free_symbols and expression.exp != 0:
        return _power_core(expression.base)
    return expression


def _positive_laurent_sum_tag(expression):
    core = _power_core(expression)
    expanded = sp.expand(core)
    terms = sp.Add.make_args(expanded)
    variables = sorted(expanded.free_symbols, key=str)
    eligible = len(terms) >= 2
    for term in terms:
        coefficient, variable_part = term.as_independent(*variables, as_Add=False)
        eligible = eligible and coefficient.is_positive is True and _monomial(variable_part)
    return {
        "matches": bool(eligible),
        "grammar": "one_whole_output_power_of_positive_coefficient_Laurent_sum_after_monomial_prefactors",
        "core_expression": str(core), "expanded_core": str(expanded),
        "boundary": "nested nonlinear summands and nonpositive expanded coefficients are outside this fine grammar",
    }


def _nested_additive_radical_tag(expression):
    for node in sp.preorder_traversal(expression):
        if not isinstance(node, sp.Add):
            continue
        has_other_monomial = any(_monomial(argument) for argument in node.args)
        has_nonlinear_radical = any(
            isinstance(descendant, sp.Pow) and not descendant.exp.free_symbols
            and descendant.exp.is_integer is not True and not _monomial(descendant.base)
            for argument in node.args for descendant in sp.preorder_traversal(argument)
        )
        if has_other_monomial and has_nonlinear_radical:
            return True
    return False


def _group_records(records, p02):
    parent = {record["eq_name"]: record["eq_name"] for record in records}
    edges = []

    def find(value):
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def unite(members, reason):
        for other in members[1:]:
            if find(members[0]) != find(other):
                parent[find(other)] = find(members[0])
            edges.append({"a": members[0], "b": other, "reason": reason})

    for key in ("canonical_structure", "constant_masked_structure", "nonlinear_core", "source_family", "power_stripped_core"):
        by_key = defaultdict(list)
        for record in records:
            by_key[record[key]].append(record["eq_name"])
        for members in by_key.values():
            unite(members, key)
    for family, members in FAMILY_MERGES.items():
        unite(members, f"manual_conservative_family:{family}")
    oscillatory = [record["eq_name"] for record in records
                   if set(record["unary_functions"]) & {"sin", "cos"}]
    unite(oscillatory, "manual_conservative_family:all_source_angular_interference")
    # All prior manual groups are authoritative evidence of exposure. A source
    # retaining a different scalar prefactor cannot form a new physical group.
    source_p02 = [task for task in p02 if task["family"] == "upstream_source"]
    by_p02 = defaultdict(list)
    for task in source_p02:
        by_p02[task["structural_group"]].append(task["task"])
    for label, members in by_p02.items():
        unite(members, f"prior_P02_group:{label}")
    direct = {task["task"]: task for task in source_p02}
    positive_sums = [record["eq_name"] for record in records if record["positive_laurent_sum_tag"]["matches"]]
    unite(positive_sums, "automatic_all120_positive_Laurent_sum_grammar")
    fine_components = defaultdict(list)
    for record in records:
        fine_components[find(record["eq_name"])].append(record)
    for index, members in enumerate(sorted(fine_components.values(), key=lambda values: min(value["eq_name"] for value in values))):
        prior = [record["eq_name"] for record in members if record["eq_name"] in direct]
        mechanism_overlap = any(record["eq_name"] in {"feynman-i.32.17", "feynman-ii.11.3"} for record in members)
        for record in members:
            record["fine_structural_group"] = f"fine-structure-{index:03d}"
            record["fine_exposure"] = "development" if prior or mechanism_overlap else "prospective_candidate"
    nested_radical = [record["eq_name"] for record in records if record["nested_additive_radical_tag"]]
    unite(positive_sums + nested_radical, "primary_conservative_sum_norm_nested_combination;fine_structure_assignment_disputed_not_algebraic_equality")
    components = defaultdict(list)
    for record in records:
        components[find(record["eq_name"])].append(record)
    groups = []
    for index, members in enumerate(sorted(components.values(), key=lambda values: min(v["eq_name"] for v in values))):
        group_id = f"physical-family-{index:03d}"
        prior = sorted(record["eq_name"] for record in members if record["eq_name"] in direct)
        designed_exposure = []
        if any(record["eq_name"] in {"feynman-i.32.17", "feynman-ii.11.3"} for record in members):
            designed_exposure = [{
                "mechanism": "P02 designed damped_resonance",
                "registered_role": "confirmation", "exposure_level": "prior_study_outcomes",
                "evidence": "prior P02 registry identity/role only; no performance read",
                "reason": "source frequency-difference resonance is a nearby mechanism",
            }]
        exposed = bool(prior or designed_exposure)
        group_level = ("prior_study_outcomes" if prior or any(item["exposure_level"] == "prior_study_outcomes" for item in designed_exposure)
                       else "mechanism_design_only" if designed_exposure else "metadata_only")
        for record in members:
            record["structural_group"] = group_id
            conservative_only = exposed and record["fine_exposure"] == "prospective_candidate"
            record["exposure"] = ("development_conservative" if conservative_only else "development") if exposed else "prospective_candidate"
            record["structural_assignment_status"] = "candidate_pending_structural_review" if conservative_only else "conservative_group_assignment_pending_full_catalog_review"
            record["source_exposure_level"] = "prior_study_outcomes" if record["eq_name"] in direct else "metadata_only"
            record["group_exposure_level"] = group_level
            record["all120_source_metadata_previously_seen"] = True
            record["exposure_evidence"] = {
                "direct_source_exposure": record["eq_name"] in direct,
                "historical32": record["historical32"],
                "P02_additional_source": record["eq_name"] in direct and not record["historical32"],
                "group_prior_source_members": prior, "prior_designed_mechanisms": designed_exposure,
                "primary_superfamily_inheritance_only": conservative_only,
                "interpretation": "Related-family exposure is conservative holdout protection; it is not evidence this individual source was fitted or an algebraic equality proof.",
            }
        groups.append({
            "structural_group": group_id, "members": sorted(record["eq_name"] for record in members),
            "exposure": "development" if exposed else "prospective_candidate",
            "group_exposure_level": group_level,
            "prior_source_members": prior, "prior_designed_mechanisms": designed_exposure,
            "scientifically_eligible_members": sorted(record["eq_name"] for record in members if record["scientific_status"] == "eligible_analytic"),
            "H1_eligible_members": sorted(record["eq_name"] for record in members if record["scientific_status"] == "eligible_analytic" and record["klass"] == "H1"),
        })
    return groups, edges


def build_registration():
    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    p02 = json.loads(P02_TASKS.read_text(encoding="utf-8"))
    source_path = Path(inventory["upstream_source"])
    if digest(source_path) != inventory["upstream_sha256"]:
        raise ValueError("Pinned source changed; require a new source census")
    mp.mp.dps, mp.iv.dps = 100, 80
    records = []
    for source in inventory["equations"]:
        record = {key: source[key] for key in (
            "eq_name", "class_name", "n_vars", "expression", "sympy_srepr",
            "upstream_docstring", "upstream_source_line", "sampling_objects",
            "historical32", "historical_class", "source_family", "canonical_structure",
            "constant_masked_structure", "nonlinear_core", "unary_functions",
        )}
        expression = sp.sympify(source["sympy_srepr"])
        record["klass"] = "H0" if source["strict_monomial_structure"] else "H1"
        record["power_stripped_core"] = canonical_structure(_power_core(expression), ignore_numbers=True)
        record["positive_laurent_sum_tag"] = _positive_laurent_sum_tag(expression)
        record["nested_additive_radical_tag"] = _nested_additive_radical_tag(expression)
        record["unit_evidence"] = "evidence_missing_explicit_units; preserve source descriptions, original numeric coordinates and constants"
        record["variable_descriptions"] = [line.strip()[2:] for line in (record["upstream_docstring"] or "").splitlines()
                                           if re.match(r"\s*- x\[\d+\]:", line)]
        record["source_variable_type_review"] = _source_variable_type_review(record)
        record["physical_source_contract_review_status"] = "pending_independent_physical_semantics_types_and_units_review"
        record["scientific_status"] = "exclusion_unresolved"
        record["inclusion"] = False
        record["inclusion_reason"] = None
        record["exclusion_reason"] = None
        record["proof_obligations"] = []
        record["execution_gates"] = []
        outer, scales, support = [], [], []
        try:
            for sampling in record["sampling_objects"]:
                bounds, scale, evidence = _positive_support(sampling)
                outer.append(bounds)
                scales.append(scale)
                support.append(evidence)
            record["outer_bounds"] = outer
            record["coordinate_scales"] = scales
            record["positive_source_support_evidence"] = support
            training, design_obligations = _training_boxes(outer, scales, protocol["common_factorial_design"]["training_width_fractions"])
            record["train_bounds_by_width"] = training
            record["proof_obligations"].extend(design_obligations)
            variables = {f"x{j}": mp.iv.mpf([str(lo), str(hi)]) for j, (lo, hi) in enumerate(outer)}
            checks, source_obligations, warnings = _source_constraints(record, variables)
            record["source_constraint_checks"] = checks
            record["source_metadata_warnings"] = warnings
            record["proof_obligations"].extend(source_obligations)
            checked_expression, identities = _safe_identity(record, expression)
            record["analytic_identities"] = identities
            try:
                lo, hi = _bounds(interval_value(checked_expression, variables))
                record["truth_interval_certificate"] = {
                    "lower": mp.nstr(lo, 45), "upper": mp.nstr(hi, 45),
                    "method": "80-decimal outward interval evaluation of the full compact box; optional symbolically checked identities",
                    "random_or_grid_point_evaluations": 0,
                }
                if not mp.isfinite(lo) or not mp.isfinite(hi):
                    record["proof_obligations"].append("Prove a real finite mathematical target on the complete unchanged box")
                elif hi <= 0:
                    record["scientific_status"] = "excluded_nonpositive_target"
                    record["exclusion_reason"] = "signed_source_target_nonpositive_on_entire_positive_input_box"
                elif lo <= 0:
                    record["proof_obligations"].append("Prove strictly positive target everywhere on the complete unchanged box; interval lower bound<=0 is not a sign certificate")
                elif not record["proof_obligations"]:
                    record["scientific_status"] = "eligible_analytic"
                    record["inclusion"] = True
                    record["inclusion_reason"] = "whole_box_positive_finite_and_source_constraints_certified"
                if lo > 0 and (lo < mp.mpf("2.2250738585072014e-308") or hi > mp.mpf("1.7976931348623157e308")):
                    record["execution_gates"].append("stable_log_truth_adapter_pending")
            except (ValueError, ZeroDivisionError, TypeError, KeyError) as exc:
                record["truth_interval_certificate"] = {"status": "not_certified", "reason": str(exc), "random_or_grid_point_evaluations": 0}
                record["proof_obligations"].append(f"Whole-box target certificate unresolved:{exc}")
        except (ValueError, OverflowError) as exc:
            record["outer_bounds"] = outer if len(outer) == record["n_vars"] else None
            record["coordinate_scales"] = scales if len(scales) == record["n_vars"] else None
            record["proof_obligations"].append(f"Original positive support/design unresolved:{exc}")
            if str(exc) == "no_original_positive_input_support":
                record["scientific_status"] = "excluded_no_positive_input_support"
                record["exclusion_reason"] = "at_least_one_original_source_feature_has_only_negative_support;no_abs_or_sign_change_permitted"
        if record["n_vars"] > protocol["source_census"]["maximum_features_scientific"]:
            record["scientific_status"] = "excluded_dimension_limit"
            record["inclusion"] = False
            record["exclusion_reason"] = "source_dimension_above_registered_maximum8"
        if record["scientific_status"] == "exclusion_unresolved":
            record["exclusion_reason"] = "; ".join(record["proof_obligations"]) or "whole_box_certificate_not_established"
        if record["n_vars"] > 4:
            record["execution_gates"].append("diagnostic_and_grid_adapter_pending")
        if any(obj["class_name"] == "IntegerSampling" for obj in record["sampling_objects"]):
            record["execution_gates"].append("original_integer_semantics_adapter_pending")
        if record.get("source_metadata_warnings"):
            record["execution_gates"].append("source_contract_review_pending")
        if any(item["status"] == "source_type_review_pending" for item in record["source_variable_type_review"]):
            record["execution_gates"].append("source_type_review_pending")
        record["execution_gates"].append("independent_physical_source_contract_review_pending")
        if "exp" in record["unary_functions"]:
            record["execution_gates"].append("stable_exponential_log_truth_audit_pending")
        if "cos" in record["unary_functions"]:
            record["execution_gates"].append("stable_trigonometric_cancellation_audit_pending")
        record["execution_gates"].append("common_grid_null_controls_and_seed_freeze_pending")
        record["execution_status"] = "await_adapter_and_registration" if record["inclusion"] else "not_admitted_pending_domain_or_sign"
        record["review_status"] = "metadata_census_complete_independent_review_pending"
        record["analytic_math_status"] = record["scientific_status"]
        records.append(record)
    if len(records) != 120 or len({record["eq_name"] for record in records}) != 120:
        raise AssertionError("Census must contain all120 unique source formulas")
    groups, edges = _group_records(records, p02)
    summary = {
        "all_equations": len(records),
        "scientific_status": dict(Counter(record["scientific_status"] for record in records)),
        "class_counts": dict(Counter(record["klass"] for record in records)),
        "all_groups": len(groups),
        "eligible_source_equations": sum(record["inclusion"] for record in records),
        "eligible_H0_source_equations": sum(record["inclusion"] and record["klass"] == "H0" for record in records),
        "eligible_H1_source_equations": sum(record["inclusion"] and record["klass"] == "H1" for record in records),
        "eligible_H1_development_source_equations": sum(record["inclusion"] and record["klass"] == "H1" and record["exposure"].startswith("development") for record in records),
        "eligible_H1_prospective_source_equations": sum(record["inclusion"] and record["klass"] == "H1" and record["exposure"] == "prospective_candidate" for record in records),
        "eligible_H1_development_groups": sum(bool(group["H1_eligible_members"]) and group["exposure"] == "development" for group in groups),
        "eligible_H1_prospective_groups": sum(bool(group["H1_eligible_members"]) and group["exposure"] == "prospective_candidate" for group in groups),
        "eligible_H1_disputed_fine_structure_candidates": sum(record["inclusion"] and record["klass"] == "H1" and record["structural_assignment_status"] == "candidate_pending_structural_review" for record in records),
        "fine_structure_groups_before_primary_superfamily_merge": len({record["fine_structural_group"] for record in records}),
        "eligible_dimension_above4": sum(record["inclusion"] and record["n_vars"] > 4 for record in records),
        "source_type_review_pending_equations": sum("source_type_review_pending" in record["execution_gates"] for record in records),
        "study_samples_generated": 0, "study_fits_generated": 0, "diagnostic_outcomes_generated": 0,
        "registration_status": "pre_execution_registration_incomplete",
    }
    return protocol, records, groups, edges, summary


def write_registration(out):
    out = Path(out).resolve()
    if out.parent != (ROOT / "runs").resolve() or not out.name.startswith("p03_physical_protocol"):
        raise ValueError("Metadata output must use a new runs/p03_physical_protocol* directory")
    if (out / "review_manifest.json").exists():
        verify_registration(out)
        return json.loads((out / "review_manifest.json").read_text(encoding="utf-8"))["summary"]
    if out.exists() and any(path.is_file() for path in out.iterdir()):
        raise ValueError("Nonempty metadata output without manifest; require a new version/directory")
    out.mkdir(parents=True, exist_ok=True)
    protocol, records, groups, edges, summary = build_registration()
    _write_json(out / "protocol_rules.json", protocol)
    _write_json(out / "candidates.json", {"status": "pre_execution_registration_incomplete", "equations": records})
    _write_json(out / "groups.json", {"groups": groups, "edges": edges})
    p02 = json.loads(P02_TASKS.read_text(encoding="utf-8"))
    mechanisms = {}
    for task in p02:
        if task["source_kind"] != "synthetic_positive_mechanism":
            continue
        mechanisms[task["family"]] = {
            "family": task["family"], "registered_role": task["split"],
            "exposure_level": "mechanism_design_only" if task["split"] == "reserve" else "prior_study_outcomes",
            "evidence_scope": "identity and registered role, not read diagnostic/prediction outcomes",
        }
    _write_json(out / "exposure_ledger.json", {
        "all120_minimum_exposure": "metadata_only",
        "prospective_definition": "absence of prior study outcomes/design overlap; not an unseen-formula claim",
        "P02_physical_sources": [task["task"] for task in p02 if task["family"] == "upstream_source"],
        "P02_designed_mechanisms": list(mechanisms.values()),
    })
    _write_json(out / "summary.json", summary)
    fields = ["eq_name", "class_name", "n_vars", "klass", "scientific_status", "inclusion",
              "inclusion_reason", "exclusion_reason", "structural_group", "exposure",
              "source_exposure_level", "group_exposure_level",
              "analytic_math_status", "physical_source_contract_review_status", "source_variable_type_review",
              "fine_structural_group", "fine_exposure", "structural_assignment_status", "positive_laurent_sum_tag",
              "execution_status", "expression", "outer_bounds", "train_bounds_by_width",
              "coordinate_scales", "proof_obligations", "execution_gates", "unit_evidence"]
    with (out / "candidates.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({field: json.dumps(record.get(field), ensure_ascii=False)
                             if isinstance(record.get(field), (list, dict)) else record.get(field)
                             for field in fields})
    lines = ["# Physical coverage source census — metadata only", "",
             "Status: pre_execution_registration_incomplete. No study samples, fits or diagnostics.", "",
             "```json", json.dumps(summary, ensure_ascii=False, indent=2), "```", "",
             "| Source equation | p | H0/H1 | Scientific decision | Exposure | Group |",
             "|---|---:|---|---|---|---|"]
    for record in records:
        lines.append(f"| {record['eq_name']} | {record['n_vars']} | {record['klass']} | {record['scientific_status']} | {record['exposure']} | {record['structural_group']} |")
    lines += ["", "Every exclusion, analytic bound, source constraint, adapter obligation and complete derived box is retained in candidates.json/.csv."]
    (out / "REGISTRATION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    inputs = {str(path): digest(path) for path in (INVENTORY, P02_TASKS, PROTOCOL,
              *sorted(PACKAGE.glob("*.py")), ROOT / "experiments/fixed_domain/inventory.py")}
    upstream = json.loads(INVENTORY.read_text(encoding="utf-8"))
    inputs[upstream["upstream_source"]] = upstream["upstream_sha256"]
    sampling_path = Path(upstream["upstream_source"]).with_name("sampling.py")
    inputs[str(sampling_path)] = digest(sampling_path)
    # Independent audits authenticate this manifest and must never be its own
    # transitive hash input; the archived pre-guard review is likewise separate.
    metadata_names = ["protocol_rules.json", "candidates.json", "groups.json", "summary.json",
                      "exposure_ledger.json", "candidates.csv", "REGISTRATION.md"]
    files = {name: digest(out / name) for name in metadata_names}
    manifest = {
        "registration_status": "pre_execution_registration_incomplete",
        "rules_status": "outcome_independent_rule_snapshot_hashed_for_review",
        "rules_sha256": digest(PROTOCOL), "candidate_status": "complete120_source_census_with_explicit_unresolved_rows",
        "independent_review_status": "pending", "source_and_input_sha256": inputs,
        "metadata_artifact_sha256": files, "summary": summary,
        "forbidden_performance_inputs_read": [],
        "allowed_inputs_only": "source formulas/constants/sampling/constraints; structure inventory; P02 identity/group manifest; local registration sources",
        "not_an_execution_freeze": True,
        "proof_environment": {"sympy": sp.__version__, "mpmath": mp.__version__,
                              "interval_decimal_precision": 80, "endpoint_decimal_precision": 100},
        "pending_before_execution": [
            "independent full-box proof/source-contract and family-exposure review",
            "resolve or retain every scientific exclusion without outcome-driven domain rescue",
            "complete source-equivalent stable log-truth adapters and high-precision certificates",
            "validate p5-8 diagnostics/full-quadratic rank and fixed external grids",
            "freeze dimension-matched null controls, exact task roles, grid/log-truth hashes and new seed namespaces",
        ],
    }
    _write_json(out / "review_manifest.json", manifest)
    return summary


def verify_registration(out):
    """Verify the metadata snapshot; never start or prepare study samples."""
    out = Path(out)
    manifest = json.loads((out / "review_manifest.json").read_text(encoding="utf-8"))
    for path, expected in manifest["source_and_input_sha256"].items():
        if digest(path) != expected:
            raise ValueError(f"Registration input drift:{path}")
    for name, expected in manifest["metadata_artifact_sha256"].items():
        if digest(out / name) != expected:
            raise ValueError(f"Metadata artifact drift:{name}")
    rows = json.loads((out / "candidates.json").read_text(encoding="utf-8"))["equations"]
    if len(rows) != 120 or len({row["eq_name"] for row in rows}) != 120:
        raise ValueError("Incomplete source census")
    if any(bool(row["inclusion"]) != (row["scientific_status"] == "eligible_analytic") for row in rows):
        raise ValueError("Scientific admission/status mismatch")
    if any(not row.get("inclusion_reason") and not row.get("exclusion_reason") for row in rows):
        raise ValueError("Missing per-source admission/exclusion reason")
    if any(row["inclusion"] and row["proof_obligations"] for row in rows):
        raise ValueError("Admitted formula has an unresolved scientific obligation")
    if any("positive_laurent_sum_tag" not in row for row in rows):
        raise ValueError("Fine grammar was not applied to the full source census")
    return {"verified": True, "equations": 120, "rules_sha256": manifest["rules_sha256"],
            "registration_status": manifest["registration_status"], "study_samples_generated": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="runs/p03_physical_protocol")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    result = verify_registration(ROOT / args.out) if args.verify else write_registration(ROOT / args.out)
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
