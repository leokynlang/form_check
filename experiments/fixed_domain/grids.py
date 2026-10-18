"""Sample-independent external grids and fail-closed region error metrics.

The geometry uses only declared domains. Truth values are checked after the
complete grid has been built; invalid truth never selects or removes a point.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict

import numpy as np
from scipy.stats import qmc

GRID_SEED = 2026101803
AXIAL_POINTS_PER_DIRECTION = 8
SOBOL_POWER = 8


def _domains(task):
    dimension = int(task["n_features"])
    outer = np.asarray(task["outer_bounds"], dtype=float)
    train = np.asarray(task["train_bounds"], dtype=float)
    scales = tuple(task["coordinate_scales"])
    if dimension < 1 or dimension > 4:
        raise ValueError("Grid supports one through four declared features")
    if outer.shape != (dimension, 2) or train.shape != (dimension, 2):
        raise ValueError("Domain shape does not match n_features")
    if len(scales) != dimension or any(s not in {"log", "linear", "integer"} for s in scales):
        raise ValueError("Unsupported or missing coordinate scale")
    if not np.isfinite(outer).all() or not np.isfinite(train).all():
        raise ValueError("Domain bounds must be finite")
    if (outer[:, 0] <= 0).any() or (outer[:, 0] >= outer[:, 1]).any():
        raise ValueError("Outer domain must have strictly positive, ordered bounds")
    if ((train[:, 0] < outer[:, 0]).any() or (train[:, 1] > outer[:, 1]).any()
            or (train[:, 0] >= train[:, 1]).any()):
        raise ValueError("Training domain must be a nondegenerate nested box")
    for j, scale in enumerate(scales):
        if scale == "integer" and np.ceil(train[j, 0]) > np.floor(train[j, 1]):
            raise ValueError("Integer training domain has no support")
    return outer, train, scales


def _from_unit(unit, outer, scales):
    unit = np.asarray(unit, dtype=float)
    result = np.empty_like(unit)
    for j, scale in enumerate(scales):
        lo, hi = outer[j]
        if scale == "log":
            result[:, j] = np.exp(np.log(lo) + unit[:, j] * (np.log(hi) - np.log(lo)))
        else:
            result[:, j] = lo + unit[:, j] * (hi - lo)
        if scale == "integer":
            result[:, j] = np.clip(np.rint(result[:, j]), np.ceil(lo), np.floor(hi))
        else:
            # Preserve the exact declared endpoints rather than exp(log(bound)).
            result[unit[:, j] == 0, j] = lo
            result[unit[:, j] == 1, j] = hi
    return result


def _to_unit(values, outer, scales):
    result = np.empty_like(values, dtype=float)
    for j, scale in enumerate(scales):
        lo, hi = outer[j]
        if scale == "log":
            result[:, j] = (np.log(values[:, j]) - np.log(lo)) / (np.log(hi) - np.log(lo))
        else:
            result[:, j] = (values[:, j] - lo) / (hi - lo)
    return result


def _external(points, train):
    return np.any((points < train[:, 0]) | (points > train[:, 1]), axis=1)


def _unique_rows(points):
    # Lexicographic ordering makes support and hashes independent of set order.
    return np.unique(np.asarray(points, dtype=float), axis=0)


def _array_hash(array):
    array = np.ascontiguousarray(array, dtype="<f8")
    shape = json.dumps(list(array.shape), separators=(",", ":")).encode("ascii")
    return hashlib.sha256(shape + b"\n<float64-le>\n" + array.tobytes()).hexdigest()


def build_grid(task, seed=GRID_SEED):
    """Build an external grid using fixed domains, never training observations.

    Half of the probability mass is divided equally across viable lower/upper
    axial directions, then equally across their unique points. The other half
    is uniform on the union of external scrambled-Sobol points and outer-box
    corners. Integer rounding and deduplication happen within each component
    before assigning its mass. Support shared across components is merged by
    adding weights; every component label is retained on the merged point.
    """
    outer, train, scales = _domains(task)
    dimension = len(scales)
    anchor = np.empty(dimension)
    for j, scale in enumerate(scales):
        if scale == "log":
            anchor[j] = np.exp((np.log(train[j, 0]) + np.log(train[j, 1])) / 2)
        else:
            anchor[j] = np.mean(train[j])
        if scale == "integer":
            anchor[j] = np.clip(np.rint(anchor[j]), np.ceil(train[j, 0]), np.floor(train[j, 1]))
    train_unit = _to_unit(train.T, outer, scales)
    anchor_unit = _to_unit(anchor[None, :], outer, scales)[0]
    components = {}
    unavailable = []
    for j in range(dimension):
        for side in ("lower", "upper"):
            label = f"axis_x{j + 1}_{side}"
            if side == "lower":
                values = np.linspace(0, train_unit[0, j], AXIAL_POINTS_PER_DIRECTION + 1)[:-1]
            else:
                values = np.linspace(train_unit[1, j], 1, AXIAL_POINTS_PER_DIRECTION + 1)[1:]
            unit = np.tile(anchor_unit, (AXIAL_POINTS_PER_DIRECTION, 1))
            unit[:, j] = values
            points = _from_unit(unit, outer, scales)
            # Fix all other coordinates to the exact registered-scale center.
            other = np.arange(dimension) != j
            points[:, other] = anchor[other]
            selected = points[:, j] < train[j, 0] if side == "lower" else points[:, j] > train[j, 1]
            points = _unique_rows(points[selected])
            if len(points):
                components[label] = points
            else:
                unavailable.append(label)
    if not components:
        raise ValueError("No representable axial points outside the training box")

    sobol_unit = qmc.Sobol(d=dimension, scramble=True, seed=int(seed)).random_base2(SOBOL_POWER)
    sobol = _unique_rows(_from_unit(sobol_unit, outer, scales))
    sobol = sobol[_external(sobol, train)]
    corner_unit = np.array([[float((i >> j) & 1) for j in range(dimension)]
                            for i in range(2 ** dimension)])
    corners = _unique_rows(_from_unit(corner_unit, outer, scales))
    corners = corners[_external(corners, train)]
    joint = _unique_rows(np.vstack([sobol, corners]))
    if not len(joint):
        raise ValueError("No joint points outside the training box")

    point_weights = defaultdict(float)
    point_labels = defaultdict(set)
    component_metadata = {}
    mass = .5 / len(components)
    for label, points in components.items():
        component_metadata[label] = {"mass": mass, "unique_points": len(points)}
        for point in points:
            key = tuple(point)
            point_weights[key] += mass / len(points)
            point_labels[key].add(label)
    for point in joint:
        key = tuple(point)
        point_weights[key] += .5 / len(joint)
        point_labels[key].add("joint")
    for label, points in (("joint_sobol", sobol), ("joint_corner", corners)):
        for point in points:
            point_labels[tuple(point)].add(label)
    keys = sorted(point_weights)
    grid = {
        "X": np.asarray(keys, dtype=float),
        "weights": np.asarray([point_weights[key] for key in keys], dtype=float),
        "directions": [sorted(point_labels[key]) for key in keys],
        "metadata": {
            "grid_version": "fixed-external-axial-sobol-v1",
            "seed": int(seed),
            "axial_points_requested_per_direction": AXIAL_POINTS_PER_DIRECTION,
            "axial_total_mass": .5,
            "axial_components": component_metadata,
            "unavailable_axial_directions": unavailable,
            "joint_total_mass": .5,
            "joint_unique_points": len(joint),
            "external_sobol_unique_points": len(sobol),
            "external_corner_unique_points": len(corners),
            "sobol_candidates": 2 ** SOBOL_POWER,
            "sobol_scramble": True,
            "joint_weight_rule": "uniform_over_unique_external_Sobol_and_corner_union",
            "integer_rule": "nearest_integer_with_numpy_rint_then_clip_to_outer_integer_support",
            "anchor": anchor.tolist(),
            "geometry_depends_on": "declared_domains_scales_dimension_and_fixed_seed_only",
            "truth_point_selection": "none; invalid truth fails validation without point removal",
        },
    }
    return grid


def _validate_weights(weights, length):
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (length,) or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError("Every grid weight must be finite and strictly positive")
    if not np.isclose(weights.sum(), 1.0, rtol=0, atol=1e-12):
        raise ValueError("Grid weights must sum to one")
    return weights


def validate_grid(task, grid, truth):
    """Validate the complete prebuilt grid; return hashes and validity facts.

    ``truth`` may be a precomputed vector or an X -> truth-vector callable.
    No failed truth point is removed or replaced. A nonempty analytic domain
    argument is required, though its mathematical soundness needs human review.
    """
    outer, train, scales = _domains(task)
    points = np.asarray(grid["X"], dtype=float)
    if points.ndim != 2 or points.shape[1] != len(scales) or not len(points):
        raise ValueError("Grid coordinates have invalid dimensions or are empty")
    if not np.isfinite(points).all() or (points <= 0).any():
        raise ValueError("Grid inputs must all be finite and strictly positive")
    if ((points < outer[:, 0]).any() or (points > outer[:, 1]).any()):
        raise ValueError("Grid lies outside the declared valid outer domain")
    if not _external(points, train).all():
        raise ValueError("Grid contains a point inside or on the training box")
    if len(_unique_rows(points)) != len(points):
        raise ValueError("Grid support contains duplicate points")
    for j, scale in enumerate(scales):
        if scale == "integer" and not np.equal(points[:, j], np.rint(points[:, j])).all():
            raise ValueError("Integer grid coordinate is not integral")
    weights = _validate_weights(grid["weights"], len(points))
    if len(grid["directions"]) != len(points) or any(not labels for labels in grid["directions"]):
        raise ValueError("Each grid point requires direction labels")
    argument = task.get("validity_argument")
    if not isinstance(argument, str) or not argument.strip():
        raise ValueError("Missing analytic validity argument for the full outer box")
    truth = np.asarray(truth(points) if callable(truth) else truth, dtype=float)
    if truth.shape != (len(points),) or not np.isfinite(truth).all() or (truth <= 0).any():
        raise ValueError("Grid truth must be finite and strictly positive at every point")
    return {
        "validated": True,
        "n_points": len(points),
        "n_features": len(scales),
        "finite_positive_inputs": True,
        "finite_positive_truth": True,
        "inside_outer_domain": True,
        "outside_train": True,
        "analyticproofpresent": True,
        "validity_argument_requires_separate_review": True,
        "truth_min": float(truth.min()),
        "truth_max": float(truth.max()),
        "weight_sum": float(weights.sum()),
        "coordinates_sha256": _array_hash(points),
        "weights_sha256": _array_hash(weights),
        "truth_sha256": _array_hash(truth),
    }


def error_metrics(pred_log, truth, weights):
    """Evaluate every fixed point, or fail the entire outcome explicitly.

    Positive compact-domain truth permits pointwise relative error without an
    arbitrary denominator floor. expm1 computes abs(pred/truth - 1) stably. An
    overflow or other nonfinite result is an unavailable outcome, never a point
    to discard. The weighted p95 is the first observed value with CDF >= .95.
    """
    truth = np.asarray(truth, dtype=float)
    pred_log = np.asarray(pred_log, dtype=float)
    if truth.ndim != 1 or not len(truth) or pred_log.shape != truth.shape:
        raise ValueError("Prediction and truth must be matching nonempty vectors")
    if not np.isfinite(truth).all() or (truth <= 0).any():
        raise ValueError("Every truth value must be finite and strictly positive")
    if not np.isfinite(pred_log).all():
        raise ValueError("Every log prediction must be finite")
    weights = _validate_weights(weights, len(truth))
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        errors = np.abs(np.expm1(pred_log - np.log(truth)))
    if not np.isfinite(errors).all():
        raise ValueError("nonfinite_relative_error: whole grid outcome unavailable")
    order = np.argsort(errors, kind="stable")
    cdf = np.cumsum(weights[order])
    p95_index = min(int(np.searchsorted(cdf, .95, side="left")), len(order) - 1)
    mean = float(weights @ errors)
    if not np.isfinite(mean):
        raise ValueError("nonfinite_region_mean: whole grid outcome unavailable")
    return {
        "weighted_mean_relative_error": mean,
        "weighted_p95_relative_error": float(errors[order[p95_index]]),
        "max_relative_error": float(errors.max()),
    }
