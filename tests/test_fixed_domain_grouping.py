"""Pre-sampling leakage controls independent of simulated performance."""
import pytest

from experiments.fixed_domain.grouping import assert_group_isolation, group_unit


def task(name, expression, split="development", klass="H1", group="g", **extra):
    return {"task": name, "expression": expression, "split": split, "klass": klass,
            "structural_group": group, "design_group": group, "source_equations": [], **extra}


def test_bootstrap_unit_uses_null_design_and_alternative_structure():
    assert group_unit(task("a", "x", klass="H0", group="null-x")) == "design:null-x"
    assert group_unit(task("b", "asinh(x)", group="shape")) == "structure:shape"


def test_manual_groups_and_syntactic_clones_cannot_cross_any_split():
    with pytest.raises(ValueError, match="group crosses"):
        assert_group_isolation([task("a", "asinh(x)", "confirmation"), task("b", "LambertW(x)", "reserve")])
    with pytest.raises(ValueError, match="twins cross"):
        assert_group_isolation([task("a", "asinh(x)", "confirmation", group="a"), task("b", "7*asinh(y)", "reserve", group="b")])


def test_monomial_prefactors_and_variable_renaming_do_not_hide_core_twins():
    with pytest.raises(ValueError, match="twins cross"):
        assert_group_isolation([task("a", "asinh(x)", "confirmation", group="a"), task("b", "z*asinh(y)", "reserve", group="b")])


def test_old_identity_alias_and_hidden_formula_clone_are_development_only():
    for record in [task("feynman-i.10.7", "asinh(x)", "confirmation"),
                   task("new-name", "asinh(x)", "confirmation", source_equations=["FeynmanICh10Eq7"]),
                   task("new-name", "x0/sqrt(1-1.1126500560536185e-17*x1**2)", "confirmation")]:
        with pytest.raises(ValueError, match="requires development"):
            assert_group_isolation([record])


def test_null_designs_isolate_roles_and_slope_clones_require_explicit_exception():
    records = [task("a", "x**0.3*y**1.7", "reference", "H0", "a"),
               task("b", "7*v**1.7*u**0.3", "test", "H0", "b")]
    with pytest.raises(ValueError, match="null slope clones"):
        assert_group_isolation(records)
    result = assert_group_isolation(records, h0_structure_exception="Predeclared shared null slope control.")
    assert result["passed"] and len(result["h0_exact_slope_cross_split_collisions"]) == 1
    records[1]["expression"] = "u**0.7*v**1.2"
    assert assert_group_isolation(records)["h0_exact_slope_cross_split_collisions"] == []


def test_stable_function_spelling_and_frozen_scalar_parameters():
    with pytest.raises(ValueError, match="twins cross"):
        assert_group_isolation([task("a", "log1p(alpha*x)", "confirmation", group="a", params={"alpha": 2}),
                                task("b", "log(1+2*y)", "reserve", group="b")])
    result = assert_group_isolation([task("a", "asinh(x)", "reserve", group="a"),
                                     task("b", "LambertW(x)", "confirmation", group="b")])
    assert result["passed"] and result["n_tasks"] == 2
