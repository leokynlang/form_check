"""Expose the source-hashed registration checks to default pytest collection."""
from experiments.physical_coverage_v1.test_registration import (
    metadata,
    test_complete_census_keeps_math_and_execution_obligations_separate,
    test_disputed_fine_combination_is_not_a_primary_new_mechanism,
    test_saved_registration_is_read_only_and_rejects_artifact_or_input_drift,
)

__all__ = [
    "metadata",
    "test_complete_census_keeps_math_and_execution_obligations_separate",
    "test_disputed_fine_combination_is_not_a_primary_new_mechanism",
    "test_saved_registration_is_read_only_and_rejects_artifact_or_input_drift",
]
