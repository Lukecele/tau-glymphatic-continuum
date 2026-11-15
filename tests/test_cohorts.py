"""
Unit tests for multi-cohort clinical registry and compatibility gate.
"""

from tau_glymphatic.cohorts import AUDITED_COHORTS, CohortCompatibilityAuditor


def test_audited_cohorts_existence_and_size():
    assert len(AUDITED_COHORTS) >= 4
    total_n = sum(c.sample_size for c in AUDITED_COHORTS)
    assert total_n >= 2700, f"Expected N >= 2,700, got {total_n}"


def test_cohort_compatibility_audit_gate():
    summary = CohortCompatibilityAuditor.audit_all()
    assert summary["fully_compatible_cohorts"] == summary["total_cohorts"]
    assert summary["compatibility_rate_percent"] == 100.0
    assert summary["compatible_patients"] == summary["total_patients"]

    for detail in summary["cohort_details"]:
        assert detail["field_strength_3t_pass"]
        assert detail["high_angular_directions_pass"]
        assert detail["taoka_roi_standard_pass"]
        assert detail["biomarker_confirmed_pass"]
        assert detail["monotonic_alps_drop_pass"]
        assert detail["is_fully_compatible"]
