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
        assert detail["category_partition_pass"]
        assert detail["is_fully_compatible"]


def test_non_assimilable_categories_partition():
    breakdown = CohortCompatibilityAuditor.get_category_breakdown()
    assert breakdown["is_partition_exact"]
    assert breakdown["total_sample_size"] == 2758
    assert breakdown["cognitively_normal_cn"]["count"] == 1075
    assert breakdown["subjective_cognitive_decline_scd"]["count"] == 248
    assert breakdown["mild_cognitive_impairment_mci"]["count"] == 891
    assert breakdown["alzheimers_disease_dementia_ad"]["count"] == 544

    # Check non-assimilable gradient: ALPS(CN) > ALPS(SCD) > ALPS(MCI) > ALPS(AD)
    alps_cn = breakdown["cognitively_normal_cn"]["weighted_mean_alps"]
    alps_scd = breakdown["subjective_cognitive_decline_scd"]["weighted_mean_alps"]
    alps_mci = breakdown["mild_cognitive_impairment_mci"]["weighted_mean_alps"]
    alps_ad = breakdown["alzheimers_disease_dementia_ad"]["weighted_mean_alps"]

    assert alps_cn > alps_scd > alps_mci > alps_ad
    assert abs(alps_mci - 1.3319) < 0.05, "MCI must track the analytical tipping point"
