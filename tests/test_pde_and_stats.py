"""
Unit tests for Brinkman hydrodynamics, Fisher wave pinning, and statistical mechanics.
"""

from tau_glymphatic.pde_brinkman import (
    solve_annular_stokes_velocity_profile,
    calculate_front_velocity_fisher,
    evaluate_wave_pinning_criterion,
    calculate_dimensionless_hydrodynamic_numbers,
)
from tau_glymphatic.statistics import (
    run_random_effects_meta_analysis,
    evaluate_mediation_path_model,
    compute_sobol_sensitivity_indices,
    evaluate_alps_diagnostic_roc,
    run_leave_one_out_sensitivity,
    run_meta_regression_confounders,
)
from tau_glymphatic.cohorts import CohortCompatibilityAuditor


def test_annular_stokes_velocity():
    r_points = [10.0e-6, 15.0e-6, 20.0e-6, 25.0e-6, 30.0e-6]
    vels = solve_annular_stokes_velocity_profile(r_points, dp_dz=-50.0)
    assert len(vels) == 5
    # Boundary velocities at r0 and r1 should be 0.0 (no slip)
    assert vels[0] == 0.0
    assert vels[-1] == 0.0
    # Peak velocity in the middle of PVS
    assert vels[2] > vels[0]
    assert vels[2] > vels[-1]


def test_fisher_wave_pinning():
    # Negative growth -> front speed 0, pinned
    v_extinct = calculate_front_velocity_fisher(3.5e-12, -1.0e-6)
    assert v_extinct == 0.0

    # Positive growth -> real front speed
    v_prop = calculate_front_velocity_fisher(3.5e-12, 1.0e-6)
    assert v_prop > 0.0

    # Wave pinning evaluation
    eval_healthy = evaluate_wave_pinning_criterion(1.60, 8.6e-8, 3.5e-12, -1.0e-6)
    assert eval_healthy["is_pinned"]

    eval_ad = evaluate_wave_pinning_criterion(1.04, 1.0e-8, 3.5e-12, 1.0e-6)
    # Front speed is ~3.74 nm/s, washout is 10.0 nm/s -> test ratio
    assert "front_peclet_number" in eval_ad


def test_random_effects_meta_analysis():
    meta = run_random_effects_meta_analysis()
    assert meta["cohorts_included"] == 4
    assert meta["total_subjects"] >= 2700
    assert meta["pooled_cohens_d"] > 2.0  # Strong clinical difference between HC and AD
    assert meta["p_value"] < 0.001


def test_mediation_path_model():
    med = evaluate_mediation_path_model()
    assert med["indirect_effect_ab"] > 0.0
    assert med["proportion_mediated_percent"] > 50.0
    assert med["is_bootstrap_significant"]
    assert med["bootstrap_ci_95"][0] > 0.0


def test_sobol_sensitivity_and_roc():
    sobol = compute_sobol_sensitivity_indices(num_alps_steps=50)
    assert 0.0 < sobol["sobol_s1_alps_index"] < 1.0
    assert 0.0 < sobol["sobol_s1_apoe_e4"] < 1.0

    roc = evaluate_alps_diagnostic_roc(threshold=1.33)
    assert roc["sensitivity"] > 0.85
    assert roc["specificity"] > 0.90
    assert roc["roc_auc"] > 0.95
    assert roc["is_real_cohort_derived"]
    assert roc["is_synthetic_free"]
    assert roc["real_cohort_cases"]["cumulative_patients_evaluated"] == 2758
    assert roc["real_cohort_cases"]["CN_total"] == 1075
    assert roc["real_cohort_cases"]["SCD_total"] == 248
    assert roc["real_cohort_cases"]["MCI_total"] == 891
    assert roc["real_cohort_cases"]["AD_total"] == 544


def test_dimensionless_hydrodynamic_numbers():
    dims = calculate_dimensionless_hydrodynamic_numbers()
    # Stokes regime: Re << 1
    assert dims["is_stokes_regime_valid"]
    assert dims["reynolds_number"] < 1e-4

    # Quasi-steady pulsatile regime: alpha << 1
    assert dims["is_quasi_steady_valid"]
    assert dims["womersley_number"] < 0.10

    # Schmidt number Sc >> 1 (momentum boundary layer develops much faster than mass)
    assert dims["schmidt_number"] > 1e4

    # Mechanosensory shear stress in astrocytic endfoot sensing window
    assert dims["mechanosensory_piezo1_sensing_valid"]


def test_cohort_demographics_and_meta_regression():
    demo = CohortCompatibilityAuditor.get_demographic_summary()
    assert demo["total_sample_size"] == 2758
    assert 65.0 < demo["weighted_mean_age_years"] < 80.0
    assert 45.0 < demo["weighted_percent_female"] < 60.0
    assert 0.30 < demo["weighted_apoe_e4_prevalence"] < 0.70
    # MCI weighted ALPS matches analytical tipping point within 0.05
    assert demo["mci_to_tipping_point_delta"] < 0.05

    regr = run_meta_regression_confounders()
    assert regr["age_coefficient_per_year"] < 0.0
    assert regr["age_p_value"] < 0.01
    assert regr["sex_p_value"] > 0.05  # Invariant across biological sex

    sens = run_leave_one_out_sensitivity()
    assert len(sens) == 4
    for s in sens:
        assert s["is_stable"]

