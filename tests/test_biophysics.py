"""
Unit tests for Neuro-Glymphatic Tau Dynamics (NGTD) engine.
"""

import pytest
import math
from tau_glymphatic import (
    DEFAULT_GEOMETRY,
    DEFAULT_PHYSICS,
    DEFAULT_GLYMPHATIC,
    calculate_effective_velocity,
    calculate_peclet_number,
    solve_alps_critical_threshold,
    calculate_microglial_strain,
    calculate_net_clearance_rate,
    simulate_radial_pvs_profile,
    run_full_biophysical_audit,
)


def test_geometric_and_physical_constants():
    assert DEFAULT_GEOMETRY.r_vessel == 10.0e-6
    assert DEFAULT_GEOMETRY.r_pvs_outer == 30.0e-6
    assert math.isclose(DEFAULT_GEOMETRY.pvs_width, 20.0e-6)
    assert DEFAULT_PHYSICS.porosity == 0.20
    assert DEFAULT_PHYSICS.tortuosity == 1.60
    assert DEFAULT_PHYSICS.d_eff_tau < DEFAULT_PHYSICS.d0_tau


def test_alps_velocity_monotonicity():
    v1 = calculate_effective_velocity(1.00)
    v2 = calculate_effective_velocity(1.20)
    v3 = calculate_effective_velocity(1.60)
    assert v1 == 0.0
    assert 0.0 < v2 < v3


def test_peclet_number_phase_transition():
    # In healthy sleep, Pe should exceed 1.0 (convection dominates)
    pe_healthy = calculate_peclet_number(DEFAULT_GLYMPHATIC.alps_healthy_ref, peak_sleep=True)
    assert pe_healthy > 1.0

    # In AD, Pe should drop significantly below 1.0 (diffusion/stasis dominates)
    pe_ad = calculate_peclet_number(DEFAULT_GLYMPHATIC.alps_ad_typical, peak_sleep=True)
    assert pe_ad < 0.5


def test_microglial_strain_genetics():
    s0 = calculate_microglial_strain(apoe_e4_alleles=0, trem2_r47h=False, nlrp3_activated=False)
    s1 = calculate_microglial_strain(apoe_e4_alleles=1, trem2_r47h=False, nlrp3_activated=False)
    s2 = calculate_microglial_strain(apoe_e4_alleles=2, trem2_r47h=False, nlrp3_activated=False)
    s3 = calculate_microglial_strain(apoe_e4_alleles=2, trem2_r47h=True, nlrp3_activated=True)

    assert s0 == 0.0
    assert s1 == 0.25
    assert s2 == 0.50
    assert s3 == 1.00  # Capped at 1.0


def test_net_clearance_bifurcation():
    # Healthy patient: negative net growth rate (active washout)
    healthy = calculate_net_clearance_rate(alps_index=1.60, apoe_e4_alleles=0)
    assert not healthy["is_pathologic_expansion"]
    assert healthy["effective_growth_rate"] < 0.0
    assert healthy["net_stasis_ratio"] < 1.0

    # Severe AD with genetic risk: positive net growth rate (pathologic seeding)
    ad = calculate_net_clearance_rate(alps_index=1.04, apoe_e4_alleles=2, nlrp3_activated=True)
    assert ad["is_pathologic_expansion"]
    assert ad["effective_growth_rate"] > 0.0
    assert ad["net_stasis_ratio"] > 1.0


def test_radial_pvs_profile_simulation():
    res = simulate_radial_pvs_profile(alps_index=1.50, num_grid_points=25)
    assert len(res["r_points_um"]) == 25
    assert len(res["c_profile"]) == 25
    assert all(0.0 <= c <= 1.0 for c in res["c_profile"])
    # Outer parenchyma concentration should be higher than perivascular drainage interface
    assert res["c_profile"][0] <= res["c_profile"][-1]


def test_full_biophysical_audit_passes():
    audit = run_full_biophysical_audit()
    assert audit["status"] == "PASSED"
    assert audit["passed_checks"] == 7
    assert audit["total_checks"] == 7
