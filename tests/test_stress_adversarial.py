"""
Stress, Boundary Condition, and Adversarial Fuzzing Tests for NGTD.
Validates numerical stability against extreme inputs, biological singularities,
and pathological boundary conditions.
"""

import math
import random
import pytest
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
)


def test_adversarial_alps_boundaries():
    """Tests zero, negative, and extreme hyper-convective ALPS indices."""
    # Complete hydraulic collapse floor (ALPS <= 1.0)
    for bad_alps in [-5.0, 0.0, 0.5, 0.99, 1.0]:
        v = calculate_effective_velocity(bad_alps)
        pe = calculate_peclet_number(bad_alps, peak_sleep=True)
        assert v == 0.0, f"Expected 0.0 velocity for ALPS={bad_alps}, got {v}"
        assert pe == 0.0, f"Expected 0.0 Péclet for ALPS={bad_alps}, got {pe}"

    # Supraphysiological ALPS (extreme upper bound clipping)
    for high_alps in [2.5, 5.0, 50.0]:
        v = calculate_effective_velocity(high_alps)
        assert not math.isnan(v)
        assert not math.isinf(v)
        # Should be bounded by maximum allowed scaling factor (1.5x reference sleep)
        assert v <= DEFAULT_GLYMPHATIC.v_sleep * 1.5


def test_sleep_fraction_extremes():
    """Tests total sleep deprivation vs continuous coma sleep."""
    alps = 1.50
    # Total insomnia (0 hours sleep)
    v_insomnia = calculate_effective_velocity(alps, sleep_fraction=0.0)
    # 24h continuous sleep
    v_coma = calculate_effective_velocity(alps, sleep_fraction=1.0)

    assert v_insomnia > 0.0
    assert v_coma > v_insomnia
    assert math.isclose(v_coma / v_insomnia, DEFAULT_GLYMPHATIC.v_sleep / DEFAULT_GLYMPHATIC.v_awake, rel_tol=1e-5)


def test_microglial_strain_saturation_and_bounds():
    """Ensures microglial strain index is strictly bounded in [0.0, 1.0]."""
    # Over-specified APOE alleles (e.g. 5 alleles) should clamp at 2
    strain_over = calculate_microglial_strain(apoe_e4_alleles=5, trem2_r47h=True, nlrp3_activated=True)
    assert strain_over == 1.0

    # Negative inputs should clamp at 0
    strain_neg = calculate_microglial_strain(apoe_e4_alleles=-2, trem2_r47h=False, nlrp3_activated=False)
    assert strain_neg == 0.0


def test_solver_grid_resolution_invariance():
    """Checks that varying radial grid density does not induce numerical divergence."""
    for n_pts in [10, 25, 50, 100, 500]:
        res = simulate_radial_pvs_profile(alps_index=1.40, apoe_e4_alleles=1, num_grid_points=n_pts)
        assert len(res["r_points_um"]) == n_pts
        assert len(res["c_profile"]) == n_pts
        # Monotonicity check: inner boundary concentration <= outer boundary
        assert res["c_profile"][0] <= res["c_profile"][-1]
        assert not any(math.isnan(c) or math.isinf(c) for c in res["c_profile"])


def test_adversarial_monte_carlo_fuzzing():
    """Fuzzes model across 500 random physiological and pathological combinations."""
    random.seed(42)
    for _ in range(500):
        rand_alps = random.uniform(0.8, 2.2)
        rand_apoe = random.choice([0, 1, 2])
        rand_trem2 = random.choice([True, False])
        rand_nlrp3 = random.choice([True, False])
        rand_sleep = random.uniform(0.1, 0.6)

        v = calculate_effective_velocity(rand_alps, sleep_fraction=rand_sleep)
        pe = calculate_peclet_number(rand_alps)
        rates = calculate_net_clearance_rate(
            alps_index=rand_alps,
            apoe_e4_alleles=rand_apoe,
            trem2_r47h=rand_trem2,
            nlrp3_activated=rand_nlrp3,
        )
        sim = simulate_radial_pvs_profile(
            alps_index=rand_alps,
            apoe_e4_alleles=rand_apoe,
            trem2_r47h=rand_trem2,
            nlrp3_activated=rand_nlrp3,
            num_grid_points=30,
        )

        assert v >= 0.0
        assert pe >= 0.0
        assert 0.0 <= rates["microglial_strain"] <= 1.0
        assert rates["k_total_clearance"] > 0.0
        assert all(0.0 <= c <= 1.0 for c in sim["c_profile"])
        assert not math.isnan(rates["effective_growth_rate"])
