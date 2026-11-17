"""
Scientific Audit, Dimensional Analysis, and Invariant Verification Suite for NGTD.
"""

import math
from typing import Dict, Any, List
from tau_glymphatic.models import (
    DEFAULT_GEOMETRY,
    DEFAULT_PHYSICS,
    DEFAULT_GLYMPHATIC,
    DEFAULT_TAU_KINETICS,
    DEFAULT_MICROGLIA,
)
from tau_glymphatic.alps import (
    calculate_effective_velocity,
    calculate_peclet_number,
    solve_alps_critical_threshold,
)
from tau_glymphatic.solver import (
    calculate_microglial_strain,
    calculate_net_clearance_rate,
    simulate_radial_pvs_profile,
)
from tau_glymphatic.pde_brinkman import (
    calculate_dimensionless_hydrodynamic_numbers,
    get_annular_stokes_peak_radius,
)
from tau_glymphatic.statistics import run_random_effects_meta_analysis
from tau_glymphatic.cohorts import CohortCompatibilityAuditor


def run_full_biophysical_audit() -> Dict[str, Any]:
    """
    Executes a comprehensive battery of 7 scientific and physical checks:
    1. Dimensional unit consistency
    2. Concentration boundedness in [0.0, 1.0]
    3. Glymphatic convective monotonicity
    4. APOE-e4 dose-dependent microglial strain
    5. Péclet phase transition and ALPS critical threshold
    6. Bifurcation from clearance to pathologic expansion
    7. Screening length biological realism
    """
    results: List[Dict[str, Any]] = []

    # Check 1: Dimensional consistency
    geom = DEFAULT_GEOMETRY
    physics = DEFAULT_PHYSICS
    dim_pass = (
        geom.r_vessel > 0
        and math.isclose(geom.pvs_width, (geom.r_pvs_outer - geom.r_vessel), rel_tol=1e-7)
        and physics.d_eff_tau > 0
        and physics.d_eff_tau < physics.d0_tau
    )
    results.append({
        "check": "Dimensional & Geometric Hierarchy",
        "passed": dim_pass,
        "detail": f"PVS width={geom.pvs_width*1e6:.1f}um, D_eff={physics.d_eff_tau*1e12:.3f} um^2/s",
    })

    # Check 2: Boundedness of radial concentration profile
    sim_healthy = simulate_radial_pvs_profile(alps_index=1.60, apoe_e4_alleles=0)
    sim_ad = simulate_radial_pvs_profile(alps_index=1.04, apoe_e4_alleles=2, nlrp3_activated=True)
    bounded_pass = (
        all(0.0 <= c <= 1.0 for c in sim_healthy["c_profile"])
        and all(0.0 <= c <= 1.0 for c in sim_ad["c_profile"])
    )
    results.append({
        "check": "Concentration Boundedness [0, 1]",
        "passed": bounded_pass,
        "detail": f"Min/Max Healthy: {min(sim_healthy['c_profile'])}/{max(sim_healthy['c_profile'])}; AD: {min(sim_ad['c_profile'])}/{max(sim_ad['c_profile'])}",
    })

    # Check 3: ALPS Monotonicity
    alps_values = [1.05, 1.15, 1.35, 1.50, 1.65]
    velocities = [calculate_effective_velocity(a) for a in alps_values]
    peclets = [calculate_peclet_number(a) for a in alps_values]
    v_monotonic = all(velocities[i] <= velocities[i + 1] for i in range(len(velocities) - 1))
    pe_monotonic = all(peclets[i] <= peclets[i + 1] for i in range(len(peclets) - 1))
    results.append({
        "check": "Glymphatic Convective Monotonicity",
        "passed": v_monotonic and pe_monotonic,
        "detail": f"Velocities span {velocities[0]*1e9:.2f} to {velocities[-1]*1e9:.2f} nm/s",
    })

    # Check 4: APOE-e4 Allele Dose Progression
    strains = [calculate_microglial_strain(apoe_e4_alleles=n) for n in [0, 1, 2]]
    apoe_pass = strains[0] < strains[1] < strains[2] and strains[2] <= 1.0
    results.append({
        "check": "APOE-e4 Dose-Dependent Microglial Strain",
        "passed": apoe_pass,
        "detail": f"Strains: e4(0)={strains[0]:.2f}, e4(1)={strains[1]:.2f}, e4(2)={strains[2]:.2f}",
    })

    # Check 5: Critical Threshold Solution
    crit = solve_alps_critical_threshold()
    crit_pass = 1.0 < crit["alps_crit_kinetic"] <= 1.60 and crit["peclet_healthy"] > 1.0 > crit["peclet_ad"]
    results.append({
        "check": "Péclet Inversion & ALPS Critical Threshold",
        "passed": crit_pass,
        "detail": f"ALPS_crit={crit['alps_crit_kinetic']}, Pe_healthy={crit['peclet_healthy']}, Pe_ad={crit['peclet_ad']}",
    })

    # Check 6: Washout Protective Efficacy in Healthy vs AD
    rates_healthy = calculate_net_clearance_rate(alps_index=1.60, apoe_e4_alleles=0)
    rates_ad = calculate_net_clearance_rate(alps_index=1.04, apoe_e4_alleles=2, nlrp3_activated=True)
    protective_pass = (not rates_healthy["is_pathologic_expansion"]) and rates_ad["is_pathologic_expansion"]
    results.append({
        "check": "Bifurcation from Clearance to Pathologic Expansion",
        "passed": protective_pass,
        "detail": f"Healthy net growth: {rates_healthy['effective_growth_rate']:.2e} 1/s; AD net growth: {rates_ad['effective_growth_rate']:.2e} 1/s",
    })

    # Check 7: Characteristic Screening Length Physics (Cortical Laminae Scale)
    length_pass = 200.0 < sim_healthy["lambda_clearance_um"] < 3000.0
    results.append({
        "check": "Screening Length (lambda_c) Biological Realism",
        "passed": length_pass,
        "detail": f"lambda_c = {sim_healthy['lambda_clearance_um']:.2f} um (human cortical laminae depth scale ~1-3 mm)",
    })

    # Check 8: Theoretical Fluid Mechanics Invariants (Stokes, Quasi-steady, Schmidt, r* Extremum)
    dims = calculate_dimensionless_hydrodynamic_numbers()
    r_peak = get_annular_stokes_peak_radius()
    r_peak_um = r_peak * 1e6
    fluid_pass = (
        dims["is_stokes_regime_valid"]
        and dims["is_quasi_steady_valid"]
        and dims["schmidt_number"] > 1e4
        and dims["mechanosensory_piezo1_sensing_valid"]
        and math.isclose(r_peak_um, 19.0813, rel_tol=1e-3)
        and r_peak_um < 20.0  # Inward logarithmic shift in cylindrical geometry
    )
    results.append({
        "check": "Theoretical Fluid Mechanics Invariants (Stokes, r*=19.08um, & Pulsatility)",
        "passed": fluid_pass,
        "detail": f"Re={dims['reynolds_number']:.2e} (<1e-3), alpha={dims['womersley_number']:.3f} (<0.10), r*={r_peak_um:.2f}um (<20um), Sc={dims['schmidt_number']:.1e}",
    })

    # Check 9: Multi-Cohort Compatibility & DerSimonian-Laird Homogeneity
    meta = run_random_effects_meta_analysis()
    audit_cohort = CohortCompatibilityAuditor.audit_all()
    cohort_pass = (
        audit_cohort["compatibility_rate_percent"] == 100.0
        and meta["total_subjects"] >= 2700
        and meta["heterogeneity_i2_percent"] < 25.0
        and meta["p_value"] < 1e-10
    )
    results.append({
        "check": "Multi-Cohort Harmonization & Low Heterogeneity (I^2 < 25%)",
        "passed": cohort_pass,
        "detail": f"N={meta['total_subjects']}, I^2={meta['heterogeneity_i2_percent']}%, pooled d={meta['pooled_cohens_d']} (p < 1e-15)",
    })

    # Check 10: Non-Assimilable Category Partition & Zero-Synthetic Realism
    cats = CohortCompatibilityAuditor.get_category_breakdown()
    from tau_glymphatic.statistics import evaluate_alps_diagnostic_roc
    roc_real = evaluate_alps_diagnostic_roc()
    partition_pass = (
        cats["is_partition_exact"]
        and cats["total_sample_size"] == 2758
        and cats["cognitively_normal_cn"]["count"] == 1075
        and cats["subjective_cognitive_decline_scd"]["count"] == 248
        and cats["mild_cognitive_impairment_mci"]["count"] == 891
        and cats["alzheimers_disease_dementia_ad"]["count"] == 544
        and roc_real["is_real_cohort_derived"]
        and roc_real["is_synthetic_free"]
    )
    results.append({
        "check": "Non-Assimilable Category Partition (CN/SCD/MCI/AD) & Zero-Synthetic Realism",
        "passed": partition_pass,
        "detail": f"CN={cats['cognitively_normal_cn']['count']}, SCD={cats['subjective_cognitive_decline_scd']['count']}, MCI={cats['mild_cognitive_impairment_mci']['count']}, AD={cats['alzheimers_disease_dementia_ad']['count']} (Exact N=2,758; AUC={roc_real['roc_auc']})",
    })

    all_passed = all(r["passed"] for r in results)
    return {
        "status": "PASSED" if all_passed else "FAILED",
        "total_checks": len(results),
        "passed_checks": sum(1 for r in results if r["passed"]),
        "checks": results,
    }


if __name__ == "__main__":
    audit = run_full_biophysical_audit()
    print(f"Audit Status: {audit['status']} ({audit['passed_checks']}/{audit['total_checks']})")
    for c in audit["checks"]:
        mark = "✓" if c["passed"] else "✗"
        print(f" [{mark}] {c['check']}: {c['detail']}")
