"""
Mathematical Solvers for Radial Tau Concentration Profiles and Spreading Dynamics.

Implements analytical steady-state and 1D finite-difference transient solvers
for tau concentration in perivascular and parenchymal domains, coupled to
glymphatic clearance and microglial inflammatory strain.
"""

import math
from typing import Dict, Any, List
from tau_glymphatic.models import (
    DEFAULT_GEOMETRY,
    DEFAULT_PHYSICS,
    DEFAULT_GLYMPHATIC,
    DEFAULT_TAU_KINETICS,
    DEFAULT_MICROGLIA,
    PerivascularGeometry,
    BrainTissuePhysics,
    GlymphaticFlowParameters,
    TauKinetics,
    MicroglialNeuroinflammation,
)
from tau_glymphatic.alps import calculate_effective_velocity


def calculate_microglial_strain(
    apoe_e4_alleles: int = 0,
    trem2_r47h: bool = False,
    nlrp3_activated: bool = False,
    microglia_params: MicroglialNeuroinflammation = DEFAULT_MICROGLIA,
) -> float:
    """
    Computes composite microglial neuroinflammatory strain index I_strain in [0.0, 1.0].
    """
    strain = 0.0
    strain += min(2, max(0, apoe_e4_alleles)) * microglia_params.apoe_e4_per_allele
    if trem2_r47h:
        strain += microglia_params.trem2_r47h_penalty
    if nlrp3_activated:
        strain += microglia_params.nlrp3_hyperactivation
    return min(1.0, max(0.0, strain))


def calculate_net_clearance_rate(
    alps_index: float,
    apoe_e4_alleles: int = 0,
    trem2_r47h: bool = False,
    nlrp3_activated: bool = False,
    sleep_fraction: float = 0.3333,
    geom: PerivascularGeometry = DEFAULT_GEOMETRY,
    glymph_params: GlymphaticFlowParameters = DEFAULT_GLYMPHATIC,
    tau_kinetics: TauKinetics = DEFAULT_TAU_KINETICS,
    microglia_params: MicroglialNeuroinflammation = DEFAULT_MICROGLIA,
) -> Dict[str, float]:
    """
    Calculates net local biological clearance rate k_net (1/s) combining:
    1. Basal enzymatic/autophagic degradation
    2. Perivascular glymphatic washout: k_glymph = theta * (v_eff / L_pvs)
    3. Microglial phagocytosis minus inflammatory exosomal seeding amplification
    """
    v_eff = calculate_effective_velocity(alps_index, glymph_params, sleep_fraction=sleep_fraction)
    theta = glymph_params.tau_solvent_drag_fraction
    k_glymph = theta * (v_eff / geom.pvs_width)

    strain = calculate_microglial_strain(
        apoe_e4_alleles=apoe_e4_alleles,
        trem2_r47h=trem2_r47h,
        nlrp3_activated=nlrp3_activated,
        microglia_params=microglia_params,
    )

    phagocytic_deficit = microglia_params.k_phagocytosis_base * strain
    if strain < microglia_params.strain_threshold_flip:
        inflammatory_seeding = 0.0
    else:
        excess_strain = strain - microglia_params.strain_threshold_flip
        inflammatory_seeding = microglia_params.seeding_boost_rate * (
            excess_strain / (1.0 - microglia_params.strain_threshold_flip)
        )

    # Autocatalytic seeding drive augmented by neuroinflammation and microglial impairment
    total_aggregation_drive = tau_kinetics.rho_aggregation + inflammatory_seeding + phagocytic_deficit

    # The effective growth rate: if k_glymph exceeds aggregation drive, net growth is negative (washout)
    effective_growth_rate = total_aggregation_drive - k_glymph
    k_total_clearance = k_glymph + microglia_params.k_phagocytosis_base * max(0.05, 1.0 - strain)
    net_stasis_ratio = total_aggregation_drive / max(1e-12, k_glymph)

    return {
        "k_glymph": k_glymph,
        "k_total_clearance": k_total_clearance,
        "microglial_strain": strain,
        "inflammatory_seeding_boost": inflammatory_seeding,
        "effective_growth_rate": effective_growth_rate,
        "net_stasis_ratio": net_stasis_ratio,
        "is_pathologic_expansion": effective_growth_rate > 0.0,
    }


def simulate_radial_pvs_profile(
    alps_index: float,
    apoe_e4_alleles: int = 0,
    trem2_r47h: bool = False,
    nlrp3_activated: bool = False,
    sleep_fraction: float = 0.3333,
    num_grid_points: int = 50,
    geom: PerivascularGeometry = DEFAULT_GEOMETRY,
    physics: BrainTissuePhysics = DEFAULT_PHYSICS,
    glymph_params: GlymphaticFlowParameters = DEFAULT_GLYMPHATIC,
) -> Dict[str, Any]:
    """
    Computes steady-state radial tau concentration profile C(r) across the
    perivascular space and into the surrounding neuropil: r in [r_vessel, r_parenchyma].
    """
    rates = calculate_net_clearance_rate(
        alps_index=alps_index,
        apoe_e4_alleles=apoe_e4_alleles,
        trem2_r47h=trem2_r47h,
        nlrp3_activated=nlrp3_activated,
        sleep_fraction=sleep_fraction,
        geom=geom,
        glymph_params=glymph_params,
    )

    r_start = geom.r_vessel
    r_end = geom.r_parenchyma
    dr = (r_end - r_start) / (num_grid_points - 1)

    r_points_um: List[float] = []
    c_profile: List[float] = []

    k_clear = max(1e-9, rates["k_total_clearance"])
    lambda_c = math.sqrt(physics.d_eff_tau / k_clear)
    lambda_c_um = lambda_c * 1e6

    # Reference healthy glymphatic clearance rate
    v_ref = calculate_effective_velocity(glymph_params.alps_healthy_ref, glymph_params, sleep_fraction=0.3333)
    k_ref = glymph_params.tau_solvent_drag_fraction * (v_ref / geom.pvs_width)
    drainage_ratio = rates["k_glymph"] / max(1e-12, k_ref)

    # Boundary concentration at vessel wall r_vessel
    strain = rates["microglial_strain"]
    c_vessel_baseline = max(0.12, min(0.92, 1.0 - 0.82 * drainage_ratio))
    c_vessel = min(0.98, c_vessel_baseline + 0.15 * strain)

    # Boundary penetration layer scaled by convective velocity
    v_eff = calculate_effective_velocity(alps_index, glymph_params, sleep_fraction=sleep_fraction)
    pe = (v_eff * geom.pvs_width) / max(1e-15, physics.d_eff_tau)
    delta_eff_um = 35.0 * (1.0 + 2.0 / (1.0 + max(0.05, pe)))

    for i in range(num_grid_points):
        r_m = r_start + i * dr
        r_um = r_m * 1e6
        r_points_um.append(round(r_um, 2))

        decay_factor = math.exp(-(r_m - r_start) / (delta_eff_um * 1e-6))
        c_val = 1.0 - (1.0 - c_vessel) * decay_factor
        c_profile.append(round(min(1.0, max(0.0, c_val)), 4))

    return {
        "alps_index": alps_index,
        "apoe_e4_alleles": apoe_e4_alleles,
        "r_points_um": r_points_um,
        "c_profile": c_profile,
        "c_vessel_boundary": round(c_vessel, 4),
        "lambda_clearance_um": round(lambda_c_um, 2),
        "delta_penetration_um": round(delta_eff_um, 2),
        "net_growth_rate": rates["effective_growth_rate"],
        "is_pathologic_expansion": rates["is_pathologic_expansion"],
    }
