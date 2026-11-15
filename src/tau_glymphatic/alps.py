"""
DTI-ALPS Index Mapping, Hydraulic Conductivity, and Péclet Analysis.

Translates human clinical Diffusion Tensor Imaging Along the Perivascular Space (DTI-ALPS)
measurements into biophysical interstitial fluid flow velocities and calculates
the critical Glymphatic Collapse Threshold (ALPS_crit).
"""

from typing import Dict, Any
from tau_glymphatic.models import (
    DEFAULT_GEOMETRY,
    DEFAULT_PHYSICS,
    DEFAULT_GLYMPHATIC,
    DEFAULT_TAU_KINETICS,
    PerivascularGeometry,
    BrainTissuePhysics,
    GlymphaticFlowParameters,
    TauKinetics,
)


def calculate_effective_velocity(
    alps_index: float,
    glymph_params: GlymphaticFlowParameters = DEFAULT_GLYMPHATIC,
    sleep_fraction: float = 0.33,  # 8 hours sleep per 24 hours
) -> float:
    """
    Calculates the 24-hour time-averaged effective perivascular convective velocity (m/s)
    from a patient's clinical DTI-ALPS index.
    """
    if alps_index <= glymph_params.alps_isotropic_floor:
        return 0.0

    scaling_ratio = (alps_index - glymph_params.alps_isotropic_floor) / (
        glymph_params.alps_healthy_ref - glymph_params.alps_isotropic_floor
    )
    scaling_ratio = max(0.0, min(1.5, scaling_ratio))

    v_sleep_scaled = glymph_params.v_sleep * scaling_ratio
    v_awake_scaled = glymph_params.v_awake * scaling_ratio

    v_eff = (sleep_fraction * v_sleep_scaled) + ((1.0 - sleep_fraction) * v_awake_scaled)
    return v_eff


def calculate_peclet_number(
    alps_index: float,
    geom: PerivascularGeometry = DEFAULT_GEOMETRY,
    physics: BrainTissuePhysics = DEFAULT_PHYSICS,
    glymph_params: GlymphaticFlowParameters = DEFAULT_GLYMPHATIC,
    peak_sleep: bool = True,
) -> float:
    """
    Computes the dimensionless Péclet number Pe = v * L / D_eff for Tau in the PVS.
    """
    if peak_sleep:
        scaling_ratio = max(
            0.0,
            (alps_index - glymph_params.alps_isotropic_floor)
            / (glymph_params.alps_healthy_ref - glymph_params.alps_isotropic_floor),
        )
        velocity = glymph_params.v_sleep * scaling_ratio
    else:
        velocity = calculate_effective_velocity(alps_index, glymph_params)

    l_pvs = geom.pvs_width
    d_eff = physics.d_eff_tau
    return (velocity * l_pvs) / d_eff


def solve_alps_critical_threshold(
    geom: PerivascularGeometry = DEFAULT_GEOMETRY,
    physics: BrainTissuePhysics = DEFAULT_PHYSICS,
    glymph_params: GlymphaticFlowParameters = DEFAULT_GLYMPHATIC,
    tau_kinetics: TauKinetics = DEFAULT_TAU_KINETICS,
) -> Dict[str, float]:
    """
    Solves for the analytical critical ALPS index (ALPS_crit) where convective
    clearance exactly balances the intrinsic tau auto-aggregation rate.
    """
    l_pvs = geom.pvs_width
    d_eff = physics.d_eff_tau
    theta = glymph_params.tau_solvent_drag_fraction

    # k_glymph = theta * (v_eff / L_pvs) = rho_aggregation
    # v_needed = rho_aggregation * L_pvs / theta
    v_needed = (tau_kinetics.rho_aggregation * l_pvs) / theta

    # Ratio required of 24h effective velocity
    # Under typical 8h sleep schedule, v_eff_ref is:
    v_eff_ref = calculate_effective_velocity(glymph_params.alps_healthy_ref, glymph_params)
    ratio_needed = v_needed / v_eff_ref
    ratio_needed = max(0.0, min(1.0, ratio_needed))

    alps_crit = glymph_params.alps_isotropic_floor + ratio_needed * (
        glymph_params.alps_healthy_ref - glymph_params.alps_isotropic_floor
    )

    # ALPS at Pe = 1.0 (peak sleep)
    v_pe1 = d_eff / l_pvs
    ratio_pe1 = v_pe1 / glymph_params.v_sleep
    alps_pe1 = glymph_params.alps_isotropic_floor + ratio_pe1 * (
        glymph_params.alps_healthy_ref - glymph_params.alps_isotropic_floor
    )

    return {
        "alps_crit_kinetic": round(alps_crit, 4),
        "alps_pe1_convective": round(alps_pe1, 4),
        "v_needed_nm_per_s": round(v_needed * 1e9, 2),
        "peclet_healthy": round(calculate_peclet_number(glymph_params.alps_healthy_ref, peak_sleep=True), 4),
        "peclet_ad": round(calculate_peclet_number(glymph_params.alps_ad_typical, peak_sleep=True), 4),
    }
