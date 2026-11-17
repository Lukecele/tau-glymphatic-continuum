"""
Brinkman-Stokes Hydrodynamics and Nucleation-Elongation PDE Engine.

Implements rigorous fluid mechanics for annular perivascular spaces (PVS) coupled to
Brinkman-extended porous neuropil flow and Knowles-Vendruscolo bi-species
(monomer-fibril) nucleated polymerization kinetics.
"""

import math
from dataclasses import dataclass
from typing import Dict, Any, List, Tuple
from tau_glymphatic.models import (
    DEFAULT_GEOMETRY,
    DEFAULT_PHYSICS,
    DEFAULT_GLYMPHATIC,
    PerivascularGeometry,
    BrainTissuePhysics,
    GlymphaticFlowParameters,
)


@dataclass(frozen=True)
class BrinkmanHydrodynamicParameters:
    """Exact hydrodynamic parameters for fluid and tissue phase."""
    viscosity_csf: float = 0.89e-3        # Dynamic viscosity of CSF at 37C: 0.89 mPa*s (Pa*s)
    permeability_parenchyma: float = 1.8e-15 # Intrinsic hydraulic permeability of cortex K_p (m^2)
    beavers_joseph_slip: float = 0.10     # Dimensionless BJS slip coefficient at glia limitans
    arterial_pulsatility_freq: float = 1.0 # Heart rate frequency 1 Hz (60 bpm)
    arterial_wall_amplitude: float = 0.05  # Radial displacement amplitude fraction (5%)
    pvs_hydraulic_attenuation: float = 0.0615 # Attenuation factor eta_pvs accounting for ECM meshwork & pulsatile duty cycle
    effective_pressure_gradient_pa_m: float = 2.46 # Effective driving gradient |dp/dz|_eff = eta_pvs * 40 Pa/m


@dataclass(frozen=True)
class NucleationElongationKinetics:
    """Knowles-Vendruscolo master equation parameters for tau fibrillization."""
    k_primary_nucleation: float = 1.0e-5  # Primary nucleation constant k_n (M^-1 s^-1)
    k_elongation: float = 2.0e4           # Monomer addition to fibril ends k_e (M^-1 s^-1)
    k_secondary_nucleation: float = 1.5e2 # Surface-catalyzed fibril nucleation k_2 (M^-2 s^-1)
    monomer_steady_conc: float = 2.0e-8   # Basal soluble tau monomer concentration (20 nM = 2e-8 M)


DEFAULT_BRINKMAN = BrinkmanHydrodynamicParameters()
DEFAULT_NUCLEATION = NucleationElongationKinetics()


def get_annular_stokes_peak_radius(geom: PerivascularGeometry = DEFAULT_GEOMETRY) -> float:
    """
    Computes exact radial coordinate r* of maximum axial velocity for Stokes flow
    in an annular duct (r0 <= r <= r1):
      r* = sqrt( (r1^2 - r0^2) / (2 * ln(r1 / r0)) )
    
    Due to logarithmic curvature in cylindrical coordinates, r* falls strictly
    inside the arithmetic midpoint (r0 + r1) / 2.
    For r0 = 10.0 um, r1 = 30.0 um: r* = sqrt(800 / (2 * ln(3))) = 19.0813 um.
    """
    r0 = geom.r_vessel
    r1 = geom.r_pvs_outer
    return math.sqrt((r1 * r1 - r0 * r0) / (2.0 * math.log(r1 / r0)))


def solve_annular_stokes_velocity_profile(
    r_points: List[float],
    dp_dz: float = -20.0,  # Axial pressure gradient in Pa/m
    geom: PerivascularGeometry = DEFAULT_GEOMETRY,
    brinkman: BrinkmanHydrodynamicParameters = DEFAULT_BRINKMAN,
    attenuation_factor: float = 1.0,  # Set to brinkman.pvs_hydraulic_attenuation (0.0615) for in vivo effective velocity
) -> List[float]:
    """
    Solves exact analytical velocity profile v_z(r) for steady laminar Stokes flow
    in an annular cylinder (r_0 <= r <= r_1):
    
      v_z(r) = (eta / 4*mu) * (-dp/dz) * [ (r_1^2 - r^2) + (r_1^2 - r_0^2) * ln(r/r_1) / ln(r_0/r_1) ]
    
    Returns axial velocity in m/s at each radial coordinate r.
    Peak velocity occurs exactly at r* = sqrt((r1^2 - r0^2)/(2*ln(r1/r0))) = 19.081 um.
    """
    r0 = geom.r_vessel
    r1 = geom.r_pvs_outer
    mu = brinkman.viscosity_csf
    grad_p = -dp_dz * attenuation_factor  # Positive effective drive

    ln_ratio = math.log(r0 / r1)
    r1_sq = r1 * r1
    r0_sq = r0 * r0

    velocities: List[float] = []
    for r in r_points:
        if r < r0 or r > r1:
            velocities.append(0.0)
            continue
        term1 = r1_sq - r * r
        term2 = (r1_sq - r0_sq) * (math.log(r / r1) / ln_ratio)
        vz = (grad_p / (4.0 * mu)) * (term1 - term2)
        velocities.append(max(0.0, vz))

    return velocities


def calculate_front_velocity_fisher(
    d_eff: float,
    effective_growth_rate: float,
) -> float:
    """
    Computes exact Kolmogorov-Petrovsky-Piskunov minimal traveling wave speed:
      v* = 2 * sqrt(D_eff * k_growth)  for k_growth > 0
    Returns 0.0 if growth rate is negative (clearance regime).
    """
    if effective_growth_rate <= 0.0:
        return 0.0
    return 2.0 * math.sqrt(d_eff * effective_growth_rate)


def evaluate_wave_pinning_criterion(
    alps_index: float,
    v_eff_perivascular: float,
    d_eff: float,
    effective_growth_rate: float,
) -> Dict[str, Any]:
    """
    Evaluates hydrodynamic wave-pinning of the Tau spreading front.
    
    The advancing seeding front has velocity v*.
    The opposing convective perivascular drainage has velocity v_pvs.
    
    If v_pvs >= v*: Front is PINNED and WASHED OUT (Clearance Stability).
    If v_pvs < v*: Front PROPAGATES into neocortex (Braak Expansion).
    """
    v_front = calculate_front_velocity_fisher(d_eff, effective_growth_rate)
    
    if v_front == 0.0:
        front_peclet = float("inf")
        is_pinned = True
        status = "COMPLETE EXTINCTION (Negative Growth Rate)"
    else:
        front_peclet = v_eff_perivascular / v_front
        is_pinned = front_peclet >= 1.0
        status = "WAVE PINNED (Convection >= Seeding Speed)" if is_pinned else "WAVE PROPAGATING (Seeding Speed > Convection)"

    return {
        "alps_index": alps_index,
        "v_front_seeding_nm_s": round(v_front * 1e9, 4),
        "v_convective_washout_nm_s": round(v_eff_perivascular * 1e9, 4),
        "front_peclet_number": round(front_peclet, 4) if front_peclet != float("inf") else 999.0,
        "is_pinned": is_pinned,
        "status": status,
    }


def calculate_dimensionless_hydrodynamic_numbers(
    v_velocity: float = 0.187e-6,  # 0.187 um/s in healthy sleep
    d_eff: float = 3.5156e-12,     # 3.5156 um^2/s
    k_growth: float = 1.0e-6,      # 1e-6 s^-1
    rho_csf: float = 1000.0,       # kg/m^3
    geom: PerivascularGeometry = DEFAULT_GEOMETRY,
    brinkman: BrinkmanHydrodynamicParameters = DEFAULT_BRINKMAN,
) -> Dict[str, Any]:
    """
    Computes all fundamental dimensionless groups in fluid mechanics and transport theory
    to rigorously demonstrate physical consistency against theoretical physics criteria:
    
    1. Reynolds Number (Re = rho * v * D_h / mu): Verifies inertia-free Stokes regime.
    2. Womersley Number (alpha = w * sqrt(omega * rho / mu)): Verifies quasi-steady pulsatile flow.
    3. Annular Peclet Number (Pe = v * w / D_eff): Verifies advection-diffusion balance.
    4. Schmidt Number (Sc = nu / D_eff): Verifies momentum vs mass boundary layers.
    5. Damkohler Number (Da = k_growth * w / v): Verifies reaction vs advection timescale.
    6. Glia Limitans Wall Shear Stress (tau_w): Verifies mechanosensory astrocytic regime.
    """
    w_pvs = geom.pvs_width
    d_hydraulic = 2.0 * w_pvs
    mu = brinkman.viscosity_csf
    nu = mu / rho_csf  # Kinematic viscosity (m^2/s)

    # 1. Reynolds number
    reynolds = (rho_csf * abs(v_velocity) * d_hydraulic) / mu

    # 2. Womersley number (pulsatile frequency f = 1 Hz -> omega = 2*pi)
    omega = 2.0 * math.pi * brinkman.arterial_pulsatility_freq
    womersley = w_pvs * math.sqrt((omega * rho_csf) / mu)

    # 3. Peclet number
    peclet = (abs(v_velocity) * w_pvs) / d_eff

    # 4. Schmidt number
    schmidt = nu / d_eff

    # 5. Damkohler number
    # Local PVS width scale:
    damkohler_local = (k_growth * w_pvs) / abs(v_velocity) if abs(v_velocity) > 1e-15 else float("inf")
    # Penetrating cortical vessel depth scale (L = 1.0 mm, elongation kinetic timescale):
    l_penetrating = 1.0e-3
    k_elongation = 2.0e-5  # Knowles-Vendruscolo elongation rate
    damkohler_penetrating = (k_elongation * l_penetrating) / abs(v_velocity) if abs(v_velocity) > 1e-15 else float("inf")

    # 6. Wall Shear Stress at glia limitans (tau = mu * dv/dr approx mu * (2 * v_avg / w_pvs))
    tau_wall = mu * (2.0 * abs(v_velocity) / w_pvs)

    return {
        "reynolds_number": reynolds,
        "is_stokes_regime_valid": reynolds < 1e-3,
        "womersley_number": round(womersley, 5),
        "is_quasi_steady_valid": womersley < 0.10,
        "annular_peclet_number": round(peclet, 4),
        "schmidt_number": round(schmidt, 1),
        "damkohler_number": round(damkohler_local, 6) if damkohler_local != float("inf") else 999.0,
        "damkohler_number_local": round(damkohler_local, 6) if damkohler_local != float("inf") else 999.0,
        "damkohler_number_penetrating": round(damkohler_penetrating, 4) if damkohler_penetrating != float("inf") else 999.0,
        "wall_shear_stress_pascal": tau_wall,
        "wall_shear_stress_millipascal": round(tau_wall * 1e3, 4),
        "mechanosensory_piezo1_sensing_valid": 1e-6 <= tau_wall <= 1e-1,
        "summary": (
            f"Stokes valid (Re={reynolds:.2e} << 1); Quasi-steady valid (alpha={womersley:.3f} << 1); "
            f"Momentum dominates mass (Sc={schmidt:.1e} >> 1)."
        ),
    }

