"""
Neuro-Glymphatic Tau Dynamics (NGTD).

A Multiscale Biophysical Reaction-Diffusion and Perivascular Clearance Engine
for Alzheimer's Disease and Neurodegenerative Tauopathies.
"""

from tau_glymphatic.models import (
    PerivascularGeometry,
    BrainTissuePhysics,
    GlymphaticFlowParameters,
    TauKinetics,
    MicroglialNeuroinflammation,
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
from tau_glymphatic.audit import run_full_biophysical_audit
from tau_glymphatic.cohorts import (
    ClinicalCohortMetadata,
    AUDITED_COHORTS,
    CohortCompatibilityAuditor,
)

__version__ = "1.1.1"
__author__ = "Luca Celebrano"
__license__ = "CC-BY-4.0"

__all__ = [
    "PerivascularGeometry",
    "BrainTissuePhysics",
    "GlymphaticFlowParameters",
    "TauKinetics",
    "MicroglialNeuroinflammation",
    "DEFAULT_GEOMETRY",
    "DEFAULT_PHYSICS",
    "DEFAULT_GLYMPHATIC",
    "DEFAULT_TAU_KINETICS",
    "DEFAULT_MICROGLIA",
    "calculate_effective_velocity",
    "calculate_peclet_number",
    "solve_alps_critical_threshold",
    "calculate_microglial_strain",
    "calculate_net_clearance_rate",
    "simulate_radial_pvs_profile",
    "run_full_biophysical_audit",
    "ClinicalCohortMetadata",
    "AUDITED_COHORTS",
    "CohortCompatibilityAuditor",
]
