"""
Biophysical Models and Parameter Definitions for Neuro-Glymphatic Tau Dynamics (NGTD).

Defines physical constants, perivascular geometries, glymphatic fluid mechanics,
Tau reaction-diffusion kinetics, and microglial inflammatory strain couplings.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class PerivascularGeometry:
    """Anatomical geometry of penetrating cortical arteriole and Virchow-Robin space."""
    r_vessel: float = 10.0e-6          # Arteriole inner radius: 10 um (m)
    r_pvs_outer: float = 30.0e-6       # PVS outer boundary (astrocyte glia limitans): 30 um (m)
    r_parenchyma: float = 150.0e-6     # Interstitial parenchymal cylinder domain: 150 um (m)
    pvs_width: float = 20.0e-6         # Annular PVS width for transversal diffusion: 20 um (m)


@dataclass(frozen=True)
class BrainTissuePhysics:
    """Parenchymal porous medium transport constants."""
    porosity: float = 0.20             # Extracellular volume fraction in cortex (dimensionless)
    tortuosity: float = 1.60           # Geometric path length tortuosity (dimensionless)
    d0_tau: float = 4.5e-11            # Free molecular diffusivity of Tau in H2O at 37C (m^2/s)
    # Effective diffusivity in cortex: D_eff = D0 * porosity / (tortuosity^2)
    # 4.5e-11 * 0.20 / (1.60^2) = 3.5156e-12 m^2/s = 3.5156 um^2/s
    d_eff_tau: float = 3.5156e-12      # Effective diffusivity in parenchymal interstitium (m^2/s)


@dataclass(frozen=True)
class GlymphaticFlowParameters:
    """Glymphatic hydraulic convective transport parameters."""
    v_sleep: float = 1.87e-7           # Peak NREM sleep convective interstitial velocity: 0.187 um/s (m/s)
    v_awake: float = 3.74e-8           # Awake state convective velocity (~80% reduction): 0.0374 um/s (m/s)
    alps_healthy_ref: float = 1.60      # Reference DTI-ALPS index in healthy young/middle-aged adults
    alps_isotropic_floor: float = 1.00  # DTI-ALPS floor: 1.00 represents zero directional perivascular bias
    alps_mci_typical: float = 1.35      # Typical DTI-ALPS in Mild Cognitive Impairment (MCI)
    alps_ad_typical: float = 1.04       # Typical DTI-ALPS in severe Alzheimer's Dementia (AD)
    
    # Kedem-Katchalsky macromolecular sieving factor (1 - sigma_tau)
    # Reflects that Tau (45-60 kDa) experiences steric hindrance across the glial-basement membrane
    tau_solvent_drag_fraction: float = 0.0005  # Effective convective transmission factor (1 - sigma)


@dataclass(frozen=True)
class TauKinetics:
    """Fisher-KPP aggregation and seeding reaction constants."""
    rho_aggregation: float = 1.2e-6     # Baseline aggregation/seeding rate: ~0.104 day^-1 (s^-1)
    k_carrying_capacity: float = 1.0    # Normalized maximum toxic tau fibrillar density
    k_basal_degradation: float = 4.0e-7 # Proteasomal/autophagic non-microglial baseline clearance (s^-1)


@dataclass(frozen=True)
class MicroglialNeuroinflammation:
    """Microglial activation, TREM2/NLRP3 strain, and APOE-e4 penalty."""
    k_phagocytosis_base: float = 6.0e-7  # Basal phagocytic clearance rate of soluble tau seeds (s^-1)
    strain_threshold_flip: float = 0.50  # Threshold above which microglia become neurotoxic/pro-aggregative
    seeding_boost_rate: float = 1.2e-6   # Exosomal/inflammatory seeding amplification when strain > flip

    # Genetic risk increments to microglial strain index:
    apoe_e4_per_allele: float = 0.25     # Strain penalty per APOE-e4 allele (0, 1, or 2)
    trem2_r47h_penalty: float = 0.35     # Strain penalty for TREM2 loss-of-function variant
    nlrp3_hyperactivation: float = 0.20  # Chronic peripheral or central inflammasome activation


DEFAULT_GEOMETRY = PerivascularGeometry()
DEFAULT_PHYSICS = BrainTissuePhysics()
DEFAULT_GLYMPHATIC = GlymphaticFlowParameters()
DEFAULT_TAU_KINETICS = TauKinetics()
DEFAULT_MICROGLIA = MicroglialNeuroinflammation()
