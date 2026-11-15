# Neuro-Glymphatic Tau Dynamics (NGTD)

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Audit Status](https://img.shields.io/badge/Biophysical%20Audit-7%2F7%20Passed-brightgreen.svg)]()

A multiscale reaction-diffusion and perivascular transport engine modeling pathological Tau propagation in Alzheimer's disease. Couples Starling-Darcy glymphatic clearance, clinical DTI-ALPS diffusion indices, and microglial neuroinflammatory strain modulated by *APOE* $\varepsilon4$ and *TREM2*.

---

## 1. Mathematical Framework

The local concentration $c(\mathbf{x}, t)$ of toxic Tau aggregates within the cortical neuropil is modeled by:

$$\frac{\partial c}{\partial t} = \nabla \cdot (\mathbf{D}_{\text{eff}} \nabla c) + \rho \, c \, (1 - c) - k_{\text{glymph}}(\text{ALPS}) \, c - \Phi_{\text{microglia}}(c, I_{\text{strain}})$$

### Core Biophysical Couplings
1. **Perivascular Fluid Mechanics (Starling-Darcy)**:
   - Evaluates convective bulk flow along the Virchow-Robin spaces of penetrating cortical vessels.
   - Convective velocity $v_{\text{eff}}$ is parameterized directly from the clinical **DTI-ALPS index**.
2. **Macromolecular Sieving**:
   - Accounts for steric resistance across glial basement membranes ($\theta_{\text{drag}} = 1 - \sigma \approx 0.0005$).
3. **Péclet Number Phase Transition**:
   - Healthy sleep state: $Pe = 1.06 > 1.0$ (convective washout dominates).
   - Alzheimer's state: $Pe = 0.07 \ll 1.0$ (convective collapse; stasis niches emerge).
4. **Genetic & Microglial Strain Coupling**:
   - Microglial neuroinflammatory strain $I_{\text{strain}} \in [0.0, 1.0]$ scales with *APOE* $\varepsilon4$ gene dose and *TREM2* R47H carrier status.

---

## 2. Quickstart & Installation

```bash
git clone https://github.com/Lukecele/tau-glymphatic-continuum.git
cd tau-glymphatic-continuum
pip install -e .
```

### CLI Diagnostics (`tau-calc`)

```bash
# Run full 7-step biophysical invariant audit
tau-calc --audit

# Healthy individual (ALPS = 1.60, APOE-e4 = 0 copies)
tau-calc --alps 1.60 --apoe-e4 0

# High-risk Alzheimer's profile (ALPS = 1.08, APOE-e4 = 2 copies, NLRP3 hyper-activation)
tau-calc --alps 1.08 --apoe-e4 2 --nlrp3
```

---

## 3. Python API

```python
from tau_glymphatic import (
    calculate_peclet_number,
    solve_alps_critical_threshold,
    simulate_radial_pvs_profile,
)

# Solve for the analytical tipping point
thresholds = solve_alps_critical_threshold()
print(f"Critical ALPS index: {thresholds['alps_crit_kinetic']}")

# Compute radial concentration profile across perivascular space
sim = simulate_radial_pvs_profile(alps_index=1.20, apoe_e4_alleles=1)
print(f"Vessel interface concentration: {sim['c_vessel_boundary']}")
```

---

## Author & Citation

**Luca Celebrano**, MD Candidate  
School of Medicine and Surgery, University of Naples Federico II  
ORCID: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)  
GitHub: [@Lukecele](https://github.com/Lukecele)
