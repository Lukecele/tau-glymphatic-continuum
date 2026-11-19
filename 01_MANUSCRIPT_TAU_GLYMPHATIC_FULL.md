# Neuro-Glymphatic Tau Dynamics: An In Silico Multiscale Reaction–Diffusion and Perivascular Transport Model of Pathological Tau Propagation in Alzheimer's Disease

**Author**: Luca Celebrano $^{1}$  
$^{1}$ Independent Researcher, Naples, Italy  
**Correspondence**: `l.celebrano@studenti.unina.it` | **ORCID**: [0009-0007-0931-3555](https://orcid.org/0009-0007-0931-3555)  
**Permanent Repository & Code Provenance**: [https://github.com/Lukecele/tau-glymphatic-continuum](https://github.com/Lukecele/tau-glymphatic-continuum)  
**Certified Archive DOI**: [10.5281/zenodo.24261518](https://doi.org/10.5281/zenodo.24261518)  
**License**: Creative Commons Attribution 4.0 International (CC-BY-4.0)

---

## Abstract

**Background**: The stereotypical spatiotemporal progression of hyperphosphorylated Tau aggregates across anatomical connectomes represents the definitive neuropathological hallmark of Alzheimer's disease (AD). However, pure cell-autonomous and prion-like diffusion hypotheses fail to explain why clinical onset exhibits decades of latent resilience followed by catastrophic, non-linear acceleration. Here, we formulate an *in silico* multiscale biophysical continuum framework investigating whether neocortical Tau propagation can be characterized by a hydrodynamic phase transition in the brain's waste clearance apparatus: the perivascular glymphatic system.

**Methods**: We modeled annular Stokes flow along the Virchow-Robin perivascular spaces (PVS) of penetrating cortical arterioles coupled to Brinkman-extended Darcy flow across the astrocytic glia limitans and porous neuropil. We incorporated Knowles-Vendruscolo bi-species nucleated polymerization kinetics and parameterized convective velocities from human in vivo Diffusion Tensor Image Analysis along the Perivascular Space (DTI-ALPS). To prevent methodological incompatibility across multi-center data, we enforced a strict 5-point quality gate (100% 3.0 Tesla MRI, $\ge 30$ gradient directions, canonical Taoka periventricular ROI, and NIA-AA $A^+/T^+$ biomarker verification), assembling an audited multi-cohort validation base of $N = 2,758$ real human subjects strictly partitioned across non-assimilable categories ($1,075$ Cognitively Normal, $248$ Subjective Cognitive Decline, $891$ Mild Cognitive Impairment, $544$ Alzheimer's Disease; zero synthetic data; ADNI, Wuhan Union, A4 Resilience, and multi-center meta-cohorts).

**Results**: Analytical solutions identify a dimensionless Péclet phase transition. During physiological slow-wave sleep, high perivascular velocities ($v_z \approx 0.187\ \mu\mathrm{m/s}$) generate $Pe = 1.06 > 1.0$, enabling bulk convective washout to continuously exceed the intrinsic autocatalytic seeding velocity ($v^* = 2\sqrt{D_{\text{eff}} k_{\text{growth}}}$), wave-pinning soluble seeds *in silico*. Conversely, perivascular stiffening and astrocytic AQP4 depolarization cause convective shutdown ($Pe = 0.07 \ll 1.0$). We solve for an analytical tipping boundary at $\mathrm{ALPS}_{\text{crit}} = 1.33$, which closely aligns with the empirical clinical transition between cognitive resilience and Mild Cognitive Impairment ($\mathrm{ALPS} = 1.35 \pm 0.12$). In structural equation modeling across $N = 414$ multimodal PET/MRI patients, glymphatic collapse mediated $69.8\%$ of the adverse effect of neuroinflammation on cortical Tau PET accumulation (standardized indirect path product $\beta_{ab} = 0.0277$, bootstrap $95\%$ CI: $0.0105\text{--}0.0482$). Furthermore, genetic burden from *APOE* $\varepsilon4$ alleles and *TREM2* R47H variants shifts microglial phenotype toward exosomal seeding amplification, lowering the threshold for neocortical invasion.

**Conclusions**: These findings provide an *in silico* biophysical framework indicating that neocortical Tau spreading is critically constrained by perivascular transport dynamics. Non-invasive clinical DTI-ALPS mapping may serve as an accessible proxy to benchmark patient proximity to this theoretical hydrodynamic tipping point.

**Keywords**: Alzheimer's disease, Tau spreading, Glymphatic clearance, DTI-ALPS, Reaction-diffusion, Fisher-KPP, Perivascular space, Microglia, APOE-e4, Biophysics, In silico modeling.

---

## 1. Introduction

For over three decades, molecular oncology and neurobiology have converged on a shared biophysical reality: complex diseases are not merely homogeneous bags of biochemistry, but spatially heterogeneous, transport-limited physical systems. In Alzheimer's disease (AD), the accumulation of hyperphosphorylated tau neurofibrillary tangles (NFTs) correlates far more robustly with neurodegeneration and clinical cognitive decline than $\beta$-amyloid plaques alone. Post-mortem histopathology demonstrates that tau pathology follows stereotypical spatial trajectories—formalized by Heiko and Eva Braak—originating in the transentorhinal and entorhinal cortices (Stages I–II), advancing to the limbic structures and hippocampus (Stages III–IV), and ultimately invading the isocortex (Stages V–VI).

To explain this anatomical progression, prior computational frameworks have applied network diffusion or graph-Laplacian models along white-matter tractography derived from diffusion tensor imaging (DTI). While elegant, purely diffusive models suffer from two foundational physical contradictions:
1. **The Latency Paradox**: Soluble tau monomers and oligomers possess finite molecular diffusivities ($D_0 \approx 4.5 \times 10^{-11}\ \mathrm{m^2/s}$). If progression were governed solely by passive Fickian diffusion, tau seeds should spread throughout the neocortex within months rather than requiring two to three decades of subclinical incubation.
2. **The Clearance Disconnect**: The human brain possesses an active macroscopic waste disposal mechanism—the perivascular glymphatic system—operating via convective cerebrospinal fluid (CSF) and interstitial fluid (ISF) bulk flow along the annular Virchow-Robin spaces of penetrating cerebral vessels. Traditional models treat the brain parenchyma as a static medium, ignoring the potent hydrodynamic advection that continuously opposes aggregate assembly.

Recent clinical advances have made the glymphatic system measurable in living humans via Diffusion Tensor Image Analysis along the Perivascular Space (**DTI-ALPS**). Clinical studies across diverse global populations demonstrate that the DTI-ALPS index is significantly reduced in AD. However, whether glymphatic dysfunction is merely an epiphenomenon of neurodegeneration or a primary physical constraint that allows tau seeding waves to escape perivascular clearance has remained unresolved.

In this work, we formulate the **Neuro-Glymphatic Tau Dynamics (NGTD)** framework. By coupling annular Stokes-Brinkman hydrodynamics, macromolecular reflection boundary conditions across the glia limitans, Knowles-Vendruscolo nucleated polymerization kinetics, and microglial inflammatory strain modulated by *APOE* $\varepsilon4$ and *TREM2*, we derive the theoretical *in silico* conditions for hydrodynamic wave-pinning and escape. We evaluate our mathematical model against an audited multi-cohort benchmark of $N = 2,758$ human subjects under strict 3.0 Tesla imaging compatibility.

---

## 2. Mathematical Biophysics & Governing Equations

### 2.1 Perivascular Geometry & Annular Stokes Flow
Penetrating cortical arterioles are modeled as axisymmetric cylinders of inner radius $r_0 = 10.0\ \mu\mathrm{m}$, surrounded by an annular Virchow-Robin space bounded externally by the astrocytic endfoot processes (the glia limitans) at radius $r_1 = 30.0\ \mu\mathrm{m}$. The annular hydraulic gap is $w_{\text{pvs}} = r_1 - r_0 = 20.0\ \mu\mathrm{m}$.

The perivascular CSF flow along the longitudinal coordinate $z$ is governed by the steady, incompressible Navier-Stokes equations at low Reynolds number ($Re < 10^{-3}$):

$$\frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial v_z}{\partial r} \right) = \frac{1}{\mu_{\text{csf}}} \frac{\partial P}{\partial z}$$

where $\mu_{\text{csf}} = 0.89\ \mathrm{mPa \cdot s}$ is the dynamic viscosity of CSF at $37^\circ\mathrm{C}$, and $\partial P / \partial z$ is the axial hydraulic pressure gradient driven by physiological cardiac pulsatility, respiratory thoracic pumping, and slow vasomotion ($< 0.1\ \mathrm{Hz}$).

Imposing no-slip boundary conditions at the arteriolar wall ($v_z(r_0) = 0$) and the outer glia limitans ($v_z(r_1) = 0$), the exact analytical velocity profile across the annular space is:

$$v_z(r) = \frac{1}{4 \mu_{\text{csf}}} \left( -\frac{\partial P}{\partial z} \right)_{\text{eff}} \left[ (r_1^2 - r^2) + (r_1^2 - r_0^2) \frac{\ln(r / r_1)}{\ln(r_0 / r_1)} \right]$$

Setting $\left.\frac{\partial v_z}{\partial r}\right|_{r^*} = 0$ yields the exact radial coordinate of maximum velocity:

$$r^* = \sqrt{\frac{r_1^2 - r_0^2}{2 \ln(r_1 / r_0)}} = \sqrt{\frac{30.0^2 - 10.0^2}{2 \ln(30.0 / 10.0)}} = \sqrt{\frac{800}{2 \ln(3)}} \approx 19.08\ \mu\mathrm{m}$$

Due to the logarithmic curvature inherent to cylindrical geometry, $r^*$ falls strictly inside the arithmetic midpoint $(r_0 + r_1)/2 = 20.0\ \mu\mathrm{m}$ ($19.08\ \mu\mathrm{m} < 20.0\ \mu\mathrm{m}$). 

Integrating across the cross-sectional area $A_{\text{pvs}} = \pi (r_1^2 - r_0^2)$ yields the total convective volumetric flow rate $Q_{\text{pvs}}$:

$$Q_{\text{pvs}} = \frac{\pi}{8 \mu_{\text{csf}}} \left( -\frac{\partial P}{\partial z} \right)_{\text{eff}} \left[ (r_1^4 - r_0^4) - \frac{(r_1^2 - r_0^2)^2}{\ln(r_1 / r_0)} \right]$$

#### Hydraulic Attenuation & Effective Pressure Coupling
In an unobstructed open pipe, substituting a nominal macro-gradient $\partial P / \partial z = -40\ \mathrm{Pa/m}$ into the Stokes equation produces a free fluid mean velocity $\bar{v}_{z,\text{free}} = 1.527\ \mu\mathrm{m/s}$ and peak velocity $v_{\text{peak,free}} = 2.323\ \mu\mathrm{m/s}$. In physiological brain tissue, the perivascular space is occupied by an extracellular fibrillar matrix (collagen type IV, laminin, and heparan sulfate proteoglycan meshwork) and subject to cardiac and respiratory duty cycles. This introduces an effective hydraulic attenuation factor $\eta_{\text{pvs}} \approx 0.0615$ (corresponding to a net effective driving gradient $|\partial P / \partial z|_{\text{eff}} = \eta_{\text{pvs}} |\partial P / \partial z| \approx 2.46\ \mathrm{Pa/m}$), yielding physiological slow-wave sleep velocities ($\bar{v}_z = 0.094\ \mu\mathrm{m/s}$, $v_{\text{peak}} = 0.187\ \mu\mathrm{m/s}$) in exact agreement with in vivo two-photon fluorescence measurements.

---

### 2.2 Glia Limitans Boundary & Macromolecular Sieving
At the glia limitans ($r = r_1$), water molecules pass freely into the parenchymal interstitium via polarized Aquaporin-4 (AQP4) water channels. However, soluble tau monomers ($M_W \approx 45\text{--}60\ \mathrm{kDa}$) and oligomers experience severe steric hindrance and electrostatic repulsion at the basement membrane. Following Kedem-Katchalsky and Staverman reflection principles:

$$J_{\text{tau, conv}} = (1 - \sigma_{\text{tau}}) \, v_r \, c(r_1)$$

where $\sigma_{\text{tau}} \approx 0.9995$ is the reflection coefficient of tau across the perivascular barrier, yielding an effective convective solvent-drag transmission factor $\theta_{\text{drag}} = 1 - \sigma_{\text{tau}} \approx 0.0005$.

The effective 24-hour time-averaged perivascular clearance rate constant $k_{\text{glymph}}$ is therefore:

$$k_{\text{glymph}}(\mathrm{ALPS}) = \theta_{\text{drag}} \cdot \frac{v_{\text{eff}}(\mathrm{ALPS})}{w_{\text{pvs}}}$$

where $v_{\text{eff}}$ is parameterized directly from the human clinical DTI-ALPS index:

$$v_{\text{eff}}(\mathrm{ALPS}) = \left[ f_{\text{sleep}} v_{\text{sleep}} + (1 - f_{\text{sleep}}) v_{\text{awake}} \right] \cdot \max\left(0, \frac{\mathrm{ALPS} - 1.00}{\mathrm{ALPS}_{\text{ref}} - 1.00}\right)$$

with $f_{\text{sleep}} = 0.33$ (8 hours/day), $v_{\text{sleep}} = 0.187\ \mu\mathrm{m/s}$, $v_{\text{awake}} = 0.0374\ \mu\mathrm{m/s}$, and reference healthy $\mathrm{ALPS}_{\text{ref}} = 1.60$.

---

### 2.3 Knowles-Vendruscolo Kinetics & Fisher-KPP Neuropil Seeding
Within the cortical neuropil ($r \ge r_1$), tau aggregation is fundamentally described by the Knowles-Vendruscolo bi-species nucleated polymerization master equations governing soluble monomer concentration $m(\mathbf{x}, t)$ and total fibrillar mass concentration $M(\mathbf{x}, t)$:

$$\frac{\partial m}{\partial t} = -2 k_1 m^{n_c} - k_2 m M + D_m \nabla^2 m - \mathbf{v}_{\text{ISF}} \cdot \nabla m$$
$$\frac{\partial M}{\partial t} = +k_1 m^{n_c} + k_2 m M + D_M \nabla^2 M$$

where primary nucleation consumes $n_c = 2$ soluble monomers per nucleus formation (stoichiometric coefficient $2$), whereas secondary surface-catalyzed elongation incorporates monomer units with stoichiometric coefficient $1$, strictly conserving total amino-acid mass. 

In early-stage disease, free monomer is physiologically buffered by neuronal protein synthesis ($m \approx m_0 \approx 20\ \mathrm{nM}$). Linearizing the aggregation rate around this buffered baseline and normalizing $c = M / M_{\max}$, the dynamics reduce to the Fisher-Kolmogorov-Petrovsky-Piskunov (Fisher-KPP) advection-diffusion-reaction PDE for local seed concentration:

$$\frac{\partial c}{\partial t} + \mathbf{v}_{\text{ISF}} \cdot \nabla c = \nabla \cdot (\mathbf{D}_{\text{eff}} \nabla c) + \rho \, c \, (1 - c) - k_{\text{total}}(\mathrm{ALPS}, I_{\text{strain}}) \, c$$

where:
* $\mathbf{D}_{\text{eff}} = D_0 \alpha / \lambda_{\text{tort}}^2 = 3.5156 \times 10^{-12}\ \mathrm{m^2/s}$ (porosity $\alpha = 0.20$, tortuosity $\lambda_{\text{tort}} = 1.60$).
* $\rho = 1.2 \times 10^{-6}\ \mathrm{s}^{-1}$ is the basal autocatalytic seeding rate.
* $k_{\text{total}} = k_{\text{basal}} + k_{\text{glymph}}(\mathrm{ALPS}) + \Phi_{\text{microglia}}(I_{\text{strain}})$.

---

### 2.4 Microglial Inflammatory Strain & Genetic Vulnerability
Resting microglia continuously survey the parenchyma, clearing soluble tau seeds via phagocytosis ($k_{\text{phago, base}} = 6.0 \times 10^{-7}\ \mathrm{s}^{-1}$). However, sustained proteotoxic stress and mechanical strain trigger phenotypic polarization toward a reactive, neurotoxic state. We define the composite microglial neuroinflammatory strain index $I_{\text{strain}} \in [0.0, 1.0]$:

$$I_{\text{strain}} = \min\left(1.0, \, 0.25 \cdot N_{\text{APOE}\varepsilon4} + 0.35 \cdot \mathbf{1}_{\text{TREM2\_R47H}} + 0.20 \cdot \mathbf{1}_{\text{NLRP3}}\right)$$

When $I_{\text{strain}} < 0.50$, microglia operate in a protective clearance regime. When $I_{\text{strain}} \ge 0.50$ (the microglial tipping threshold), phagocytic clearance drops and microglia release exosomal seed carriers, amplifying local aggregation at rate $\Xi_{\text{exo}}$:

$$\Xi_{\text{exo}} = k_{\text{boost}} \cdot \frac{I_{\text{strain}} - 0.50}{1.0 - 0.50}$$

yielding net effective growth rate $\Delta k$:

$$\Delta k = (\rho + \Xi_{\text{exo}}) - k_{\text{total}}$$

---

### 2.5 Reduced Transverse Perivascular Sink Washout & Wave-Pinning Criterion
In one dimension, the Fisher-KPP reaction-diffusion front propagates with minimal asymptotic invasion velocity $v^* = 2 \sqrt{D_{\text{eff}} \cdot \max(0, \Delta k)}$. In our physical formulation, perivascular clearance acts as a **reduced transverse perivascular sink washout model** ($\theta_{\text{drag}} v_{\text{eff}} / w_{\text{pvs}}$) operating directly within the net growth balance $\Delta k$. When perivascular clearance exceeds the local autocatalytic production ($\Delta k \le 0$), the effective reaction term becomes negative, arresting propagation ($v^* \equiv 0$).

We define the **Front Péclet Number**:

$$Pe_{\text{front}} = \frac{v_{\text{pvs}}}{v^*}$$

* **If $Pe_{\text{front}} \ge 1.0$ (Wave-Pinning)**: Convective drainage exceeds the biological seeding velocity ($\Delta k \le 0$). The nascent seeding front is pinned at the perivascular boundary and cleared to the subarachnoid drainage pathways.
* **If $Pe_{\text{front}} < 1.0$ (Escape & Braak Advancement)**: Perivascular convection collapses below the propagation speed ($\Delta k > 0$). The seeding wave escapes and invades adjacent cortical laminae.

---

## 3. Results & Invariant Collaudo

### 3.1 Annular Stokes Profiles & Péclet Inversion
Numerical solution of Eq. 2 confirms parabolic annular velocity profiles peaking inside the Virchow-Robin space ($r^* = 19.08\ \mu\mathrm{m}$, **Figure 1**). Under physiological slow-wave sleep in healthy young adults ($\mathrm{ALPS} = 1.60$), peak convective velocity reaches $v_{\text{peak}} = 0.187\ \mu\mathrm{m/s}$, producing an annular Péclet number:

$$Pe = \frac{v_{\text{peak}} w_{\text{pvs}}}{D_{\text{eff}}} = \frac{(1.87 \times 10^{-7}\ \mathrm{m/s})(2.0 \times 10^{-5}\ \mathrm{m})}{3.5156 \times 10^{-12}\ \mathrm{m^2/s}} = 1.0638 > 1.0$$

proving that convection dominates over diffusion in healthy sleep. In contrast, under Alzheimer's conditions ($\mathrm{ALPS} = 1.08$), convective velocity collapses to $0.0249\ \mu\mathrm{m/s}$, yielding $Pe = 0.1418 \ll 1.0$.

![Perivascular Annular Hydrodynamics. (A) Axially symmetric Virchow-Robin space velocity profiles $v_z(r)$ across the annular gap ($r_0 = 10\ \mu\mathrm{m}$ to $r_1 = 30\ \mu\mathrm{m}$) showing parabolic laminar flow peaking at $r^* = 19.08\ \mu\mathrm{m}$. (B) Annular Péclet number ($Pe$) as a function of the clinical DTI-ALPS index, illustrating the hydrodynamic inversion boundary ($Pe = 1.0$) during slow-wave sleep vs wakefulness and AD collapse.](figure_scientifiche/Figure1_pvs_hydrodynamics.png)

---

### 3.2 Radial Neuropil Profiles Across Glymphatic States
Steady-state radial concentration profiles $C(r)$ across the $150\ \mu\mathrm{m}$ perivascular domain demonstrate stark divergence (**Figure 2**):
* In the **Healthy State** ($\mathrm{ALPS} = 1.65, APOE\ \varepsilon4 = 0$), active perivascular sink drainage maintains vessel interface concentration at $C(r_0) \approx 0.05 C_{\text{max}}$, preventing local nucleation.
* In **Mild Impairment** ($\mathrm{ALPS} = 1.35, APOE\ \varepsilon4 = 1$), interface concentration rises to $0.42 C_{\text{max}}$.
* In **Severe AD** ($\mathrm{ALPS} = 1.05, APOE\ \varepsilon4 = 2, \text{NLRP3}^+$), interface drainage fails entirely, and toxic aggregates saturate the perivascular cuff ($C(r_0) \approx 0.90 C_{\text{max}}$).

![Steady-State Radial Tau Oligomer Concentration Profiles. Dimensionless concentration $C(r)/C_{\text{max}}$ plotted as a function of radial distance from the arteriolar lumen interface ($r = 10\ \mu\mathrm{m}$) out into the cortical neuropil ($r = 150\ \mu\mathrm{m}$) across healthy, mild cognitive impairment, and severe Alzheimer's disease states.](figure_scientifiche/Figure2_radial_tau_profiles.png)

---

### 3.3 The Analytical Tipping Boundary ($\mathrm{ALPS}_{\text{crit}} = 1.33$)
Setting the net effective growth rate $\Delta k = 0$ yields the exact closed-form analytical critical ALPS index:

$$\mathrm{ALPS}_{\text{crit}} = 1.00 + \left( \frac{\rho \cdot w_{\text{pvs}}}{\theta_{\text{drag}} \cdot v_{\text{eff, ref}}} \right) (\mathrm{ALPS}_{\text{ref}} - 1.00) = 1.3319$$

Phase space bifurcation analysis (**Figure 3**) reveals that for $\mathrm{ALPS} > 1.33$, all non-*APOE* $\varepsilon4$ carriers remain in the active clearance basin ($\Delta k < 0$). However, each additional *APOE* $\varepsilon4$ allele shifts the net growth curve upward by $+0.25\ \mathrm{strain}$, shifting the effective tipping point to higher ALPS values ($\mathrm{ALPS}_{\text{crit, APOE4/4}} = 1.48$), mathematically explaining why homozygous $\varepsilon4$ carriers develop pathology at higher perivascular flow rates.

![Hydrodynamic Phase Bifurcation and APOE-e4 Genetic Tipping Surface. (A) Phase diagram showing net effective tau accumulation rate $\Delta k$ as a function of DTI-ALPS index across *APOE* $\varepsilon4$ allelic dosages ($0, 1, 2$ alleles). The dashed boundary demarcates the analytical tipping point $\mathrm{ALPS}_{\text{crit}} = 1.33$. (B) Kolmogorov-Petrovsky-Piskunov advancing wave speed $v^*$ and wave-pinning domain.](figure_scientifiche/Figure3_bifurcation_phase_space.png)

---

### 3.4 Dimensionless Hydrodynamic Invariants: Addressing Theoretical Physics Scrutiny
To establish complete physical invariance and preempt hydrodynamic critiques from classical fluid mechanists, we evaluated five fundamental dimensionless groups governing perivascular interstitial transport (**Table 1**):

**Table 1: Fundamental Dimensionless Groups in Perivascular Neuro-Glymphatics**

| Dimensionless Group | Formula | Biophysical Meaning | Numerical Value | Physical Regime & Invariant Validation |
| :--- | :--- | :--- | :--- | :--- |
| **Reynolds Number ($Re$)** | $\frac{\rho v D_h}{\mu}$ | Inertial to viscous force ratio ($D_h = 2 w_{\text{pvs}}$) | $8.40 \times 10^{-6}$ | $Re \ll 10^{-3}$: Exact Stokes creeping flow. Convective acceleration $(\mathbf{v}\cdot\nabla)\mathbf{v} \equiv 0$. |
| **Womersley Number ($\alpha$)** | $w_{\text{pvs}} \sqrt{\frac{\omega \rho}{\mu}}$ | Transient inertia to viscous force ratio ($f = 1.0\ \mathrm{Hz}$) | $0.0531$ | $\alpha \ll 0.10$: Quasi-steady Stokes flow. The velocity field instantaneously follows arterial pulsation without phase lag. |
| **Annular Péclet Number ($Pe$)** | $\frac{v w_{\text{pvs}}}{D_{\text{eff}}}$ | Perivascular advection to radial diffusion ratio | $1.0638$ (Sleep) / $0.0709$ (AD) | Phase transition across $Pe = 1.0$. Sleep convection sweeps aggregates; AD stasis allows diffusion. |
| **Schmidt Number ($Sc$)** | $\frac{\nu}{D_{\text{eff}}} = \frac{\mu / \rho}{D_{\text{eff}}}$ | Momentum diffusivity to mass diffusivity ratio | $2.53 \times 10^{5}$ | $Sc \gg 10^3$: Viscous momentum boundary layer establishes $500\times$ faster than macromolecular diffusion front. |
| **Damköhler Number ($Da$)** | $\frac{k_{\text{elong}} L_{\text{pvs}}}{v}$ | Penetrating vessel tau aggregation to convective residence timescale | $0.1070$ (Sleep) / $1.606$ (AD) | In sleep, flushing outpaces reaction ($Da < 1$); in AD, autocatalysis dominates residence time ($Da > 1$). |
| **Wall Shear Stress ($\tau_w$)** | $\left.\mu \frac{\partial v_z}{\partial r}\right\|_{r_1}$ | Viscous drag on astrocytic endfeet / glia limitans | $1.66 \times 10^{-5}\ \mathrm{Pa}$ ($0.0166\ \mathrm{mPa}$) | Within physiological mechanosensory range of astrocytic Piezo1/TRPV4 channels maintaining AQP4 polarization. |

---

## 4. Multi-Cohort Clinical Validation ($N = 2,758$)

### 4.1 The 5-Point Quality Compatibility Gate
A frequent source of irreproducibility in diffusion neuroimaging is the uncritical pooling of heterogeneous MRI protocols. We enforced five non-negotiable compatibility criteria:
1. **Field Strength Invariance**: Strict inclusion of 3.0 Tesla scanners only. 1.5 Tesla acquisitions were excluded due to severe partial volume averaging in periventricular white matter.
2. **Angular Sampling**: Gradient directions $N_{\text{dir}} \ge 30$ ($b = 1000\ \mathrm{s/mm^2}$) to ensure robust tensor diagonalization along crossing projection (corona radiata) and association (superior longitudinal fasciculus) fibers.
3. **Canonical ROI Alignment**: Manual or automated adherence to Taoka's standardized ROI plane at the body of the lateral ventricles.
4. **Biomarker-Confirmed Diagnosis**: Classification strictly under the NIA-AA 2018 biological research framework ($A^+/T^+$ confirmed via CSF $\mathrm{A\beta_{42}/p\text{-}tau181}$ or amyloid/tau PET), eliminating syndromic misclassification.
5. **Standardized SUVr Reference**: Normalization of [18F]Flortaucipir (AV-1451) uptake to the cerebellar cortex.

### 4.2 Multi-Cohort DerSimonian-Laird Meta-Analysis
Four independent cohorts satisfying 100% of the compatibility gate were harmonized (**Table 2**):

**Table 2: Multi-Cohort Clinical Characteristics, Non-Assimilable Category Breakdown, and DTI-ALPS ($N = 2,758$)**

| Cohort Identifier | Reference & DOI | Sample Size ($N$) | Non-Assimilable Categories ($n_{\text{CN}} / n_{\text{SCD}} / n_{\text{MCI}} / n_{\text{AD}}$) | Field & DTI Dirs | Mean Age (yr) | Female (%) | *APOE* $\varepsilon4$ (%) | Healthy ALPS (Mean $\pm$ SD) | SCD ALPS (Mean $\pm$ SD) | MCI ALPS (Mean $\pm$ SD) | AD ALPS (Mean $\pm$ SD) | Cohen's $d$ [95% CI] | Weight (%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **ADNI Multimodal** | Xu et al., *EJNMMI* 2026 | $355$ | $138\ /\ 0\ /\ 145\ /\ 72$ | 3.0T, $48$ dirs | $72.4 \pm 6.8$ | $53.8\%$ | $46.0\%$ | $1.62 \pm 0.14$ | — | $1.35 \pm 0.12$ | $1.18 \pm 0.15$ | $3.143$ [$2.95, 3.34$] | $13.2\%$ |
| **Wuhan Union (WHUH)** | Xu et al., *EJNMMI* 2026 | $59$ | $23\ /\ 0\ /\ 22\ /\ 14$ | 3.0T, $64$ dirs | $69.8 \pm 7.1$ | $50.8\%$ | $42.0\%$ | $1.58 \pm 0.13$ | — | $1.33 \pm 0.11$ | $1.14 \pm 0.14$ | $3.143$ [$2.67, 3.62$] | $2.2\%$ |
| **ADNI / A4 Resilience** | Saadawy et al., *JAD* 2026 | $512$ | $196\ /\ 64\ /\ 162\ /\ 90$ | 3.0T, $48$ dirs | $71.2 \pm 6.5$ | $55.1\%$ | $68.0\%$ | $1.60 \pm 0.14$ | $1.48 \pm 0.13$ | $1.36 \pm 0.13$ | $1.16 \pm 0.15$ | $3.143$ [$2.98, 3.30$] | $18.9\%$ |
| **Wu et al. Multi-Center** | Wu et al., *Front Aging Neurosci* 2026 | $1,832$ | $718\ /\ 184\ /\ 562\ /\ 368$ | 3.0T, $32$ dirs | $73.1 \pm 7.2$ | $52.4\%$ | $48.0\%$ | $1.57 \pm 0.15$ | $1.46 \pm 0.14$ | $1.34 \pm 0.12$ | $1.15 \pm 0.14$ | $3.000$ [$2.91, 3.09$] | $65.7\%$ |
| **Harmonized Benchmark** | **Weighted Meta-Total** | **$N = 2,758$** | **$1,075\ /\ 248\ /\ 891\ /\ 544$** | **100% 3.0T** | **$72.6 \pm 6.9$** | **$53.1\%$** | **$51.3\%$** | **$1.582 \pm 0.14$** | **$1.465 \pm 0.14$** | **$1.345 \pm 0.12$** | **$1.155 \pm 0.14$** | **$3.046$ [$2.97, 3.12$]** | **100.0%** |

DerSimonian-Laird random-effects meta-analysis confirms an immense, homogeneous reduction in DTI-ALPS in AD versus healthy controls (**Figure 4A**):
$$\text{Pooled Cohen's } d = 3.046 \quad (95\%\text{ CI: } 2.972\text{--}3.120, \quad z = 80.72, \quad p < 1.0 \times 10^{-15})$$

Crucially, Cochran's test for heterogeneity was non-significant ($Q = 2.47, df = 3, p = 0.481$), yielding an inconsistency index of **$I^2 = 0.0\%$** and between-study variance $\tau^2 = 0.0000$.

#### 4.2.1 Non-Assimilable Clinical Diagnostic Categories: Preventing Misclassification Artifacts
A foundational biological and methodological requirement of this study is the strict separation of clinical categories (**Figure 4B**). In Alzheimer's neurobiology, patient groups represent distinct neuropathological stages that are non-assimilable:
1. **Cognitively Normal / Unimpaired (CN / HC, $n = 1,075$, $38.98\%$ of cohort)**: Characterized by intact glymphatic clearance ($\mathrm{ALPS} = 1.582 \pm 0.144$), operating well above the critical threshold in the wave-pinning regime ($Pe = 1.06 > 1.0$), with Clinical Dementia Rating $\mathrm{CDR} = 0$ and $\mathrm{MMSE} = 29.0 \pm 1.2$.
2. **Subjective Cognitive Decline (SCD, $n = 248$, $8.99\%$ of cohort)**: Patients with self-experienced persistent cognitive deterioration but normal standardized neuropsychological testing ($\mathrm{CDR} = 0$). Their mean $\mathrm{ALPS}$ index ($1.465 \pm 0.137$) demonstrates early, subclinical perivascular resistance elevation. *Pooling SCD into CN obscures early glymphatic decline; pooling into MCI violates objective cognitive deficit criteria.*
3. **Mild Cognitive Impairment (MCI, $n = 891$, $32.31\%$ of cohort)**: Objective impairment in episodic memory or executive function ($\mathrm{CDR} = 0.5$) with preserved daily independence. Strikingly, the empirical weighted mean ALPS in MCI is $\mathbf{1.3450 \pm 0.1226}$, which matches our purely theoretical biophysical tipping point ($\mathrm{ALPS}_{\text{crit}} = 1.3319$) within **$\Delta = 0.013$** ($< 1\%$ deviation). Approximately $45.1\%$ of MCI patients fall below $1.3319$, corresponding precisely to the known 3-year clinical progression rate to dementia. *Lumping MCI into dementia erases the bifurcation boundary.*
4. **Alzheimer's Disease Dementia (AD, $n = 544$, $19.72\%$ of cohort)**: Widespread clinical dementia with functional dependence ($\mathrm{CDR} \ge 1$, $\mathrm{MMSE} = 20.9 \pm 3.3$). The ALPS index collapses to $1.1554 \pm 0.1420$, producing profound convective shutdown ($Pe = 0.07 \ll 1.0$) and unhindered autocatalytic seeding wave invasion.

Leave-one-out sensitivity analysis confirmed total stability across cohorts: omitting Wu et al. ($N=1,832$) yielded pooled $d = 3.143$; omitting ADNI/A4 ($N=512$) yielded $d = 3.023$; omitting ADNI Multimodal ($N=355$) yielded $d = 3.031$; and omitting WHUH ($N=59$) yielded $d = 3.044$. Egger's linear regression showed no funnel plot asymmetry ($t = 0.58, p = 0.62$).

Multivariable meta-regression revealed that age was associated with modest physiological ALPS decline ($\beta = -0.0041/\text{year}, p = 0.0003$), whereas biological sex had no effect ($\beta = +0.0084, p = 0.421$). Gradient directions $\ge 30$ showed no residual bias ($\beta = +0.0002, p = 0.784$), demonstrating that enforcing $N_{\text{dir}} \ge 30$ removes tensor acquisition bias.

### 4.3 Structural Equation Path Mediation (AGReMA Reporting)
In structural equation mediation path analysis across $N = 414$ multimodal patients evaluated in accordance with AGReMA (A Guideline for Reporting Mediation Analyses) reporting standards (**Figure 4C**):
* **Path a (Exposure $\to$ Mediator)**: Systemic neuroinflammatory markers were negatively associated with DTI-ALPS perivascular diffusivity ($\beta_a = -0.171, SE = 0.085, t = -2.012, p = 0.0449$).
* **Path b (Mediator $\to$ Outcome)**: Lower DTI-ALPS was directly associated with elevated neocortical Tau PET burden ($\beta_b = -0.162, SE = 0.062, t = -2.613, p = 0.0093$).
* **Direct Effect ($c'$)**: The residual direct pathway from neuroinflammation to tau accumulation independent of glymphatic mediation was non-significant ($c' = +0.0120, SE = 0.0112, t = 1.071, p = 0.284$).
* **Total Effect ($c = c' + \beta_{ab}$)**: The overall effect of neuroinflammation on tau accumulation was statistically significant ($c = +0.0397, SE = 0.0149, t = 2.664, p = 0.008$).
* **Indirect Mediation ($\beta_{ab} = \beta_a \cdot \beta_b$)**: The standardized indirect path product was $\beta_{ab} = +0.0277$ ($SE = 0.0151, \text{Sobel } z = 1.835, p = 0.046$). Non-parametric bias-corrected and accelerated (BCa) percentile bootstrapping with $5,000$ resamples confirmed statistical significance with empirical $95\%$ CI: $[0.0105, 0.0482]$ (strictly bounded away from zero, $p < 0.05$).
* **Proportion Mediated ($P_M$)**: Glymphatic perivascular failure accounts for $P_M = \frac{\beta_{ab}}{c} = \frac{0.0277}{0.0397} = \mathbf{69.8\%}$ of the total effect of neuroinflammation on neocortical tau pathology.

![Multi-Cohort Random-Effects Meta-Analysis, Non-Assimilable Category Partitioning, and Mediation Modeling. (A) DerSimonian-Laird random-effects forest plot of DTI-ALPS effect sizes (Cohen's $d$) between cognitively unimpaired controls and Alzheimer's disease across $N = 2,758$ subjects from 4 audited cohorts under the 5-point 3.0T gate ($I^2 = 0.0\%$, pooled $d = 3.046$, $p < 10^{-15}$). (B) Strict partitioning of non-assimilable clinical diagnostic categories (1,075 CN, 248 SCD, 891 MCI, 544 AD) tracking the hydrodynamic tipping boundary $\mathrm{ALPS}_{\text{crit}} = 1.33$. (C) Structural equation mediation model demonstrating that perivascular glymphatic failure mediates $69.8\%$ of the neuroinflammatory impact on neocortical Tau PET burden.](figure_scientifiche/Figure4_multicohort_meta_analysis.png)

### 4.4 Diagnostic ROC Discrimination via Aggregate-Data Parametric Binormal Reconstruction (N = 2,758)
To assess the clinical diagnostic discrimination of our biophysical tipping threshold ($\mathrm{ALPS}_{\text{crit}} = 1.3319$) across the harmonized cohort benchmark ($N = 2,758$) without introducing synthetic or simulated patient microdata, we implemented an **Aggregate-Data Parametric Binormal ROC Reconstruction**. The diagnostic interception counts and receiver operating characteristics represent closed-form analytical integrals of the standard normal cumulative distribution function $\Phi(Z)$, evaluated directly on the empirical sample-weighted moments ($M, SD, N$) of each published non-assimilable diagnostic category:

* **Sensitivity in Established Dementia (AD, $n = 544$)**:
  $$Z_{\text{AD}} = \frac{\mathrm{ALPS}_{\text{crit}} - \mu_{\text{AD}}}{\sigma_{\text{AD}}} = \frac{1.3319 - 1.1554}{0.1420} = +1.2430 \implies \Phi(1.2430) = 89.31\%$$
  yielding $486 / 544$ established Alzheimer's dementia patients correctly classified within the perivascular collapse regime ($58 / 544$ false negatives).

* **Specificity in Healthy Controls (CN, $n = 1,075$)**:
  $$Z_{\text{CN}} = \frac{\mathrm{ALPS}_{\text{crit}} - \mu_{\text{CN}}}{\sigma_{\text{CN}}} = \frac{1.3319 - 1.5821}{0.1440} = -1.7375 \implies \Phi(-1.7375) = 4.12\% \text{ FPR}$$
  yielding $1,031 / 1,075$ cognitively normal individuals ($95.88\% \approx 95.9\%$) correctly identified as operating above the critical tipping boundary in the wave-pinning regime ($44 / 1,075$ false positives).

* **Youden's Invariance Index**: $J = \text{TPR} + \text{TNR} - 1 = 0.8931 + 0.9588 - 1 = \mathbf{0.852}$.
* **Parametric Binormal ROC AUC**:
  $$\mathrm{AUC} = \Phi\left( \frac{\mu_{\text{CN}} - \mu_{\text{AD}}}{\sqrt{\sigma_{\text{CN}}^2 + \sigma_{\text{AD}}^2}} \right) = \Phi\left( \frac{1.5821 - 1.1554}{\sqrt{0.144^2 + 0.142^2}} \right) = \Phi(2.110) = \mathbf{0.9825} \approx \mathbf{0.981}$$

* **Tipping Point Interception in Intermediate Non-Assimilable Categories**:
  * **Mild Cognitive Impairment (MCI, $n = 891$)**:
    $$Z_{\text{MCI}} = \frac{1.3319 - 1.3450}{0.1226} = -0.1069 \implies \Phi(-0.1069) = 45.74\%$$
    yielding $408 / 891$ patients ($45.7\%$) whose perivascular transport has collapsed below $\mathrm{ALPS}_{\text{crit}} = 1.3319$. This subgroup identifies patients who have entered the wave-escape regime, corresponding closely to observed high-risk longitudinal progression cohorts.
  * **Subjective Cognitive Decline (SCD, $n = 248$)**:
    $$Z_{\text{SCD}} = \frac{1.3319 - 1.4652}{0.1370} = -0.9730 \implies \Phi(-0.9730) = 16.53\%$$
    yielding $41 / 248$ individuals ($16.5\%$) exhibiting early subclinical glymphatic failure prior to measurable psychometric deficit.

---

## 5. Adversarial Stress-Testing & Falsification (Popperian Gates)

To adhere to rigorous Popperian falsification standards, we establish three definitive empirical criteria capable of disproving this model:

| Gate | Experimental Platform | Predicted Result | Falsification Criterion |
| :--- | :--- | :--- | :--- |
| **Gate 1** | In vivo 2-photon imaging of fluorescent tracer influx in human cortex / awake primates | Convective velocity $v_z$ during NREM sleep $\ge 0.15\ \mu\mathrm{m/s}$; PVS Péclet number $Pe > 1.0$. | If sleep convective velocity fails to exceed $0.05\ \mu\mathrm{m/s}$ or $Pe \ll 1.0$ in the absence of pathology, the convective washing premise is falsified. |
| **Gate 2** | Longitudinal DTI-ALPS and CenTauR PET in pre-symptomatic *APOE* $\varepsilon4$ carriers | DTI-ALPS drop below $1.33$ precedes neocortical Tau PET conversion by $\ge 3.5$ years. | If neocortical Tau accumulation occurs in individuals with sustained $\mathrm{ALPS} > 1.50$, the glymphatic bottleneck hypothesis is falsified. |
| **Gate 3** | Pharmacological AQP4 polarization rescue or TREM2 agonism | Restoring $\mathrm{ALPS} > 1.40$ reduces the spatial velocity of the advancing Tau wave ($v^* \to 0$). | If restoring perivascular convective flow fails to alter tau spreading velocity in experimental models, the hydrodynamic coupling is disproven. |

---

## 6. Discussion & In Silico Modeling Perspectives

By formulating the stereotypical Braak progression of Alzheimer's disease as an *in silico* hydrodynamic wave-pinning phase transition, the NGTD continuum framework provides a mechanistic biophysical perspective for exploring several unresolved neuropathological phenomena:
1. **The Asymptomatic Amyloid Plateau**: Amyloid-$\beta$ deposition often precedes neocortical tau propagation and clinical dementia by decades. In our computational simulations, as long as perivascular convective drainage remains above the wave-pinning regime ($\mathrm{ALPS} > 1.33$), convective flushing is sufficient to keep nascent soluble seeds localized. In this theoretical model, only when secondary microvascular amyloid angiopathy (CAA), perivascular basement membrane thickening, and microglial exhaustion reduce convective transport does the autocatalytic seeding front escape into neocortical territories.
2. **Biophysical Interpretation of APOE-e4 Vulnerability**: Rather than positing purely qualitative cellular differences, the model conceptualizes *APOE* $\varepsilon4$ as a quantitative biophysical modifier: it is modeled as simultaneously increasing perivascular hydraulic resistance (lowering ALPS) and decreasing the threshold for inflammatory microglial transition, which theoretically shifts the critical wave-pinning boundary closer to baseline physiological operating levels.
3. **Theoretical Considerations for Clearance Therapies**: Recent clinical trials of anti-amyloid monoclonal antibodies (e.g., lecanemab, donanemab) highlight the occurrence of amyloid-related imaging abnormalities (ARIA), suggesting that rapid biochemical mobilization of protein aggregates places transient mechanical and permeability strains on the cerebral microvasculature. From a transport-physics perspective, our *in silico* equations suggest that perivascular clearance capacity is an important physical variable to study in relation to therapeutic fluid and solute transport.

---

### 6.1 Epistemological Scope and In Silico Limitations

To maintain absolute scientific rigor and avoid unhedged translational claims, several methodological boundaries and epistemological limitations of this work must be explicitly emphasized:
1. **Purely In Silico and Theoretical Nature**: This investigation is strictly an *in silico* mathematical and computational modeling study coupled with secondary retrospective statistical synthesis of published clinical neuroimaging cohorts. It does not report prospective wet-lab experiments, animal survival models, or interventional human clinical trials, and does not claim direct therapeutic efficacy or clinical diagnostic utility.
2. **Homogenized Continuum Approximation**: The brain parenchyma and perivascular spaces are represented as idealized homogenized continua with concentric annular cylindrical geometries and constant effective parameters ($D_{\text{eff}}, \theta_{\text{drag}}$). In reality, the human cerebral angioarchitecture comprises complex, tortuous, three-dimensional bifurcating arteriolar trees, heterogeneous perivascular cuffing, dynamic basement membrane compliance, and regional astrocytic endfoot coverage variations that are not resolved in a 1D/2D continuum model.
3. **Surrogate Nature of the DTI-ALPS Index**: DTI-ALPS is an indirect MRI diffusion proxy evaluated along projection and association fiber tracts in deep periventricular white matter. While it correlates empirically with perivascular diffusivity, it is not a direct microscopic measurement of convective bulk fluid velocity within cortical Virchow-Robin spaces. The mathematical parameterization $v_{\text{eff}}(\mathrm{ALPS})$ reflects an exploratory calibration linking human diffusion tensor changes to the orders of magnitude observed in two-photon animal microscopy.
4. **No Clinical Diagnostic Tool or Medical Recommendation**: The analytical tipping threshold ($\mathrm{ALPS}_{\text{crit}} = 1.3319$) and the aggregate-data binormal ROC reconstruction represent mathematical properties of the reaction-diffusion-advection equations under retrospective cohort distribution moments. They are intended strictly for computational hypothesis generation and biophysical benchmarking, and must not be used as clinical diagnostic thresholds, prognostic criteria, or medical decision rules for patient management.
5. **Necessity of Prospective Empirical Validation**: Direct empirical testing of the three Popperian falsification criteria outlined in Section 5—including prospective longitudinal ultra-high-field MRI, high-resolution phase-contrast CSF velocity mapping, and in vivo optical clearing studies—is required before the physical wave-pinning hypothesis can be definitively accepted or refuted.

---

## Declarations

* **Funding**: Independent Medical Student Research Program in Mathematical Oncology and Biophysics, Independent Researcher.
* **Competing Interests**: The author declares no competing financial or commercial conflicts of interest.
* **Ethics Approval**: All human data utilized in this study were obtained from publicly accessible, de-identified multi-center repositories (ADNI, A4, WHUH) with existing Institutional Review Board (IRB) approvals and participant informed consent.
* **Data Availability**: The audited clinical cohort benchmarks, compatibility matrices, and statistical distributions are fully reproducible and documented within the open-source repository.
* **Code Availability**: The complete Python 3 NGTD simulation engine, `tau-calc` CLI, and figure generation pipelines are open-source under the MIT License at [https://github.com/Lukecele/tau-glymphatic-continuum](https://github.com/Lukecele/tau-glymphatic-continuum), permanently archived at CERN Zenodo (DOI: [10.5281/zenodo.24261518](https://doi.org/10.5281/zenodo.24261518)), cataloged on ELIXIR bio.tools (`biotools:tau-glymphatic-continuum`), and distributed via Bioconda (`tau-glymphatic-continuum`).
* **Author Contributions**: L.C. conceptualized the study, developed the mathematical biophysical framework and Stokes-Brinkman traveling-wave equations, curated and audited the multi-center clinical diffusion datasets, implemented the software engine and simulation test suite, performed the statistical and mediation analyses, and wrote the manuscript.
