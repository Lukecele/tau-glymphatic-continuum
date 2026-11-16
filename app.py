"""
Neuro-Glymphatic Tau Dynamics (NGTD) — Interactive Bedside Clinical Simulator.

A multiscale biophysical reaction-diffusion and perivascular transport dashboard
for translational neurologists, biophysicists, and neuroscientists.

Author: Luca Celebrano, MD Candidate (Università degli Studi di Napoli Federico II / TIGEM)
Manuscript ID: BIORXIV/2026/766809
Certified DOI: 10.5281/zenodo.24261518
"""

import sys
import os
import math
from pathlib import Path

# Add src to python path for tau_glymphatic engine
current_dir = Path(__file__).resolve().parent
src_dir = current_dir / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

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
from tau_glymphatic.models import (
    DEFAULT_GEOMETRY,
    DEFAULT_PHYSICS,
    DEFAULT_GLYMPHATIC,
    DEFAULT_TAU_KINETICS,
    DEFAULT_MICROGLIA,
)

# Page configuration
st.set_page_config(
    page_title="NGTD — Neuro-Glymphatic Tau Simulator",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom Styling
st.markdown(
    """
    <style>
    .main-title {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E3A8A;
        margin-bottom: 0.2rem;
    }
    .sub-title {
        font-size: 1.1rem;
        color: #4B5563;
        margin-bottom: 1.2rem;
    }
    .card-safe {
        background-color: #ECFDF5;
        border-left: 5px solid #10B981;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
    }
    .card-danger {
        background-color: #FEF2F2;
        border-left: 5px solid #EF4444;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
    }
    .citation-box {
        background-color: #F3F4F6;
        padding: 0.75rem;
        border-radius: 0.375rem;
        font-size: 0.85rem;
        color: #374151;
        margin-top: 1rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# Title & Metadata
st.markdown('<div class="main-title">🧠 Neuro-Glymphatic Tau Dynamics (NGTD)</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="sub-title">Interactive Multiscale Reaction-Diffusion & Perivascular Clearance Engine | Bedside Risk Simulator</div>',
    unsafe_allow_html=True,
)

# Sidebar: Patient Biomarkers
st.sidebar.header("🔬 Patient Clinical Profile")

# Presets
preset = st.sidebar.selectbox(
    "Clinical Cohort Preset:",
    ["Custom Parameters", "Healthy Young Adult", "Healthy Elderly (Control)", "Prodromal MCI", "Clinical Alzheimer's Disease"],
)

if preset == "Healthy Young Adult":
    default_alps = 1.75
    default_apoe = 0
    default_sleep = 8.0
    default_trem2 = False
    default_nlrp3 = False
elif preset == "Healthy Elderly (Control)":
    default_alps = 1.55
    default_apoe = 0
    default_sleep = 7.0
    default_trem2 = False
    default_nlrp3 = False
elif preset == "Prodromal MCI":
    default_alps = 1.34
    default_apoe = 1
    default_sleep = 6.0
    default_trem2 = False
    default_nlrp3 = True
elif preset == "Clinical Alzheimer's Disease":
    default_alps = 1.15
    default_apoe = 2
    default_sleep = 5.0
    default_trem2 = True
    default_nlrp3 = True
else:
    default_alps = 1.45
    default_apoe = 0
    default_sleep = 7.5
    default_trem2 = False
    default_nlrp3 = False

st.sidebar.subheader("1. Diffusion Tensor Imaging (DTI)")
alps_index = st.sidebar.slider(
    "DTI-ALPS Index (Along Perivascular Space):",
    min_value=0.80,
    max_value=2.00,
    value=float(default_alps),
    step=0.01,
    help="Ratio of diffusivity along projection/association fibers at the level of the lateral ventricle body. Healthy > 1.45; AD < 1.25.",
)

sleep_hours = st.sidebar.slider(
    "Slow-Wave Sleep Duration (hours/night):",
    min_value=3.0,
    max_value=10.0,
    value=float(default_sleep),
    step=0.5,
    help="Perivascular interstitial fluid (ISF) bulk flow accelerates dramatically during slow-wave N3 NREM sleep.",
)

st.sidebar.subheader("2. Genetic & Inflammatory Modifiers")
apoe_e4_count = st.sidebar.radio(
    "APOE-ε4 Allele Burden:",
    options=[0, 1, 2],
    index=int(default_apoe),
    format_func=lambda x: f"{x} Allele(s)" + (" (E3/E3 neutral)" if x == 0 else " (Heterozygous E3/E4)" if x == 1 else " (Homozygous E4/E4)"),
)

trem2_r47h = st.sidebar.checkbox(
    "TREM2 R47H Loss-of-Function Variant",
    value=default_trem2,
    help="Impairs microglial phagocytic clearance of oligomeric tau seeds.",
)

nlrp3_activated = st.sidebar.checkbox(
    "NLRP3 Inflammasome Hyperactivation (Neuroinflammation)",
    value=default_nlrp3,
    help="Elevated systemic/central neuroinflammation suppresses baseline microglial degradation efficiency.",
)

# Core Computations
sleep_fraction = sleep_hours / 24.0
v_eff = calculate_effective_velocity(alps_index, sleep_fraction=sleep_fraction)
v_eff_nm_s = v_eff * 1e9

peclet_peak = calculate_peclet_number(alps_index, peak_sleep=True)
peclet_mean = calculate_peclet_number(alps_index, peak_sleep=False)

crit_info = solve_alps_critical_threshold()
alps_crit = crit_info["alps_crit_kinetic"]

rates = calculate_net_clearance_rate(
    alps_index=alps_index,
    apoe_e4_alleles=apoe_e4_count,
    trem2_r47h=trem2_r47h,
    nlrp3_activated=nlrp3_activated,
)

gamma_eff_per_year = rates["effective_growth_rate"] * (365.25 * 86400)
is_pathologic = rates["is_pathologic_expansion"]

# Top KPI Metric Cards
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        label="DTI-ALPS Index",
        value=f"{alps_index:.2f}",
        delta=f"{alps_index - alps_crit:+.2f} vs ALPS_crit" if alps_index != alps_crit else "At critical point",
        delta_color="normal" if not is_pathologic else "inverse",
    )

with col2:
    st.metric(
        label="Critical Collapse Threshold",
        value=f"{alps_crit:.4f}",
        help="Analytical bifurcation threshold: below this value, convective clearance fails to balance auto-catalytic tau aggregation.",
    )

with col3:
    st.metric(
        label="Annular Péclet Number (Pe)",
        value=f"{peclet_peak:.3f}",
        delta="Convective regime (Pe > 1)" if peclet_peak >= 1.0 else "Diffusive stasis (Pe < 1)",
        delta_color="normal" if peclet_peak >= 1.0 else "inverse",
    )

with col4:
    st.metric(
        label="Net Tau Expansion Rate",
        value=f"{gamma_eff_per_year:+.2f} yr⁻¹",
        delta="Clearance Dominant" if not is_pathologic else "Seeding Accumulation",
        delta_color="normal" if not is_pathologic else "inverse",
    )

# Clinical Status Banner
if not is_pathologic:
    st.markdown(
        f"""
        <div class="card-safe">
            <strong>🟢 CLINICAL CLEARANCE STABILITY: SAFE CONVECTIVE REGIME</strong><br>
            Patient ALPS index (<b>{alps_index:.2f}</b>) exceeds the critical bifurcation threshold (<b>{alps_crit:.4f}</b>). 
            Bulk perivascular convective flux (effective velocity <b>{v_eff_nm_s:.1f} nm/s</b>) successfully clears monomeric and oligomeric tau species before 
            surface-catalyzed secondary nucleation takes over.
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        f"""
        <div class="card-danger">
            <strong>🔴 HYDRODYNAMIC CLEARANCE COLLAPSE: PATHOLOGICAL SEEDING REGIME</strong><br>
            Patient ALPS index (<b>{alps_index:.2f}</b>) has breached the critical threshold (<b>{alps_crit:.4f}</b>).
            Hydrodynamic convective failure causes interstitial solute stagnation. The rate of Knowles-Vendruscolo autocatalytic tau nucleation 
            exceeds total clearance by <b>{abs(gamma_eff_per_year):.2f} yr⁻¹</b>, predicting exponential neocortical Braak stage propagation.
        </div>
        """,
        unsafe_allow_html=True,
    )

# Visualizations
tab1, tab2, tab3, tab4 = st.tabs([
    "📈 Radial Tau Spreading Profile",
    "⚡ Bifurcation Phase Diagram",
    "⏳ 10-Year Tau-PET Trajectory",
    "📐 Mathematical Foundations & Equations",
])

with tab1:
    st.subheader("Steady-State Radial Tau Concentration: Perivascular Annulus to Neuropil")
    st.caption("Coupled solution of Stokes-Brinkman advection and Fisher-KPP reaction-diffusion with microglial inflammatory strain.")

    profile_data = simulate_radial_pvs_profile(
        alps_index=alps_index,
        apoe_e4_alleles=apoe_e4_count,
        trem2_r47h=trem2_r47h,
        nlrp3_activated=nlrp3_activated,
        num_grid_points=80,
    )

    r_um = np.array(profile_data["r_points_um"])
    c_profile = np.array(profile_data["c_profile"])

    # Comparison baseline (Healthy young adult ALPS = 1.70, no risk alleles)
    healthy_data = simulate_radial_pvs_profile(alps_index=1.70, apoe_e4_alleles=0, num_grid_points=80)
    c_healthy = np.array(healthy_data["c_profile"])

    # Severe AD baseline (ALPS = 1.10, APOE e4/e4)
    ad_data = simulate_radial_pvs_profile(alps_index=1.10, apoe_e4_alleles=2, trem2_r47h=True, nlrp3_activated=True, num_grid_points=80)
    c_ad = np.array(ad_data["c_profile"])

    fig, ax = plt.subplots(figsize=(10, 4.5), dpi=150)
    ax.plot(r_um, c_healthy, color="#10B981", linestyle="--", linewidth=1.8, label="Healthy Control Reference (ALPS 1.70, E3/E3)")
    ax.plot(r_um, c_ad, color="#EF4444", linestyle=":", linewidth=1.8, label="Severe AD Baseline (ALPS 1.10, E4/E4)")
    ax.plot(r_um, c_profile, color="#1E3A8A", linewidth=2.8, label=f"Current Patient Simulation (ALPS {alps_index:.2f}, {apoe_e4_count} E4)")

    ax.axvspan(r_um[0], r_um[0] + DEFAULT_GEOMETRY.pvs_width * 1e6, color="#E0F2FE", alpha=0.4, label="Perivascular Space (PVS)")
    ax.axvspan(r_um[0] + DEFAULT_GEOMETRY.pvs_width * 1e6, r_um[-1], color="#F3F4F6", alpha=0.4, label="Brain Parenchyma (Neuropil)")

    ax.set_xlabel("Radial Distance from Vessel Center r (μm)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Normalized Soluble Tau Concentration C(r) / C_max", fontsize=11, fontweight="bold")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    plt.tight_layout()
    st.pyplot(fig)

    st.markdown(
        f"**Hydrodynamic Penetration Length**: $\\lambda_c = \\sqrt{{D_{{\\text{{eff}}}} / k_{{\\text{{clear}}}}}} = \\mathbf{{{profile_data['lambda_clearance_um']:.1f}\\,\\mu\\text{{m}}}}$. "
        f"Drainage boundary sink efficiency at the glia limitans: $\\mathbf{{{100*(1.0 - profile_data['c_vessel_boundary']):.1f}\\%}}$."
    )

with tab2:
    st.subheader("Hydrodynamic Bifurcation: Effective Growth Rate vs. DTI-ALPS Index")
    st.caption("Visualizing the critical transition point where convective clearance breaks down.")

    alps_range = np.linspace(0.90, 1.90, 120)
    growth_rates = []
    for a in alps_range:
        r = calculate_net_clearance_rate(
            alps_index=float(a),
            apoe_e4_alleles=apoe_e4_count,
            trem2_r47h=trem2_r47h,
            nlrp3_activated=nlrp3_activated,
        )
        growth_rates.append(r["effective_growth_rate"] * (365.25 * 86400))

    fig2, ax2 = plt.subplots(figsize=(10, 4.5), dpi=150)
    ax2.plot(alps_range, growth_rates, color="#3B82F6", linewidth=2.5, label="Net Growth Rate γ_eff(ALPS)")
    ax2.axhline(0.0, color="#6B7280", linestyle="-", linewidth=1.2)
    ax2.axvline(alps_crit, color="#DC2626", linestyle="--", linewidth=2.0, label=f"ALPS_crit = {alps_crit:.4f}")

    # Highlight patient point
    ax2.scatter(
        [alps_index],
        [gamma_eff_per_year],
        color="#DC2626" if is_pathologic else "#10B981",
        s=120,
        zorder=5,
        label=f"Patient Operating Point ({alps_index:.2f}, {gamma_eff_per_year:+.2f} yr⁻¹)",
    )

    # Shaded regions
    ax2.axvspan(0.90, alps_crit, color="#FEE2E2", alpha=0.5, label="Pathological Expansion Zone (γ_eff > 0)")
    ax2.axvspan(alps_crit, 1.90, color="#D1FAE5", alpha=0.5, label="Clearance Stability Zone (γ_eff ≤ 0)")

    ax2.set_xlabel("Clinical DTI-ALPS Index", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Net Expansion Rate γ_eff (year⁻¹)", fontsize=11, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", frameon=True, fontsize=9)
    plt.tight_layout()
    st.pyplot(fig2)

with tab3:
    st.subheader("10-Year Neocortical Tau-PET SUVR Accumulation Projection")
    st.caption("Longitudinal trajectory based on Fisher-KPP logistic saturation with current vs. optimized glymphatic clearance.")

    years = np.linspace(0, 10, 100)
    base_suvr = 1.10  # Baseline sub-threshold tau PET SUVR (entorhinal cortex)
    max_suvr = 2.80   # Braak Stage VI saturation level

    # Patient trajectory
    growth_param = max(-0.05, rates["effective_growth_rate"] * (365.25 * 86400))
    if growth_param > 0:
        trajectory_patient = max_suvr / (1.0 + ((max_suvr - base_suvr) / base_suvr) * np.exp(-growth_param * years * 0.35))
    else:
        trajectory_patient = np.maximum(1.02, base_suvr * np.exp(growth_param * years * 0.15))

    # Optimized lifestyle intervention scenario (Sleep optimized to 8.5h + ALPS improved by +0.25 via hypertension control)
    alps_opt = min(1.80, alps_index + 0.25)
    rates_opt = calculate_net_clearance_rate(
        alps_index=alps_opt,
        apoe_e4_alleles=apoe_e4_count,
        trem2_r47h=False,
        nlrp3_activated=False,
    )
    growth_opt = max(-0.08, rates_opt["effective_growth_rate"] * (365.25 * 86400))
    if growth_opt > 0:
        trajectory_opt = max_suvr / (1.0 + ((max_suvr - base_suvr) / base_suvr) * np.exp(-growth_opt * years * 0.35))
    else:
        trajectory_opt = np.maximum(1.02, base_suvr * np.exp(growth_opt * years * 0.15))

    fig3, ax3 = plt.subplots(figsize=(10, 4.5), dpi=150)
    ax3.plot(years, trajectory_patient, color="#EF4444" if is_pathologic else "#3B82F6", linewidth=2.5, label="Current Patient Trajectory")
    ax3.plot(years, trajectory_opt, color="#10B981", linestyle="--", linewidth=2.0, label=f"Optimized Protocol (Sleep 8.5h + ALPS +0.25)")
    ax3.axhline(1.35, color="#F59E0B", linestyle=":", linewidth=1.5, label="Tau-PET Positivity Threshold (SUVR = 1.35)")

    ax3.set_xlabel("Follow-up Duration (Years)", fontsize=11, fontweight="bold")
    ax3.set_ylabel("Neocortical Tau-PET SUVR", fontsize=11, fontweight="bold")
    ax3.set_ylim(0.95, 3.00)
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.legend(loc="upper left", frameon=True, fontsize=9)
    plt.tight_layout()
    st.pyplot(fig3)

with tab4:
    st.subheader("Mathematical Model & Exact Analytical Derivations")
    st.markdown(
        r"""
        ### 1. Coupled Stokes-Brinkman & Porous Annulus Mechanics
        Flow within the annular perivascular space (PVS, radius $r_1 \le r \le r_2$) is governed by momentum transport with glia limitans Beavers-Joseph-Saffman slip:
        $$\mu \nabla^2 \mathbf{u} - \nabla p - \frac{\mu}{K_p} \mathbf{u} = 0$$

        ### 2. Knowles-Vendruscolo Dual-Species Nucleated Polymerization
        The kinetic evolution of soluble monomer $m(t)$ and fibrillar polymer mass $M(t)$ obeys:
        $$\frac{\partial M}{\partial t} = 2 k_e m(t) P(t) + k_2 m(t) M(t)^2 - k_{\text{glymph}} M(t) - k_{\text{microglia}} M(t)$$

        ### 3. Closed-Form Critical Collapse Threshold ($\mathrm{ALPS}_{\text{crit}}$)
        At the bifurcation threshold $\gamma_{\text{eff}} = 0$, convective sink clearance exactly equals autocatalytic aggregation:
        $$k_{\text{glymph}}(\mathrm{ALPS}_{\text{crit}}) = \rho_{\text{agg}} = 2 k_e m_0 P_0 + k_2 m_0^2$$
        $$\mathrm{ALPS}_{\text{crit}} = \mathrm{ALPS}_{\text{floor}} + \left(\frac{\rho_{\text{agg}} L_{\text{pvs}}}{\theta \cdot v_{\text{ref}}}\right) \cdot (\mathrm{ALPS}_{\text{healthy}} - \mathrm{ALPS}_{\text{floor}}) = \mathbf{1.3319}$$
        """
    )

# Footer Citation & Links
st.markdown("---")
st.markdown(
    """
    <div class="citation-box">
        <strong>Citation Metadata:</strong><br>
        Celebrano, L. (2026). <em>Neuro-Glymphatic Tau Dynamics: An In Silico Multiscale Reaction Diffusion and Perivascular Transport Model of Pathological Tau Spreading in Alzheimer's Disease</em>. 
        bioRxiv preprint <code>BIORXIV/2026/766809</code>. Certified Archive DOI: <a href="https://doi.org/10.5281/zenodo.24261518" target="_blank">10.5281/zenodo.24261518</a>.
    </div>
    """,
    unsafe_allow_html=True,
)
