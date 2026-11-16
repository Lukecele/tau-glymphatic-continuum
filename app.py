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

# Initialize session state for reactive inputs
if "alps_val" not in st.session_state:
    st.session_state.alps_val = 1.45
if "sleep_val" not in st.session_state:
    st.session_state.sleep_val = 7.5
if "apoe_val" not in st.session_state:
    st.session_state.apoe_val = 0
if "trem2_val" not in st.session_state:
    st.session_state.trem2_val = False
if "nlrp3_val" not in st.session_state:
    st.session_state.nlrp3_val = False

def on_preset_change():
    sel = st.session_state.preset_select
    if sel == "Healthy Young Adult":
        st.session_state.alps_val = 1.75
        st.session_state.sleep_val = 8.0
        st.session_state.apoe_val = 0
        st.session_state.trem2_val = False
        st.session_state.nlrp3_val = False
    elif sel == "Healthy Elderly (Control)":
        st.session_state.alps_val = 1.55
        st.session_state.sleep_val = 7.0
        st.session_state.apoe_val = 0
        st.session_state.trem2_val = False
        st.session_state.nlrp3_val = False
    elif sel == "Prodromal MCI":
        st.session_state.alps_val = 1.28
        st.session_state.sleep_val = 6.0
        st.session_state.apoe_val = 1
        st.session_state.trem2_val = False
        st.session_state.nlrp3_val = True
    elif sel == "Clinical Alzheimer's Disease":
        st.session_state.alps_val = 1.12
        st.session_state.sleep_val = 5.0
        st.session_state.apoe_val = 2
        st.session_state.trem2_val = True
        st.session_state.nlrp3_val = True

preset = st.sidebar.selectbox(
    "Clinical Cohort Preset:",
    ["Custom Parameters", "Healthy Young Adult", "Healthy Elderly (Control)", "Prodromal MCI", "Clinical Alzheimer's Disease"],
    key="preset_select",
    on_change=on_preset_change,
)

st.sidebar.subheader("1. Diffusion Tensor Imaging (DTI)")
alps_index = st.sidebar.slider(
    "DTI-ALPS Index (Along Perivascular Space):",
    min_value=0.80,
    max_value=2.00,
    step=0.01,
    key="alps_val",
    help="Diffusivity ratio along projection/association fibers. Threshold ALPS_crit = 1.3319; Healthy > 1.45; AD < 1.25.",
)

sleep_hours = st.sidebar.slider(
    "Slow-Wave Sleep Duration (hours/night):",
    min_value=3.0,
    max_value=10.0,
    step=0.5,
    key="sleep_val",
    help="Perivascular interstitial fluid (ISF) convective flux accelerates dramatically during slow-wave N3 NREM sleep.",
)

st.sidebar.subheader("2. Genetic & Inflammatory Modifiers")
apoe_e4_count = st.sidebar.radio(
    "APOE-ε4 Allele Burden:",
    options=[0, 1, 2],
    format_func=lambda x: f"{x} Allele(s)" + (" (E3/E3 neutral)" if x == 0 else " (Heterozygous E3/E4)" if x == 1 else " (Homozygous E4/E4)"),
    key="apoe_val",
)

trem2_r47h = st.sidebar.checkbox(
    "TREM2 R47H Loss-of-Function Variant",
    key="trem2_val",
    help="Impairs microglial phagocytic clearance of oligomeric tau seeds.",
)

nlrp3_activated = st.sidebar.checkbox(
    "NLRP3 Inflammasome Hyperactivation (Neuroinflammation)",
    key="nlrp3_val",
    help="Elevated neuroinflammation triggers exosomal seed release and suppresses degradation.",
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
    sleep_fraction=sleep_fraction,
)

gamma_eff_per_year = rates["effective_growth_rate"] * (365.25 * 86400)
is_pathologic = rates["is_pathologic_expansion"]

# Top KPI Metric Cards
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        label="DTI-ALPS Index",
        value=f"{alps_index:.2f}",
        delta=f"{alps_index - alps_crit:+.2f} vs ALPS_crit" if alps_index != alps_crit else "At critical bifurcation",
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
        delta="Convective Washout (Pe > 1)" if peclet_peak >= 1.0 else "Diffusive Stasis (Pe < 1)",
        delta_color="normal" if peclet_peak >= 1.0 else "inverse",
    )

with col4:
    st.metric(
        label="Net Tau Expansion Rate",
        value=f"{gamma_eff_per_year:+.2f} yr⁻¹",
        delta="Clearance Dominant" if not is_pathologic else "Pathological Accumulation",
        delta_color="normal" if not is_pathologic else "inverse",
    )

# Clinical Status Banner
if not is_pathologic:
    st.markdown(
        f"""
        <div class="card-safe">
            <strong>🟢 CLINICAL CLEARANCE STABILITY: SAFE CONVECTIVE REGIME</strong><br>
            Patient ALPS index (<b>{alps_index:.2f}</b>) provides sufficient perivascular flux (effective velocity <b>{v_eff_nm_s:.1f} nm/s</b>) 
            to continuously clear monomeric and oligomeric tau species (Net growth: <b>{gamma_eff_per_year:+.2f} yr⁻¹</b>). 
            Protective slow-wave convective drainage suppresses secondary autocatalytic nucleation.
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        f"""
        <div class="card-danger">
            <strong>🔴 HYDRODYNAMIC CLEARANCE COLLAPSE: PATHOLOGICAL SEEDING REGIME</strong><br>
            Patient ALPS index (<b>{alps_index:.2f}</b>) is in hydrodynamic failure (effective velocity <b>{v_eff_nm_s:.1f} nm/s</b>).
            Convective clearance has collapsed below autocatalytic secondary nucleation (Net expansion: <b>{gamma_eff_per_year:+.2f} yr⁻¹</b>).
            Perivascular solute stagnation accelerates neurofibrillary tangle propagation across neocortical networks.
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
    st.caption("Real-time solution of perivascular advection-diffusion coupled to boundary sink drainage.")

    profile_data = simulate_radial_pvs_profile(
        alps_index=alps_index,
        apoe_e4_alleles=apoe_e4_count,
        trem2_r47h=trem2_r47h,
        nlrp3_activated=nlrp3_activated,
        sleep_fraction=sleep_fraction,
        num_grid_points=80,
    )

    r_um = np.array(profile_data["r_points_um"])
    c_profile = np.array(profile_data["c_profile"])

    # Reference baselines
    healthy_data = simulate_radial_pvs_profile(alps_index=1.70, apoe_e4_alleles=0, sleep_fraction=8.0/24.0, num_grid_points=80)
    c_healthy = np.array(healthy_data["c_profile"])

    ad_data = simulate_radial_pvs_profile(alps_index=1.10, apoe_e4_alleles=2, trem2_r47h=True, nlrp3_activated=True, sleep_fraction=5.0/24.0, num_grid_points=80)
    c_ad = np.array(ad_data["c_profile"])

    fig, ax = plt.subplots(figsize=(10, 4.5), dpi=150)
    ax.plot(r_um, c_healthy, color="#10B981", linestyle="--", linewidth=1.8, label="Healthy Reference (ALPS 1.70, E3/E3, 8h Sleep)")
    ax.plot(r_um, c_ad, color="#EF4444", linestyle=":", linewidth=1.8, label="Severe AD Reference (ALPS 1.10, E4/E4, 5h Sleep)")
    
    # Active patient curve with dynamic styling
    patient_color = "#DC2626" if is_pathologic else "#1E3A8A"
    ax.plot(r_um, c_profile, color=patient_color, linewidth=2.8, label=f"Current Patient (ALPS {alps_index:.2f}, {apoe_e4_count} E4, {sleep_hours:.1f}h Sleep)")
    ax.scatter([r_um[0]], [c_profile[0]], color=patient_color, s=80, zorder=6)

    # Shaded anatomic zones
    pvs_outer_um = r_um[0] + DEFAULT_GEOMETRY.pvs_width * 1e6
    ax.axvspan(r_um[0], pvs_outer_um, color="#E0F2FE", alpha=0.45, label="Perivascular Space (PVS Annulus)")
    ax.axvspan(pvs_outer_um, r_um[-1], color="#F3F4F6", alpha=0.35, label="Brain Parenchyma (Neuropil)")

    ax.set_xlabel("Radial Distance from Vessel Center r (μm)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Normalized Soluble Tau Concentration C(r) / C_max", fontsize=11, fontweight="bold")
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right", frameon=True, fontsize=9)
    plt.tight_layout()
    st.pyplot(fig)
    plt.close(fig)

    sink_efficiency = 100.0 * (1.0 - profile_data["c_vessel_boundary"])
    st.info(
        f"💡 **Real-Time Biophysical Telemetry**: Perivascular Wall Sink Efficiency = **{sink_efficiency:.1f}%** | "
        f"Boundary Concentration $C(r_{{\\text{{vessel}}}})$ = **{profile_data['c_vessel_boundary']:.3f}** | "
        f"Characteristic Drainage Boundary Layer $\\delta_{{\\text{{eff}}}}$ = **{profile_data['delta_penetration_um']:.1f}\\,\\mu\\text{{m}}**."
    )

with tab2:
    st.subheader("Hydrodynamic Bifurcation: Effective Growth Rate vs. DTI-ALPS Index")
    st.caption("Visualizing the exact bifurcation point where convective clearance breaks down.")

    alps_range = np.linspace(0.90, 1.90, 120)
    growth_rates = []
    for a in alps_range:
        r = calculate_net_clearance_rate(
            alps_index=float(a),
            apoe_e4_alleles=apoe_e4_count,
            trem2_r47h=trem2_r47h,
            nlrp3_activated=nlrp3_activated,
            sleep_fraction=sleep_fraction,
        )
        growth_rates.append(r["effective_growth_rate"] * (365.25 * 86400))

    fig2, ax2 = plt.subplots(figsize=(10, 4.5), dpi=150)
    ax2.plot(alps_range, growth_rates, color="#3B82F6", linewidth=2.5, label="Net Growth Rate γ_eff(ALPS)")
    ax2.axhline(0.0, color="#6B7280", linestyle="-", linewidth=1.2)
    ax2.axvline(alps_crit, color="#DC2626", linestyle="--", linewidth=1.8, label=f"ALPS_crit = {alps_crit:.4f}")

    # Dynamic crosshair for current patient operating point
    ax2.axvline(alps_index, color=patient_color, linestyle=":", alpha=0.6)
    ax2.axhline(gamma_eff_per_year, color=patient_color, linestyle=":", alpha=0.6)
    ax2.scatter(
        [alps_index],
        [gamma_eff_per_year],
        color=patient_color,
        s=140,
        zorder=5,
        edgecolor="black",
        label=f"Patient Operating Point ({alps_index:.2f}, {gamma_eff_per_year:+.2f} yr⁻¹)",
    )

    # Shaded stability vs collapse zones
    ax2.axvspan(0.90, alps_crit, color="#FEE2E2", alpha=0.45, label="Pathological Expansion (γ_eff > 0)")
    ax2.axvspan(alps_crit, 1.90, color="#D1FAE5", alpha=0.45, label="Clearance Stability (γ_eff ≤ 0)")

    ax2.set_xlabel("Clinical DTI-ALPS Index", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Net Expansion Rate γ_eff (year⁻¹)", fontsize=11, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", frameon=True, fontsize=9)
    plt.tight_layout()
    st.pyplot(fig2)
    plt.close(fig2)

with tab3:
    st.subheader("10-Year Neocortical Tau-PET SUVR Accumulation Projection")
    st.caption("Longitudinal trajectory based on Fisher-KPP logistic saturation with current vs. optimized glymphatic clearance.")

    years = np.linspace(0, 10, 100)
    base_suvr = 1.10  # Baseline sub-threshold tau PET SUVR (entorhinal cortex)
    max_suvr = 2.80   # Braak Stage VI saturation level

    # Patient trajectory
    if gamma_eff_per_year > 0:
        rate_k = min(1.2, 0.04 * gamma_eff_per_year)
        trajectory_patient = max_suvr / (1.0 + ((max_suvr - base_suvr) / base_suvr) * np.exp(-rate_k * years))
    else:
        clearance_rate = min(0.10, 0.005 * abs(gamma_eff_per_year))
        trajectory_patient = np.maximum(1.02, base_suvr * np.exp(-clearance_rate * years))

    # Optimized lifestyle intervention scenario (Sleep optimized to 8.5h + ALPS improved by +0.25 via vascular care)
    alps_opt = min(1.85, alps_index + 0.25)
    rates_opt = calculate_net_clearance_rate(
        alps_index=alps_opt,
        apoe_e4_alleles=apoe_e4_count,
        trem2_r47h=False,
        nlrp3_activated=False,
        sleep_fraction=8.5 / 24.0,
    )
    growth_opt = rates_opt["effective_growth_rate"] * (365.25 * 86400)
    if growth_opt > 0:
        rate_opt = min(1.2, 0.04 * growth_opt)
        trajectory_opt = max_suvr / (1.0 + ((max_suvr - base_suvr) / base_suvr) * np.exp(-rate_opt * years))
    else:
        clearance_opt = min(0.10, 0.005 * abs(growth_opt))
        trajectory_opt = np.maximum(1.02, base_suvr * np.exp(-clearance_opt * years))

    fig3, ax3 = plt.subplots(figsize=(10, 4.5), dpi=150)
    ax3.plot(years, trajectory_patient, color="#EF4444" if is_pathologic else "#3B82F6", linewidth=2.5, label="Current Patient Trajectory")
    ax3.plot(years, trajectory_opt, color="#10B981", linestyle="--", linewidth=2.0, label="Optimized Protocol (Sleep 8.5h + ALPS +0.25)")
    ax3.axhline(1.35, color="#F59E0B", linestyle=":", linewidth=1.5, label="Tau-PET Positivity Threshold (SUVR = 1.35)")

    ax3.set_xlabel("Follow-up Duration (Years)", fontsize=11, fontweight="bold")
    ax3.set_ylabel("Neocortical Tau-PET SUVR", fontsize=11, fontweight="bold")
    ax3.set_ylim(0.95, 3.00)
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.legend(loc="upper left", frameon=True, fontsize=9)
    plt.tight_layout()
    st.pyplot(fig3)
    plt.close(fig3)

    # Time to positivity calculation
    if is_pathologic:
        crossing_idx = np.where(trajectory_patient >= 1.35)[0]
        if len(crossing_idx) > 0:
            t_cross = years[crossing_idx[0]]
            st.warning(f"⚠️ **Clinical Alert**: Projected Tau-PET Positivity (SUVR $\\ge 1.35$) reached in **{t_cross:.1f} years** under current perivascular stasis.")
        else:
            st.warning("⚠️ **Clinical Alert**: Pathological accumulation detected; Tau-PET positivity imminent within 10 years.")
    else:
        st.success("🛡️ **Prognostic Assessment**: Convective clearance maintains long-term tau suppression below diagnostic PET thresholds (>10 years protected).")

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
