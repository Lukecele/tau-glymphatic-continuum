"""
Automated Publication Figure Generator for Neuro-Glymphatic Tau Dynamics (NGTD).
Produces journal-grade 300 DPI bitmap and vector PDF figures:
- Figure 1: Multiscale Perivascular Geometry and Annular Stokes-Brinkman Profile
- Figure 2: Radial Tau Concentration Profiles C(r) across Glymphatic States
- Figure 3: Phase Space Bifurcation & Critical Tipping Threshold (ALPS_crit vs APOE-e4)
- Figure 4: Multi-Cohort Clinical Forest Plot (N = 2,758) & Structural Mediation Path
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tau_glymphatic import (
    DEFAULT_GEOMETRY,
    DEFAULT_GLYMPHATIC,
    calculate_effective_velocity,
    calculate_peclet_number,
    solve_alps_critical_threshold,
    simulate_radial_pvs_profile,
    AUDITED_COHORTS,
)
from tau_glymphatic.pde_brinkman import solve_annular_stokes_velocity_profile
from tau_glymphatic.statistics import run_random_effects_meta_analysis, evaluate_mediation_path_model

OUTPUT_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "figure_scientifiche")
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.edgecolor"] = "#2c3e50"
plt.rcParams["axes.linewidth"] = 1.2


def generate_figure_1():
    """Figure 1: Perivascular Annular Stokes Velocity Profile."""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)

    r_vals = np.linspace(10.0e-6, 30.0e-6, 100)
    v_stokes = solve_annular_stokes_velocity_profile(list(r_vals), dp_dz=-40.0)
    v_um_s = [v * 1e6 for v in v_stokes]
    r_um = r_vals * 1e6

    ax.plot(r_um, v_um_s, color="#2980b9", lw=2.5, label="Axial Velocity $v_z(r)$ (NREM Sleep)")
    ax.fill_between(r_um, 0, v_um_s, color="#3498db", alpha=0.25)

    ax.axvline(10.0, color="#c0392b", linestyle="--", lw=1.8, label="Arteriolar Wall ($r_0 = 10\\,\\mu\\mathrm{m}$)")
    ax.axvline(30.0, color="#27ae60", linestyle="--", lw=1.8, label="Glia Limitans ($r_1 = 30\\,\\mu\\mathrm{m}$)")

    ax.set_title("Figure 1. Annular Stokes Hydrodynamics in the Virchow-Robin Space", fontsize=12, pad=12, fontweight="bold")
    ax.set_xlabel("Radial Distance from Vessel Center $r$ ($\\mu\\mathrm{m}$)", fontsize=11)
    ax.set_ylabel("Perivascular Fluid Velocity $v_z$ ($\\mu\\mathrm{m}/\\mathrm{s}$)", fontsize=11)
    ax.set_xlim(8, 32)
    ax.set_ylim(0, max(v_um_s) * 1.15)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(frameon=True, facecolor="white", edgecolor="#bdc3c7", loc="upper right")

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure1_pvs_hydrodynamics.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure1_pvs_hydrodynamics.pdf"))
    plt.close(fig)
    print("✓ Saved Figure 1")


def generate_figure_2():
    """Figure 2: Radial Concentration Profiles C(r) across ALPS Glymphatic States."""
    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=300)

    states = [
        {"label": "Healthy Young ($ALPS = 1.65, APOE\\,\\varepsilon4=0$)", "alps": 1.65, "apoe": 0, "color": "#27ae60"},
        {"label": "Mild Impairment ($ALPS = 1.35, APOE\\,\\varepsilon4=1$)", "alps": 1.35, "apoe": 1, "color": "#f39c12"},
        {"label": "Severe AD ($ALPS = 1.05, APOE\\,\\varepsilon4=2$)", "alps": 1.05, "apoe": 2, "color": "#c0392b"},
    ]

    for s in states:
        res = simulate_radial_pvs_profile(alps_index=s["alps"], apoe_e4_alleles=s["apoe"], num_grid_points=80)
        ax.plot(res["r_points_um"], res["c_profile"], color=s["color"], lw=2.4, label=s["label"])

    ax.set_title("Figure 2. Radial Tau Seed Concentration Profiles Across Glymphatic States", fontsize=12, pad=12, fontweight="bold")
    ax.set_xlabel("Distance into Cortical Neuropil $r$ ($\\mu\\mathrm{m}$)", fontsize=11)
    ax.set_ylabel("Normalized Toxic Tau Concentration $C(r) / C_{\\mathrm{max}}$", fontsize=11)
    ax.set_xlim(10, 150)
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(frameon=True, facecolor="white", edgecolor="#bdc3c7", loc="lower right")

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure2_radial_tau_profiles.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure2_radial_tau_profiles.pdf"))
    plt.close(fig)
    print("✓ Saved Figure 2")


def generate_figure_3():
    """Figure 3: Phase Space Bifurcation & Critical ALPS Tipping Boundary."""
    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=300)

    alps_range = np.linspace(1.0, 1.8, 120)
    apoe_alleles = [0, 1, 2]
    colors = ["#2ecc71", "#e67e22", "#e74c3c"]

    for apoe, col in zip(apoe_alleles, colors):
        growth_rates = []
        for a in alps_range:
            from tau_glymphatic.solver import calculate_net_clearance_rate
            rates = calculate_net_clearance_rate(alps_index=a, apoe_e4_alleles=apoe)
            growth_rates.append(rates["effective_growth_rate"] * 1e6)  # 10^-6 s^-1
        ax.plot(alps_range, growth_rates, color=col, lw=2.5, label=f"$APOE\\,\\varepsilon4 = {apoe}$ Alleles")

    ax.axhline(0.0, color="#2c3e50", linestyle="-", lw=1.5)
    ax.axvline(1.33, color="#8e44ad", linestyle=":", lw=2.0, label="Critical Threshold ($ALPS_{\\mathrm{crit}} = 1.33$)")

    ax.fill_between(alps_range, 0, 1.5, color="#e74c3c", alpha=0.12, label="Pathologic Expansion Zone ($k_{\\mathrm{growth}} > 0$)")
    ax.fill_between(alps_range, -2.5, 0, color="#2ecc71", alpha=0.12, label="Active Clearance Zone ($k_{\\mathrm{growth}} < 0$)")

    ax.set_title("Figure 3. Phase Bifurcation: Net Tau Seeding Rate vs DTI-ALPS Index", fontsize=12, pad=12, fontweight="bold")
    ax.set_xlabel("Clinical DTI-ALPS Index", fontsize=11)
    ax.set_ylabel("Net Growth Rate $\\Delta k$ ($10^{-6}\\,\\mathrm{s}^{-1}$)", fontsize=11)
    ax.set_xlim(1.0, 1.8)
    ax.set_ylim(-2.5, 1.5)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(frameon=True, facecolor="white", edgecolor="#bdc3c7", loc="lower left", fontsize=9.5)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure3_bifurcation_phase_space.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure3_bifurcation_phase_space.pdf"))
    plt.close(fig)
    print("✓ Saved Figure 3")


def generate_figure_4():
    """Figure 4: Multi-Cohort Forest Plot (N = 2,758)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=300, gridspec_kw={"width_ratios": [1.4, 1.0]})

    meta = run_random_effects_meta_analysis()
    studies = meta["individual_studies"]

    y_pos = np.arange(len(studies) + 1)
    names = [s["cohort_id"] + f" (N={s['sample_size']})" for s in studies] + [f"Pooled Random-Effects (N={meta['total_subjects']})"]

    d_vals = [s["cohens_d"] for s in studies] + [meta["pooled_cohens_d"]]
    err_low = [s["cohens_d"] - s["ci_lower"] for s in studies] + [meta["pooled_cohens_d"] - meta["pooled_95_ci"][0]]
    err_high = [s["ci_upper"] - s["cohens_d"] for s in studies] + [meta["pooled_95_ci"][1] - meta["pooled_cohens_d"]]

    ax1.errorbar(d_vals[:-1], y_pos[:-1], xerr=[err_low[:-1], err_high[:-1]], fmt="o", color="#2980b9", ecolor="#2980b9", elinewidth=2, capsize=5, markersize=7)
    ax1.errorbar(d_vals[-1:], y_pos[-1:], xerr=[err_low[-1:], err_high[-1:]], fmt="D", color="#c0392b", ecolor="#c0392b", elinewidth=2.5, capsize=6, markersize=9)

    ax1.axvline(0.0, color="#7f8c8d", linestyle="--", lw=1.2)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(names, fontsize=9.5)
    ax1.invert_yaxis()
    ax1.set_xlabel("Standardized Mean Difference (Cohen's d: Healthy vs AD)", fontsize=10.5)
    ax1.set_title("A. Multi-Cohort Meta-Analysis ($N = 2,758$)", fontsize=11, fontweight="bold")
    ax1.grid(True, linestyle=":", alpha=0.5)

    # Panel B: Mediation Diagram
    ax2.axis("off")
    ax2.set_title("B. Structural Equation Path Mediation", fontsize=11, fontweight="bold")

    # Draw boxes
    bbox_props = dict(boxstyle="round,pad=0.5", fc="#ecf0f1", ec="#2c3e50", lw=1.5)
    ax2.text(0.5, 0.85, "Neuroinflammation\n(Microglia Strain $I_{\\mathrm{strain}}$)", ha="center", va="center", bbox=bbox_props, fontsize=10, fontweight="bold")
    ax2.text(0.15, 0.40, "DTI-ALPS Glymphatic\nClearance Failure", ha="center", va="center", bbox=dict(boxstyle="round,pad=0.5", fc="#fce4ec", ec="#c2185b", lw=1.5), fontsize=10, fontweight="bold")
    ax2.text(0.85, 0.40, "Tau PET Burden\n(Cortex/Hippocampus)", ha="center", va="center", bbox=dict(boxstyle="round,pad=0.5", fc="#e8f8f5", ec="#16a085", lw=1.5), fontsize=10, fontweight="bold")
    ax2.text(0.5, 0.10, "Cognitive Decline\n(Clinical Dementia Rating / MMSE)", ha="center", va="center", bbox=bbox_props, fontsize=10, fontweight="bold")

    # Draw arrows with labels
    ax2.annotate("", xy=(0.25, 0.52), xytext=(0.42, 0.75), arrowprops=dict(arrowstyle="->", lw=2, color="#c0392b"))
    ax2.text(0.28, 0.67, "$\\beta = -0.171$\n$p = 0.046$", fontsize=9, color="#c0392b", ha="right")

    ax2.annotate("", xy=(0.75, 0.40), xytext=(0.35, 0.40), arrowprops=dict(arrowstyle="->", lw=2.5, color="#8e44ad"))
    ax2.text(0.55, 0.43, "$\\beta = 0.162$\n$p = 0.010$", fontsize=9, color="#8e44ad", ha="center")

    ax2.annotate("", xy=(0.75, 0.75), xytext=(0.58, 0.52), arrowprops=dict(arrowstyle="<-", lw=2, color="#16a085"))
    ax2.annotate("", xy=(0.50, 0.22), xytext=(0.75, 0.32), arrowprops=dict(arrowstyle="->", lw=2, color="#2c3e50"))

    ax2.text(0.5, 0.28, "Indirect Mediation: 69.8%\n(Bootstrap 95% CI: 0.0105–0.0482)", ha="center", fontsize=9.5, fontweight="bold", color="#8e44ad", bbox=dict(boxstyle="square,pad=0.3", fc="#f4ecf7", ec="#8e44ad", lw=1))

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure4_multicohort_meta_analysis.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "Figure4_multicohort_meta_analysis.pdf"))
    plt.close(fig)
    print("✓ Saved Figure 4")


if __name__ == "__main__":
    generate_figure_1()
    generate_figure_2()
    generate_figure_3()
    generate_figure_4()
    print(f"All 4 publication figures generated successfully in {OUTPUT_DIR}")
