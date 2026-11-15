"""
Command-Line Interface for Neuro-Glymphatic Tau Dynamics (NGTD).
Provides interactive biophysical risk assessment, DTI-ALPS translation, and audit execution.
"""

import sys
import argparse
import json
from tau_glymphatic.models import (
    DEFAULT_GLYMPHATIC,
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


def print_banner():
    banner = """
================================================================================
  🧠 NGTD v1.0.0 — Neuro-Glymphatic Tau Dynamics Calculator
  A Multiscale Biophysical Reaction-Diffusion & Perivascular Clearance Engine
  Author: Luca Celebrano, MD Candidate | University of Naples Federico II
================================================================================
"""
    print(banner.strip())


def main():
    parser = argparse.ArgumentParser(
        description="NGTD: Biophysical reaction-diffusion and perivascular glymphatic clearance calculator."
    )
    parser.add_argument(
        "--alps",
        type=float,
        default=DEFAULT_GLYMPHATIC.alps_healthy_ref,
        help="Patient's clinical DTI-ALPS index (e.g. 1.60 for healthy, 1.15 for AD).",
    )
    parser.add_argument(
        "--apoe-e4",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help="Number of APOE-e4 high-risk alleles (0, 1, or 2).",
    )
    parser.add_argument(
        "--trem2-r47h",
        action="store_true",
        help="Flag indicating carrier of TREM2 R47H loss-of-function variant.",
    )
    parser.add_argument(
        "--nlrp3",
        action="store_true",
        help="Flag indicating chronic central or peripheral NLRP3 inflammasome hyper-activation.",
    )
    parser.add_argument(
        "--audit",
        action="store_true",
        help="Execute the full 9-step biophysical invariant and dimensional audit suite.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw structured JSON for integration into clinical pipelines.",
    )

    args = parser.parse_args()

    if args.audit:
        print_banner()
        print("\n🔬 Executing Full Biophysical Invariant & Dimensional Audit Battery...\n")
        audit_res = run_full_biophysical_audit()
        for c in audit_res["checks"]:
            mark = "✓" if c["passed"] else "✗"
            print(f"  [{mark}] {c['check']}: {c['detail']}")
        print(f"\nAudit Conclusion: {audit_res['status']} ({audit_res['passed_checks']}/{audit_res['total_checks']} Passed)")
        sys.exit(0 if audit_res["status"] == "PASSED" else 1)

    # Perform patient simulation
    rates = calculate_net_clearance_rate(
        alps_index=args.alps,
        apoe_e4_alleles=args.apoe_e4,
        trem2_r47h=args.trem2_r47h,
        nlrp3_activated=args.nlrp3,
    )
    peclet_sleep = calculate_peclet_number(args.alps, peak_sleep=True)
    peclet_avg = calculate_peclet_number(args.alps, peak_sleep=False)
    v_eff = calculate_effective_velocity(args.alps)
    crit_info = solve_alps_critical_threshold()

    data = {
        "inputs": {
            "dti_alps_index": args.alps,
            "apoe_e4_alleles": args.apoe_e4,
            "trem2_r47h": args.trem2_r47h,
            "nlrp3_activated": args.nlrp3,
        },
        "biophysical_diagnostics": {
            "perivascular_velocity_nm_s": round(v_eff * 1e9, 2),
            "peclet_number_sleep": round(peclet_sleep, 4),
            "peclet_number_24h": round(peclet_avg, 4),
            "alps_critical_tipping_threshold": crit_info["alps_crit_kinetic"],
            "microglial_strain_index": round(rates["microglial_strain"], 3),
            "net_growth_rate_1_s": rates["effective_growth_rate"],
            "net_stasis_ratio": round(rates["net_stasis_ratio"], 3),
            "regime": "PATHOLOGIC EXPONENTIAL EXPANSION (STASIS NICHES)" if rates["is_pathologic_expansion"] else "ACTIVE GLYMPHATIC CLEARANCE (EQUILIBRIUM)",
        },
    }

    if args.json:
        print(json.dumps(data, indent=2))
        return

    print_banner()
    print("\n📋 CLINICAL & BIOPHYSICAL INPUTS:")
    print(f"   • DTI-ALPS Index:               {args.alps:.2f} (Reference: Healthy ≥ 1.50, AD ≤ 1.25)")
    print(f"   • APOE-ε4 Allele Burden:        {args.apoe_e4} copies")
    print(f"   • TREM2 R47H Carrier:           {'YES (Loss of Function)' if args.trem2_r47h else 'NO (Normal Phagocytosis)'}")
    print(f"   • NLRP3 Inflammasome Status:    {'HYPER-ACTIVATED' if args.nlrp3 else 'BASAL RESTING'}")

    print("\n⚡ BIOPHYSICAL CLEARANCE DYNAMICS:")
    print(f"   • 24h Interstitial Velocity:    {data['biophysical_diagnostics']['perivascular_velocity_nm_s']} nm/s")
    print(f"   • PVS Péclet Number (Sleep):    {data['biophysical_diagnostics']['peclet_number_sleep']} ({'Convection Dominates' if peclet_sleep > 1.0 else 'Diffusion & Stasis Dominate'})")
    print(f"   • Critical ALPS Tipping Point:  {crit_info['alps_crit_kinetic']} (Below this value, washout fails)")
    print(f"   • Microglial Strain Index:      {data['biophysical_diagnostics']['microglial_strain_index']} / 1.000")
    print(f"   • Solute Stasis Ratio:          {data['biophysical_diagnostics']['net_stasis_ratio']} (Aggregation / Clearance)")

    status_color = "🚨" if rates["is_pathologic_expansion"] else "🛡️"
    print(f"\n{status_color} PHENOTYPE REGIME: {data['biophysical_diagnostics']['regime']}")
    if rates["is_pathologic_expansion"]:
        print("   ⚠️  WARNING: Glymphatic convective clearance is insufficient to counteract autocatalytic")
        print("       Tau seeding. High risk of neocortical Braak stage advancement.")
    else:
        print("   ✅ STABLE: Convective perivascular bulk flow maintains clearance equilibrium.")
        print("       Tau oligomers are actively purged before cross-synaptic spreading fronts stabilize.")
    print()


if __name__ == "__main__":
    main()
