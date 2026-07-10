#!/usr/bin/env python3
"""
Reproduce M-Star Fig. 25 with a minimal three-line comparison.

https://docs.mstarcfd.com/1b_HowToGuides/predicting-mass-transfer-kLA.html

The plot shows tank-averaged dissolved-O2 concentration ⟨C_L⟩ (mol/L) vs
time (s) on M-Star Fig. 25 axes (x ∈ [0, 35], y ∈ [0, 8e-5]).  Three
model curves are overlaid on top of our simulated ⟨C_L⟩(t) trace:

    1. M-Star Fig. 25 linear fit
           y = 1.5e-6 * t + 2.9e-6      => k_La = 3.78 /hr
       (from the M-Star documentation, valid over t ∈ [10, 35] s)

    2. Our method [2] log-linear fit over the WHOLE simulation
           ln(C_sat - <C_L>) = ln(C_sat) - k_La * t
       (uses every sample from t ≥ 0.2 s to t_end; single k_La per run)

    3. Van't Riet (1979) theoretical estimate for coalescing (aqueous)
       sparged agitated tanks
           k_La = 0.026 * (P/V)^0.4 * u_g^0.5      [1/s]
       with P/V = Np * ρ * N³ * D⁵ / V, u_g = Q / A_tank, and Np = 1.5.
       Van't Riet is well known to over-predict at low aeration rates
       (~2x here), so this line is included as an upper-bound reference,
       not as a competing "best guess".
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from plot_kla import parse_bubble_stats, C_SAT_MOL_M3  # noqa: E402


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
# Primary run for the fit — the completed nu=0.016 (SGS off) benchmark run,
# whose data was archived to resultsJuly8/ on 2026-07-08 when we started
# the halved-viscosity experiment in the same directory.
CSV_PRIMARY   = SCRIPT_DIR / "resultsJuly8" / "bubble_stats.csv"
LABEL_PRIMARY = "kla_bioreactor  (nu=0.016, no SGS, C_eq=1.427, resultsJuly8)"

# Secondary run for overlay comparison.
CSV_SECONDARY   = SCRIPT_DIR.parent / "kla_bioreactor_surface_bc" / "bubble_stats.csv"
LABEL_SECONDARY = "kla_bioreactor_surface_bc  (nu=0.016, no SGS, C_eq=0.275)"

# Current run — the live LES benchmark with SGS in collision + surface BC.
# Same nu=0.016 as the primary/secondary but adds Smagorinsky eddy viscosity
# in the BGK relaxation and the desorption-only surface O2 BC.
CSV_TERTIARY   = SCRIPT_DIR / "bubble_stats.csv"
LABEL_TERTIARY = "kla_bioreactor  (nu=0.016, SGS on, C_eq=0.275, live)"

# M-Star Fig. 22 (probe-based constrained exponential fit).
MSTAR_KLA_FIG22_HR = 4.1        # /hr — the "headline" M-Star benchmark

# M-Star Fig. 25 published linear fit.
MSTAR_SLOPE_MOL_L_S   = 1.5e-6
MSTAR_INTERCEPT_MOL_L = 2.9e-6
MSTAR_KLA_FIG25_HR    = MSTAR_SLOPE_MOL_L_S / (C_SAT_MOL_M3 * 1e-3) * 3600

# Applikon experimental measurement quoted in the M-Star docs
# ("which agrees well with the experimentally measured value of 4.5 1/hr").
# M-Star does not cite the source of this number, so we do not name it.
EXPERIMENTAL_KLA_HR   = 4.5

# Zakrzewski (2020) experimental k_La value cited in Thomas et al. 2021
# (this benchmark's primary reference).  Same benchtop bioreactor
# geometry (5 L, 400 RPM, 0.4 L/min sparge) but a formally published
# measurement whose provenance we CAN cite in a paper.
ZAKRZEWSKI_KLA_HR     = 5.1

# Van't Riet (1979) — coalescing/aqueous form:
#     k_La = 0.026 * (P/V)^0.4 * u_g^0.5     [1/s]
# with P computed from the impeller power number.
VR_NP       = 1.5           # power number for pitched-blade turbine (M-Star docs)
VR_RHO      = 1000.0        # kg/m^3
VR_N_RPS    = 400.0 / 60.0  # rev/s  (400 RPM)
VR_D_IMP_M  = 0.058         # m      (impeller diameter)
VR_T_TANK_M = 0.18          # m      (tank diameter)
VR_V_LIQ_M3 = 3.9e-3        # m^3    (fluid volume)
VR_Q_M3_S   = 0.4 / 60000.0 # 0.4 L/min = 6.67e-6 m^3/s

def vant_riet_kla_hr():
    A_tank = np.pi * (VR_T_TANK_M / 2.0) ** 2
    u_g    = VR_Q_M3_S / A_tank                              # m/s
    P      = VR_NP * VR_RHO * VR_N_RPS**3 * VR_D_IMP_M**5    # W
    P_over_V = P / VR_V_LIQ_M3                               # W/m^3
    kla_s = 0.026 * (P_over_V ** 0.4) * (u_g ** 0.5)         # 1/s
    return kla_s * 3600.0, P_over_V, u_g

VR_KLA_HR, VR_POV, VR_UG = vant_riet_kla_hr()

# Fig. 25 axes
XMAX_S = 35.0
YMAX_MOL_L = 8.0e-5

OUT_PNG = SCRIPT_DIR / "kla_fig25.png"


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main() -> int:
    C_sat_mol_L = C_SAT_MOL_M3 * 1e-3

    def load_and_fit(path):
        """Load run, return (t, C_mol_L, method2 fit, method1 fit)."""
        b = parse_bubble_stats(path)
        t = b["phys_time_s"]
        C_mol_L = b["C_L_mol_m3"] * 1e-3
        keep = np.isfinite(C_mol_L) & np.isfinite(t)
        t, C_mol_L = t[keep], C_mol_L[keep]
        # Method [2]: log-linear regression on ln(C_sat - C_L) vs t.
        m = (t >= 0.2) & (C_mol_L > 0) & (C_mol_L < C_sat_mol_L)
        b_ll, a_ll = np.polyfit(t[m], np.log(C_sat_mol_L - C_mol_L[m]), 1)
        kla2_s = -b_ll
        # Method [1]: constrained exponential fit
        #   C_L(t) = C_sat * (1 - exp(-kLa*t))
        # single-parameter least-squares on the raw C_L(t) values.
        model = lambda tt, k: C_sat_mol_L * (1.0 - np.exp(-k * tt))
        popt, _ = curve_fit(model, t[m], C_mol_L[m],
                             p0=[kla2_s], bounds=([1e-6], [1.0]),
                             maxfev=20000)
        kla1_s = float(popt[0])
        return t, C_mol_L, kla1_s, kla1_s*3600, kla2_s, kla2_s*3600

    t1, C1, kla1a_s, kla1a_hr, kla1b_s, kla1b_hr = load_and_fit(CSV_PRIMARY)
    t2, C2, kla2a_s, kla2a_hr, kla2b_s, kla2b_hr = load_and_fit(CSV_SECONDARY)
    # Tertiary is loaded only if the CSV exists AND has at least a few
    # samples past 0.2 s (the log-linear fit's start time).  Otherwise we
    # print a warning and skip the third series so the plot still works
    # from a fresh run directory.
    have_tertiary = False
    if CSV_TERTIARY.exists():
        try:
            t3, C3, kla3a_s, kla3a_hr, kla3b_s, kla3b_hr = load_and_fit(CSV_TERTIARY)
            have_tertiary = (len(t3) > 20 and t3[-1] > 1.0)
            if not have_tertiary:
                print(f"NOTE: {CSV_TERTIARY.name} exists but has <20 samples or "
                      f"t_max <= 1 s; skipping tertiary curve.")
        except Exception as e:
            print(f"NOTE: could not fit {CSV_TERTIARY.name} ({e}); skipping.")
    else:
        print(f"NOTE: {CSV_TERTIARY} not found; skipping tertiary curve.")

    print(f"C_sat                 = {C_SAT_MOL_M3:.4f} mol/m^3 = {C_sat_mol_L:.4e} mol/L")
    print(f"M-Star Fig. 22 fit    : k_La = {MSTAR_KLA_FIG22_HR:.3f} /hr  (probe constrained exp — 'headline')")
    print(f"M-Star Fig. 25 fit    : k_La = {MSTAR_KLA_FIG25_HR:.3f} /hr  "
          f"(y = {MSTAR_SLOPE_MOL_L_S:.2e}*t + {MSTAR_INTERCEPT_MOL_L:.2e})")
    print(f"Experiment (M-Star)   : k_La = {EXPERIMENTAL_KLA_HR:.3f} /hr  (unattributed in the docs)")
    print(f"Experiment (Zakrzewski 2020) : k_La = {ZAKRZEWSKI_KLA_HR:.3f} /hr  (cited by Thomas 2021)")
    print(f"Primary   ({LABEL_PRIMARY})")
    print(f"  method [1] (constrained exp) : k_La = {kla1a_hr:.3f} /hr  "
          f"({kla1a_hr/MSTAR_KLA_FIG22_HR:.2f}x M-Star Fig.22)")
    print(f"  method [2] (log-linear)      : k_La = {kla1b_hr:.3f} /hr  "
          f"({kla1b_hr/MSTAR_KLA_FIG25_HR:.2f}x M-Star Fig.25)")
    print(f"Secondary ({LABEL_SECONDARY})")
    print(f"  method [1] (constrained exp) : k_La = {kla2a_hr:.3f} /hr  "
          f"({kla2a_hr/MSTAR_KLA_FIG22_HR:.2f}x M-Star Fig.22)")
    print(f"  method [2] (log-linear)      : k_La = {kla2b_hr:.3f} /hr  "
          f"({kla2b_hr/MSTAR_KLA_FIG25_HR:.2f}x M-Star Fig.25)")
    if have_tertiary:
        print(f"Tertiary  ({LABEL_TERTIARY})")
        print(f"  t_max                        = {t3[-1]:.2f} s  ({len(t3)} samples)")
        print(f"  method [1] (constrained exp) : k_La = {kla3a_hr:.3f} /hr  "
              f"({kla3a_hr/MSTAR_KLA_FIG22_HR:.2f}x M-Star Fig.22)")
        print(f"  method [2] (log-linear)      : k_La = {kla3b_hr:.3f} /hr  "
              f"({kla3b_hr/MSTAR_KLA_FIG25_HR:.2f}x M-Star Fig.25)")

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.6), constrained_layout=True,
                             sharey=True)
    tfit_mstar = np.linspace(10.0, XMAX_S, 200)
    tfull = np.linspace(0.0, XMAX_S, 400)

    # Panel specs: each panel picks the M-Star reference that matches its
    # own fit convention (exponential vs linear-slope), plus the experimental
    # 4.5 /hr line shown as an exponential on both panels for direct visual
    # comparison of the underlying asymptotic behaviour.
    panel_specs = [
        {
            "title":  "Method [2] — log-linear on  ln(C_sat − C_L)  vs t",
            "our_kla_s":  kla1b_s,  "our_kla_hr":  kla1b_hr,
            "sec_kla_s":  kla2b_s,  "sec_kla_hr":  kla2b_hr,
            "ter_kla_s":  (kla3b_s  if have_tertiary else None),
            "ter_kla_hr": (kla3b_hr if have_tertiary else None),
            "mstar_kind": "linear",   # matches Fig. 25 (linear slope 1.5e-6)
        },
        {
            "title":  "Method [1] — constrained exp:  C_sat·(1 − exp(−k_La·t))",
            "our_kla_s":  kla1a_s,  "our_kla_hr":  kla1a_hr,
            "sec_kla_s":  kla2a_s,  "sec_kla_hr":  kla2a_hr,
            "ter_kla_s":  (kla3a_s  if have_tertiary else None),
            "ter_kla_hr": (kla3a_hr if have_tertiary else None),
            "mstar_kind": "exp_fig22",   # matches Fig. 22 (exp fit, 4.1 /hr)
        },
    ]

    for ax, sp in zip(axes, panel_specs):
        # Data traces.
        ax.plot(t1, C1, "-", color="0.75", lw=1.0,
                label=f"⟨C_L⟩(t)  primary  (t up to {t1[-1]:.1f} s)")
        ax.plot(t2, C2, "-", color="#f7b2a8", lw=1.0,
                label=f"⟨C_L⟩(t)  surface_bc  (t up to {t2[-1]:.1f} s)")
        if have_tertiary:
            ax.plot(t3, C3, "-", color="#98df8a", lw=1.0,
                    label=f"⟨C_L⟩(t)  LES live  (t up to {t3[-1]:.1f} s)")

        # M-Star reference (per panel).
        if sp["mstar_kind"] == "linear":
            ax.plot(tfit_mstar,
                    MSTAR_SLOPE_MOL_L_S * tfit_mstar + MSTAR_INTERCEPT_MOL_L,
                    "-", color="black", lw=2.0,
                    label=f"M-Star Fig. 25 (linear fit):  k_La = {MSTAR_KLA_FIG25_HR:.2f} /hr")
        else:
            ax.plot(tfull,
                    C_sat_mol_L * (1.0 - np.exp(-MSTAR_KLA_FIG22_HR/3600.0 * tfull)),
                    "-", color="black", lw=2.0,
                    label=f"M-Star Fig. 22 (exp fit):  k_La = {MSTAR_KLA_FIG22_HR:.2f} /hr")

        # Experimental references — plotted on both panels as constrained
        # exponential curves at the reported k_La so they compare directly
        # to our method [1] fits.  Two experiments:
        #   * unattributed 4.5 /hr from the M-Star how-to page (no source
        #     cited in the docs, so we leave it unattributed here too)
        #   * Zakrzewski (2020) 5.1 /hr cited by Thomas et al. 2021 for the
        #     same 5 L benchtop bioreactor geometry at 400 RPM / 0.4 L/min.
        ax.plot(tfull,
                C_sat_mol_L * (1.0 - np.exp(-EXPERIMENTAL_KLA_HR/3600.0 * tfull)),
                "--", color="#7f7f7f", lw=1.6,
                label=f"experiment (M-Star docs):  k_La = {EXPERIMENTAL_KLA_HR:.2f} /hr")
        ax.plot(tfull,
                C_sat_mol_L * (1.0 - np.exp(-ZAKRZEWSKI_KLA_HR/3600.0 * tfull)),
                "--", color="#9467bd", lw=1.6,
                label=f"experiment (Zakrzewski 2020):  k_La = {ZAKRZEWSKI_KLA_HR:.2f} /hr")

        # Our fits (primary + secondary always; tertiary if available).
        ax.plot(tfull, C_sat_mol_L * (1.0 - np.exp(-sp["our_kla_s"] * tfull)),
                "-", color="#1f77b4", lw=2.2,
                label=f"primary fit:  k_La = {sp['our_kla_hr']:.2f} /hr")
        ax.plot(tfull, C_sat_mol_L * (1.0 - np.exp(-sp["sec_kla_s"] * tfull)),
                "-", color="#d62728", lw=2.2,
                label=f"surface_bc fit:  k_La = {sp['sec_kla_hr']:.2f} /hr")
        if have_tertiary:
            ax.plot(tfull, C_sat_mol_L * (1.0 - np.exp(-sp["ter_kla_s"] * tfull)),
                    "-", color="#2ca02c", lw=2.2,
                    label=f"LES live fit:  k_La = {sp['ter_kla_hr']:.2f} /hr")

        ax.set_title(sp["title"], fontsize=10)
        ax.set_xlabel("Time  (s)")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
        ax.set_xlim(0.0, XMAX_S)
        ax.set_ylim(0.0, YMAX_MOL_L)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    axes[0].set_ylabel(r"$\langle C_L\rangle$   Tank-average DO concentration  (mol/L)")
    fig.suptitle("k_La comparison — M-Star Fig. 25 axes  (methods [1] and [2] side-by-side)",
                 fontsize=11)

    fig.savefig(OUT_PNG, dpi=140)
    print(f"\nsaved -> {OUT_PNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
