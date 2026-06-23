#!/usr/bin/env python3
"""
k_La verification for the kLa-bioreactor run against the M-Star benchmark
https://docs.mstarcfd.com/1b_HowToGuides/predicting-mass-transfer-kLA.html
target k_La ≈ 4.1 /hr at 400 RPM, 0.4 L/min air, dilute aqueous O₂.

Inputs (all relative to this script's directory):
  resultsJune21full/marbles_14431463.out   — main run log (steps 0..534400)
  resultsJune21full/marbles_14516702.out   — restart run log (steps 520000+)
  resultsJune21full/bubble_stats.csv       — bubble-side mass balance

Approach
========

(A) rho_O2_max trajectory (cheap, qualitative).
    The .out file logs `[O2_debug step=N] rho_O2_after=X` where X is the
    max-norm of the dissolved-O2 component lattice (peak cell value).
    This is NOT a volume average, but its time-constant matches kLa
    because every cell relaxes toward the same saturation under the
    same k_L·a → max(C_L) follows the same exponential shape.
    Fit:  C_max(t) = C_inf × (1 - exp(-kLa × t))   (free C_inf, kLa)

(B) Mean dissolved-O2 from bubble mass balance (quantitative cross-check).
    Conservation:
      n_O2_dissolved(t) = n_O2_injected(t) - n_O2_in_bubbles(t)
                         - n_O2_vented_out(t)
    With no surface model for venting, treat (n_O2_injected - n_O2_in_bubbles)
    as an upper bound on dissolved.  Plot that envelope and compare to
    the theoretical exponential trajectory.

Both convert LB → SI using parameters from kla_bioreactor.inp:
  dt_phys     = 1.19365e-5 s/step    (omega_LB / omega_phys = 0.0005/41.888)
  C_ref       = 100.0 mol/m³ per LB-rho unit
  V_liquid    = 0.16 m × π × (0.09 m)² = 4.07e-3 m³  (tank radius 0.09 m,
                liquid height 0.16 m)
  C_sat       = S × C_g = 0.032 × 44.6 = 1.428 mol/m³  (Henry sat at STP)
  injection   = 0.4 L/min = 6.67e-6 m³/s gas at STP
                ⇒ n_inj_rate = 6.67e-6 × 44.6 = 2.974e-4 mol/s
                ⇒ d(n_inj)/dstep = 2.974e-4 × dt_phys = 3.551e-9 mol/step

The M-Star benchmark page gives k_La ≈ 4.1 /hr = 1.139e-3 /s.
Marbles target the same value; this script reports the fitted k_La with
ASCII confidence interval and a residual plot.
"""

import csv
import math
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# --------------------------------------------------------------------------
# 1. CONFIG — keep in sync with kla_bioreactor.inp
# --------------------------------------------------------------------------
DT_PHYS     = 1.19365e-5     # s / LB step
C_REF       = 100.0          # mol/m³ per LB-rho unit (lbm.O2_concentration_reference)
C_SAT_MOL_M3 = 0.032 * 44.6  # = 1.428 mol/m³  Henry-saturated liquid
V_LIQUID_M3 = math.pi * (0.09 ** 2) * 0.16   # ≈ 4.07e-3 m³  (cylindrical tank,
                                              # H=0.16 m, R=T/2=0.09 m — neglects
                                              # impeller-/sparger-occupied volume)
N_INJ_RATE_MOL_S = 6.67e-6 * 44.6            # 2.974e-4 mol/s O₂ injected
KLA_BENCHMARK_S  = 4.1 / 3600.0              # 1.139e-3 /s (M-Star reference)

RESULTS_DIR = Path(__file__).resolve().parent / "resultsJune21full"
OUTFILES = [
    RESULTS_DIR / "marbles_14431463.out",
    RESULTS_DIR / "marbles_14516702.out",
]
BUBBLE_CSV = RESULTS_DIR / "bubble_stats.csv"


# --------------------------------------------------------------------------
# 2. Parse [O2_debug step=N] rho_O2_after=X across both log files; stitch.
# --------------------------------------------------------------------------
def parse_rho_o2_max():
    """Return arrays (step, rho_O2_max_LB) from concatenated logs.

    Stitching strategy: keep run-A entries up to its last step, then drop
    any run-B entries with step <= that boundary (handles the overlap
    between the original run and the restart cleanly).
    """
    pat = re.compile(r"\[O2_debug step=(\d+)\]\s+o2_src\.norm0=\S+\s+rho_O2_before=(\S+)\s*$", re.M)
    pat_after = re.compile(r"\[O2_debug step=(\d+)\]\s+rho_O2_after=(\S+)\s*$", re.M)

    by_step = {}
    for path in OUTFILES:
        if not path.exists():
            print(f"warn: missing {path}", file=sys.stderr)
            continue
        text = path.read_text(errors="replace")
        for step_s, rho_s in pat_after.findall(text):
            step = int(step_s)
            rho = float(rho_s)
            # If a duplicate step appears (overlap), the later log wins —
            # but only if its value is finite; restart-init artefacts can
            # produce zeros, so we trust the original log for shared steps.
            if step in by_step:
                continue
            by_step[step] = rho

    steps = np.array(sorted(by_step.keys()))
    rhos  = np.array([by_step[s] for s in steps])
    return steps, rhos


# --------------------------------------------------------------------------
# 3. Parse bubble_stats.csv
# --------------------------------------------------------------------------
def parse_bubble_stats():
    """Return dict-of-arrays from bubble_stats.csv.

    Handles two schemas:
      legacy 8-column : step, phys_time_s, n_bubbles, d_mean_mm, d_min_mm,
                        d_max_mm, n_O2_total_mol, dn_O2_step_mol_per_s
      new   10-column : ... + C_L_mol_m3, V_liq_m3
    Missing columns are absent from the returned dict so callers can
    detect-and-branch on them.
    """
    cols = {}
    with open(BUBBLE_CSV, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for k, v in row.items():
                if k is None or v is None or v == "":
                    continue
                try:
                    val = float(v)
                except ValueError:
                    continue
                cols.setdefault(k, []).append(val)
    return {k: np.array(v) for k, v in cols.items()}


# --------------------------------------------------------------------------
# 4. k_La model & fit
# --------------------------------------------------------------------------
def C_t(t, C_inf, kLa):
    """C(t) = C_inf · (1 − exp(−kLa·t))  with C(0) = 0."""
    return C_inf * (1.0 - np.exp(-kLa * t))


def fit_kla(t_s, C_mol_m3, fit_start_s=0.5):
    """Fit C_t against the post-spin-up portion of the trajectory.

    The first ~0.5 s is the impeller ramp + bubble inventory build-up,
    not the kLa-controlled mass-transfer phase.
    """
    mask = t_s >= fit_start_s
    if mask.sum() < 5:
        raise RuntimeError("not enough points after fit_start_s")
    t = t_s[mask]
    C = C_mol_m3[mask]
    # Sensible initial guess
    p0 = (max(C.max(), C_SAT_MOL_M3), KLA_BENCHMARK_S)
    bounds = ([1e-4, 1e-6], [10 * C_SAT_MOL_M3, 1.0])
    popt, pcov = curve_fit(C_t, t, C, p0=p0, bounds=bounds, maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr, t, C


# --------------------------------------------------------------------------
# 5. Main
# --------------------------------------------------------------------------
def main():
    print("=== k_La verification: kla_bioreactor vs M-Star benchmark ===")
    print(f"  benchmark target : k_La = 4.1 /hr = {KLA_BENCHMARK_S:.4e} /s")
    print(f"  C_sat (Henry)    : {C_SAT_MOL_M3:.3f} mol/m³")
    print(f"  V_liquid (est.)  : {V_LIQUID_M3*1000:.2f} L")
    print(f"  n_inj_rate       : {N_INJ_RATE_MOL_S:.3e} mol/s")
    print()

    # ---- (A) rho_O2_max trajectory ----
    steps, rho_max_LB = parse_rho_o2_max()
    if len(steps) == 0:
        sys.exit("no rho_O2_after entries found in logs")
    t_s = steps * DT_PHYS
    C_max = rho_max_LB * C_REF                  # mol/m³

    print(f"[A] rho_O2_max  : {len(steps)} samples, t = [{t_s[0]:.3f}, {t_s[-1]:.3f}] s")
    print(f"    max(rho_O2_max_LB) = {rho_max_LB.max():.4e} → {C_max.max():.3f} mol/m³")

    try:
        (C_inf_max, kLa_max), (err_C, err_kLa), tA, CA = fit_kla(t_s, C_max)
        kLa_max_hr = kLa_max * 3600.0
        err_hr     = err_kLa * 3600.0
        ratio_max  = kLa_max_hr / 4.1
        print(f"    fit: C_inf = {C_inf_max:.3f} ± {err_C:.3f} mol/m³")
        print(f"         k_La  = {kLa_max:.3e} ± {err_kLa:.1e} /s  =  {kLa_max_hr:.2f} ± {err_hr:.2f} /hr")
        print(f"         k_La / k_La_benchmark = {ratio_max:.2f}× (benchmark = 4.1 /hr)")
    except Exception as ex:
        print(f"    fit failed: {ex}")
        C_inf_max = kLa_max = float("nan")
        kLa_max_hr = ratio_max = float("nan")

    # ---- (B) Mean dissolved-O2 envelope from bubble mass balance ----
    print()
    print("[B] bubble mass balance:")
    b = parse_bubble_stats()
    b_t = b["phys_time_s"]
    b_n_total = b["n_O2_total_mol"]   # moles still inside all bubbles
    n_inj_cum = N_INJ_RATE_MOL_S * b_t  # cumulative moles injected (linear in t)
    # Upper bound on dissolved (assumes nothing has vented through the surface):
    n_dissolved_upper = n_inj_cum - b_n_total
    n_dissolved_upper = np.maximum(n_dissolved_upper, 0.0)
    C_mean_upper = n_dissolved_upper / V_LIQUID_M3  # mol/m³

    print(f"    n_inj_cum(end)   = {n_inj_cum[-1]:.4e} mol")
    print(f"    n_total_bub(end) = {b_n_total[-1]:.4e} mol")
    print(f"    n_diss_upper(end)= {n_dissolved_upper[-1]:.4e} mol")
    print(f"    <C_L>_upper(end) = {C_mean_upper[-1]:.4f} mol/m³  ({100*C_mean_upper[-1]/C_SAT_MOL_M3:.1f}% of C_sat)")

    try:
        (C_inf_avg, kLa_avg), (err_C2, err_kLa2), tB, CB = fit_kla(b_t, C_mean_upper)
        kLa_avg_hr = kLa_avg * 3600.0
        err_hr_avg = err_kLa2 * 3600.0
        ratio_avg  = kLa_avg_hr / 4.1
        print(f"    fit: C_inf = {C_inf_avg:.3f} ± {err_C2:.3f} mol/m³")
        print(f"         k_La  = {kLa_avg:.3e} ± {err_kLa2:.1e} /s  =  {kLa_avg_hr:.2f} ± {err_hr_avg:.2f} /hr")
        print(f"         k_La / k_La_benchmark = {ratio_avg:.2f}× (benchmark = 4.1 /hr)")
    except Exception as ex:
        print(f"    fit failed: {ex}")
        C_inf_avg = kLa_avg = float("nan")
        kLa_avg_hr = ratio_avg = float("nan")

    # ---- (C) Direct C_L from bubble_stats.csv (preferred when present) ----
    # If the CSV has the new schema (commit XXX), it contains the
    # volume-averaged dissolved-O₂ concentration <C_L>_mol_m3 computed
    # on-the-fly by the LBM solver over CELL_LIQUID + φ·CELL_INTERFACE.
    # This is the cleanest source for a kLa fit — no max-norm bias, no
    # bubble-balance closure assumption.  Only run this branch when the
    # column exists and contains usable (non-NaN) data.
    have_direct_CL = ("C_L_mol_m3" in b) and np.isfinite(b["C_L_mol_m3"]).any()
    print()
    if have_direct_CL:
        b_CL = b["C_L_mol_m3"]
        mask = np.isfinite(b_CL) & (b_CL > 0)
        t_C  = b_t[mask]
        C_C  = b_CL[mask]
        print(f"[C] direct <C_L> from bubble_stats.csv: {len(t_C)} samples")
        print(f"    <C_L>(end) = {C_C[-1]:.4f} mol/m³  ({100*C_C[-1]/C_SAT_MOL_M3:.1f}% of C_sat)")
        try:
            (C_inf_dir, kLa_dir), (err_C3, err_kLa3), tC, CC = fit_kla(t_C, C_C)
            kLa_dir_hr = kLa_dir * 3600.0
            err_hr_dir = err_kLa3 * 3600.0
            ratio_dir  = kLa_dir_hr / 4.1
            print(f"    fit: C_inf = {C_inf_dir:.3f} ± {err_C3:.3f} mol/m³")
            print(f"         k_La  = {kLa_dir:.3e} ± {err_kLa3:.1e} /s  =  {kLa_dir_hr:.2f} ± {err_hr_dir:.2f} /hr")
            print(f"         k_La / k_La_benchmark = {ratio_dir:.2f}× (benchmark = 4.1 /hr)")
        except Exception as ex:
            print(f"    fit failed: {ex}")
            C_inf_dir = kLa_dir = float("nan")
            kLa_dir_hr = ratio_dir = float("nan")
    else:
        print("[C] direct <C_L> column not in CSV; skipping (old-schema bubble_stats.csv)")
        C_inf_dir = kLa_dir = float("nan")
        kLa_dir_hr = ratio_dir = float("nan")

    # ---- Plot ----
    n_panels = 3 if have_direct_CL else 2
    fig, axes = plt.subplots(n_panels, 1, figsize=(9, 4 * n_panels), sharex=True)

    ax = axes[0]
    ax.plot(t_s, C_max, ".", ms=2, alpha=0.4, label="max C_L (lattice norm0)")
    if math.isfinite(kLa_max):
        t_fit = np.linspace(0, t_s[-1], 400)
        ax.plot(t_fit, C_t(t_fit, C_inf_max, kLa_max), "-",
                label=f"fit: k_La = {kLa_max_hr:.2f} /hr (C∞ = {C_inf_max:.2f})")
    ax.axhline(C_SAT_MOL_M3, ls="--", color="gray", label=f"C_sat = {C_SAT_MOL_M3:.2f} mol/m³")
    # Reference: theoretical M-Star curve toward C_sat with k_La = 4.1/hr
    t_ref = np.linspace(0, t_s[-1], 400)
    ax.plot(t_ref, C_t(t_ref, C_SAT_MOL_M3, KLA_BENCHMARK_S), ":", color="k",
            label="M-Star: k_La=4.1/hr, C_inf=C_sat")
    ax.set_ylabel("max C_L  (mol/m³)")
    ax.set_title("(A) peak dissolved-O₂ trajectory (lattice max-norm)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")

    ax = axes[1]
    ax.plot(b_t, C_mean_upper, ".-", ms=2, label="<C_L> upper bound from bubble mass balance")
    if math.isfinite(kLa_avg):
        t_fit = np.linspace(0, b_t[-1], 400)
        ax.plot(t_fit, C_t(t_fit, C_inf_avg, kLa_avg), "-",
                label=f"fit: k_La = {kLa_avg_hr:.2f} /hr (C∞ = {C_inf_avg:.2f})")
    ax.axhline(C_SAT_MOL_M3, ls="--", color="gray", label=f"C_sat = {C_SAT_MOL_M3:.2f} mol/m³")
    ax.plot(t_ref, C_t(t_ref, C_SAT_MOL_M3, KLA_BENCHMARK_S), ":", color="k",
            label="M-Star: k_La=4.1/hr, C_inf=C_sat")
    ax.set_ylabel("<C_L>  (mol/m³)")
    ax.set_title("(B) volume-averaged dissolved O₂ from bubble mass balance (upper bound)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")

    if have_direct_CL:
        ax = axes[2]
        ax.plot(t_C, C_C, ".-", ms=2, color="C2",
                label="<C_L> from m_component_lattices[0] (direct)")
        if math.isfinite(kLa_dir):
            t_fit = np.linspace(0, t_C[-1], 400)
            ax.plot(t_fit, C_t(t_fit, C_inf_dir, kLa_dir), "-",
                    label=f"fit: k_La = {kLa_dir_hr:.2f} /hr (C∞ = {C_inf_dir:.2f})")
        ax.axhline(C_SAT_MOL_M3, ls="--", color="gray", label=f"C_sat = {C_SAT_MOL_M3:.2f} mol/m³")
        ax.plot(t_ref, C_t(t_ref, C_SAT_MOL_M3, KLA_BENCHMARK_S), ":", color="k",
                label="M-Star: k_La=4.1/hr, C_inf=C_sat")
        ax.set_ylabel("<C_L>  (mol/m³)")
        ax.set_title("(C) volume-averaged dissolved O₂ from solver field reduction (preferred)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")

    axes[-1].set_xlabel("physical time (s)")
    fig.suptitle(
        f"kla_bioreactor: 400 RPM, 0.4 L/min, target k_La ≈ 4.1 /hr  (M-Star benchmark)",
        fontsize=10,
    )
    fig.tight_layout()
    out = Path(__file__).parent / "kla_verification.png"
    fig.savefig(out, dpi=140)
    print()
    print(f"saved → {out}")


if __name__ == "__main__":
    main()
