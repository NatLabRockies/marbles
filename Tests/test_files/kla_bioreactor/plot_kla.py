#!/usr/bin/env python3
"""
k_La verification for the kLa-bioreactor run against the Reference Solver benchmark
(reference solver documentation)
target k_La ≈ 4.1 /hr at 400 RPM, 0.4 L/min air, dilute aqueous O₂.

Focus: Method C only — direct volume-averaged ⟨C_L⟩(t) reduced from the
solver's component-0 lattice over CELL_LIQUID + φ·CELL_INTERFACE.  This is
written every bubble.stats_int steps into bubble_stats.csv as the
C_L_mol_m3 column (along with V_liq_m3).

Three estimators of k_La, in order of physical defensibility:

  (1) Constrained exponential fit with C∞ ≡ C_sat = 1.428 mol/m³.
      Matches the Reference Solver "gassing-out" measurement convention exactly:
          C_L(t) = C_sat · (1 − exp(−k_La · t))
      Single-parameter fit, robust even when the trajectory is
      mostly in the linear-rise regime (curvature near saturation
      not yet visible).

  (2) Log-linear regression on  ln(C_sat − C_L)  vs  t.
      Same model, but solved via a linear least squares so the
      uncertainty is dominated by point scatter rather than the
      optimiser bracket.  Slope = −k_La directly.

  (3) Local slope at the trajectory's tail:
          k_La,local ≈ (dC_L/dt) / (C_sat − C_L)
      Most direct — no model assumption beyond C∞ = C_sat — but
      noisy at short times.  Reported only for sanity.

When the trajectory is in the very low-C regime (⟨C_L⟩ ≪ C_sat), free
fits are degenerate: every (C∞, k_La) pair satisfying k_La·C∞ = (slope)
gives the same residual.  This is why previous unconstrained fits
returned C∞ ≈ 14 mol/m³ (10× C_sat — unphysical).  Imposing
C∞ = C_sat picks the right branch.
"""

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# --------------------------------------------------------------------------
# 1. CONFIG — keep in sync with kla_bioreactor.inp
# --------------------------------------------------------------------------
DT_PHYS          = 1.19365e-5     # s / LB step
C_REF            = 100.0          # mol/m³ per LB-rho unit (lbm.O2_concentration_reference)
C_SAT_MOL_M3     = 0.032 * 44.6   # = 1.428 mol/m³  Henry-saturated liquid
KLA_BENCHMARK_S  = 4.1 / 3600.0   # 1.139e-3 /s (Reference Solver reference)

SCRIPT_DIR = Path(__file__).resolve().parent
WORK_CSV   = SCRIPT_DIR / "bubble_stats.csv"


# --------------------------------------------------------------------------
# 2. Parse bubble_stats.csv (Method C source)
# --------------------------------------------------------------------------
def parse_bubble_stats(path):
    """Return dict-of-arrays keyed by column name.

    Handles the 8-col (legacy) and 10-col schemas.  Also handles restart-
    induced duplicate rows: the bubble stats writer opens the CSV with
    ``std::ios::app`` on restart, so when the run resumes from a
    checkpoint whose step index precedes the previous end-of-file step,
    the overlapping window appears twice.  We deduplicate by the ``step``
    column keeping the LAST occurrence for each step (which is the row
    written by the most recent restart, i.e. the freshest physics).

    Post-blowup truncation: if a run enters a catastrophic FSLBM state
    (rare in DOUBLE, more common in FLOAT — the LBM component-lattice
    can hit a NaN cascade at the bottom-wall boundary after O(1e6)
    steps), the tail of the CSV holds sentinel-value junk (n_bubbles=0,
    d_min=1e30, C_L~1e-21).  Downstream fits crash on these values.
    Truncate at the first row where the bubble population has fully
    collapsed OR C_L drops below zero or above 2*C_sat — all valid
    physics rows are preserved intact.
    """
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        for row in reader:
            parsed = {}
            for k, v in row.items():
                if k is None or v is None or v == "":
                    continue
                try:
                    parsed[k] = float(v)
                except ValueError:
                    continue
            if parsed:
                rows.append(parsed)

    # Dedupe by step (last-write-wins), then sort ascending.
    if rows and "step" in rows[0]:
        by_step = {}
        for r in rows:
            by_step[r["step"]] = r
        rows = [by_step[s] for s in sorted(by_step)]

    # Truncate at first post-blowup row.  A row is "bad" if:
    #   * n_bubbles == 0 AND the previous row had n_bubbles > 100
    #     (population totally collapsed — physical run cannot bring it
    #     back to hundreds within one stats interval)
    #   * OR C_L_mol_m3 outside a physically-plausible band
    #     [-C_SAT*0.01, 2*C_SAT], catching Inf / negative-1e30 sentinels
    if rows:
        n_col = "n_bubbles"
        c_col = "C_L_mol_m3"
        prev_n = None
        cutoff = None
        for i, r in enumerate(rows):
            n = r.get(n_col)
            c = r.get(c_col)
            bad = False
            if (n is not None and n == 0 and prev_n is not None
                    and prev_n > 100):
                bad = True
            if c is not None and (c < -0.01 * C_SAT_MOL_M3
                                  or c > 2.0 * C_SAT_MOL_M3
                                  or not np.isfinite(c)):
                bad = True
            if bad:
                cutoff = i
                break
            prev_n = n
        if cutoff is not None:
            rows = rows[:cutoff]

    cols = {}
    for r in rows:
        for k, v in r.items():
            cols.setdefault(k, []).append(v)
    return {k: np.array(v) for k, v in cols.items()}


# --------------------------------------------------------------------------
# 3. Estimator (1): constrained exponential C∞ = C_sat
# --------------------------------------------------------------------------
def C_t_constrained(t, kLa):
    """C(t) = C_sat · (1 − exp(−kLa·t)) with C∞ fixed at saturation."""
    return C_SAT_MOL_M3 * (1.0 - np.exp(-kLa * t))


def fit_constrained(t, C, fit_start_s=0.2):
    """Single-parameter k_La fit with C∞ ≡ C_sat."""
    mask = t >= fit_start_s
    if mask.sum() < 5:
        raise RuntimeError("not enough points after fit_start_s")
    t_fit, C_fit = t[mask], C[mask]
    popt, pcov = curve_fit(C_t_constrained, t_fit, C_fit,
                           p0=[KLA_BENCHMARK_S],
                           bounds=([1e-6], [1.0]),
                           maxfev=20000)
    return popt[0], float(np.sqrt(pcov[0, 0])), t_fit, C_fit


# --------------------------------------------------------------------------
# 4. Estimator (2): log-linear regression on ln(C_sat − C_L) vs t
# --------------------------------------------------------------------------
def fit_loglinear(t, C, fit_start_s=0.2):
    """ln(C_sat − C_L) = ln(C_sat) − k_La · t.  Returns (k_La, sigma_kLa)."""
    mask = (t >= fit_start_s) & (C < C_SAT_MOL_M3) & (C > 0)
    if mask.sum() < 5:
        raise RuntimeError("not enough points for log-linear fit")
    t_fit = t[mask]
    y     = np.log(C_SAT_MOL_M3 - C[mask])
    # Linear fit y = a + b·t  ⇒  k_La = −b
    coeffs, cov = np.polyfit(t_fit, y, 1, cov=True)
    b, a = coeffs
    kLa_loglin = -b
    sigma_kLa  = float(np.sqrt(cov[0, 0]))
    return kLa_loglin, sigma_kLa, t_fit, y


# --------------------------------------------------------------------------
# 5. Estimator (3): tail-local slope
# --------------------------------------------------------------------------
def estimate_local(t, C, tail_frac=0.2):
    """k_La,local ≈ (dC/dt) / (C_sat − C) averaged over the last `tail_frac`."""
    n_tail = max(5, int(tail_frac * len(t)))
    t_tail, C_tail = t[-n_tail:], C[-n_tail:]
    slope, _ = np.polyfit(t_tail, C_tail, 1)
    C_bar = float(np.mean(C_tail))
    if C_bar >= C_SAT_MOL_M3:
        return float("nan"), C_bar, slope
    kLa_local = slope / (C_SAT_MOL_M3 - C_bar)
    return kLa_local, C_bar, slope


# --------------------------------------------------------------------------
# 6. Main
# --------------------------------------------------------------------------
def main():
    print("=== k_La verification (Method C only) ===")
    print(f"  benchmark target : k_La = 4.1 /hr = {KLA_BENCHMARK_S:.4e} /s")
    print(f"  C_sat (Henry)    : {C_SAT_MOL_M3:.3f} mol/m³")
    print(f"  CSV path         : {WORK_CSV}")
    print()

    b = parse_bubble_stats(WORK_CSV)
    if "C_L_mol_m3" not in b or not np.isfinite(b["C_L_mol_m3"]).any():
        raise SystemExit("bubble_stats.csv has no C_L_mol_m3 column — old schema?")

    t = b["phys_time_s"]
    C = b["C_L_mol_m3"]
    keep = np.isfinite(C) & np.isfinite(t)
    t, C = t[keep], C[keep]

    print(f"data: {len(t)} samples, t ∈ [{t[0]:.3f}, {t[-1]:.3f}] s")
    print(f"      ⟨C_L⟩(end) = {C[-1]:.4f} mol/m³  ({100*C[-1]/C_SAT_MOL_M3:.2f}% of C_sat)")
    print()

    # --- (1) Constrained exponential ---
    try:
        kLa_c, err_c, _, _ = fit_constrained(t, C)
        kLa_c_hr = kLa_c * 3600.0
        err_c_hr = err_c * 3600.0
        print(f"[1] Constrained exponential fit (C∞ ≡ C_sat = {C_SAT_MOL_M3:.3f} mol/m³)")
        print(f"    k_La = {kLa_c:.4e} ± {err_c:.1e} /s")
        print(f"         = {kLa_c_hr:.2f} ± {err_c_hr:.2f} /hr")
        print(f"         = {kLa_c_hr/4.1:.2f}× benchmark (4.1 /hr)")
    except Exception as ex:
        print(f"[1] fit failed: {ex}")
        kLa_c = kLa_c_hr = float("nan")

    print()

    # --- (2) Log-linear regression ---
    try:
        kLa_ll, err_ll, _, _ = fit_loglinear(t, C)
        kLa_ll_hr = kLa_ll * 3600.0
        err_ll_hr = err_ll * 3600.0
        print(f"[2] Log-linear regression on ln(C_sat − C_L) vs t")
        print(f"    k_La = {kLa_ll:.4e} ± {err_ll:.1e} /s")
        print(f"         = {kLa_ll_hr:.2f} ± {err_ll_hr:.2f} /hr")
        print(f"         = {kLa_ll_hr/4.1:.2f}× benchmark (4.1 /hr)")
    except Exception as ex:
        print(f"[2] fit failed: {ex}")
        kLa_ll = kLa_ll_hr = float("nan")

    print()

    # --- (3) Local slope at the tail ---
    kLa_loc, C_bar, dCdt = estimate_local(t, C)
    print(f"[3] Local slope at tail (last 20% of data)")
    print(f"    dC/dt = {dCdt:.4e} mol/m³/s, ⟨C_L⟩ = {C_bar:.4f} mol/m³")
    if math.isfinite(kLa_loc):
        print(f"    k_La,local = (dC/dt)/(C_sat − ⟨C_L⟩) = {kLa_loc:.4e} /s")
        print(f"               = {kLa_loc*3600:.2f} /hr  ({kLa_loc*3600/4.1:.2f}× benchmark)")

    print()

    # ---- Plot ----
    fig, axes = plt.subplots(2, 1, figsize=(9, 8))

    # Linear plot of ⟨C_L⟩(t).  Auto-scale to the data range — when ⟨C_L⟩ is
    # only a few % of C_sat, plotting on the full [0, C_sat] range hides the
    # rise and the three fits visually collapse onto the x-axis.  Show C_sat
    # only via the right-side "% of saturation" axis and a corner annotation.
    ax = axes[0]
    ax.plot(t, C, ".-", ms=3, color="C0", label="⟨C_L⟩ (direct solver reduction)")
    t_ref = np.linspace(max(1e-3, t[0]), t[-1], 400)
    if math.isfinite(kLa_c):
        ax.plot(t_ref, C_t_constrained(t_ref, kLa_c), "-",
                color="C2",
                label=f"[1] Constrained fit (C∞=C_sat):  k_La = {kLa_c_hr:.2f} /hr")
    if math.isfinite(kLa_ll):
        ax.plot(t_ref, C_t_constrained(t_ref, kLa_ll), "--",
                color="C3",
                label=f"[2] Log-linear:                 k_La = {kLa_ll_hr:.2f} /hr")
    ax.plot(t_ref, C_t_constrained(t_ref, KLA_BENCHMARK_S), ":",
            color="k", label=f"Reference Solver benchmark: k_La = 4.1 /hr")
    # Y-range: pad above the data so the fitted curves (which extrapolate
    # toward C_sat) and the local slope all stay visible without compressing
    # everything to a line.  Use 3× the data range as headroom.
    C_max_data = float(max(C.max(), C_t_constrained(t[-1], kLa_c)
                            if math.isfinite(kLa_c) else 0.0))
    ax.set_ylim(0.0, max(C_max_data * 3.0, 1.1 * C[-1]))
    ax.set_xlabel("physical time (s)")
    ax.set_ylabel("⟨C_L⟩  (mol/m³)")
    # Mirror axis on the right: % of C_sat
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim()[0] * 100.0 / C_SAT_MOL_M3,
                 ax.get_ylim()[1] * 100.0 / C_SAT_MOL_M3)
    ax2.set_ylabel("⟨C_L⟩ / C_sat  (%)")
    ax.set_title(
        f"Method C: liquid-volume-averaged dissolved O₂  (⟨C_L⟩(end) = "
        f"{C[-1]:.4f} mol/m³ = {100*C[-1]/C_SAT_MOL_M3:.2f}% of C_sat = {C_SAT_MOL_M3:.3f} mol/m³)"
    )
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")

    # Log-linear plot of ln(C_sat − C_L) vs t  — slope = −k_La
    ax = axes[1]
    valid = (C > 0) & (C < C_SAT_MOL_M3) & (t >= 0.2)
    if valid.any():
        y = np.log(C_SAT_MOL_M3 - C[valid])
        ax.plot(t[valid], y, ".", ms=3, color="C0", label="data (t ≥ 0.2 s)")
        if math.isfinite(kLa_ll):
            slope, intercept = np.polyfit(t[valid], y, 1)
            ax.plot(t_ref, intercept + slope * t_ref, "-", color="C3",
                    label=f"log-linear fit: −k_La = {slope:.4e} /s ⇒ k_La = {-slope*3600:.2f} /hr")
        ax.plot(t_ref, np.log(C_SAT_MOL_M3) - KLA_BENCHMARK_S * t_ref, ":",
                color="k", label="Reference Solver slope (k_La = 4.1 /hr)")
        ax.set_xlabel("physical time (s)")
        ax.set_ylabel("ln(C_sat − ⟨C_L⟩)")
        ax.set_title("Log-linear view — slope = −k_La directly (steeper = larger k_La)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower left")

    fig.suptitle("kla_bioreactor — 400 RPM, 0.4 L/min — Method C analysis",
                 fontsize=10)
    fig.tight_layout()
    out = SCRIPT_DIR / "kla_verification.png"
    fig.savefig(out, dpi=140)
    print(f"saved → {out}")


if __name__ == "__main__":
    main()
