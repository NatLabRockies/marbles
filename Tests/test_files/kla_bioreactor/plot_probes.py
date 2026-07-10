#!/usr/bin/env python3
"""
Plot the four fixed-location O2 probes from probes.csv.

Probes are the M-Star Fig. 22 analog: sample dissolved O2 at four fixed
LB-cell positions every probe.stats_int steps.  For a lumped-parameter
gassing-out fit, each probe should follow

    C_probe(t) = C_sat * (1 - exp(-k_La * t))

with the same k_La as the tank-average (⟨C_L⟩) provided the tank is well
mixed.  Non-uniform probes tell us about mixing time; matched probes tell
us about accurate k_La.

Positions (from the run header):
    p0 = (150, 90, -40)   between impeller tip and side wall (M-Star analog)
    p1 = ( 90, 90,  20)   mid-column, tank centerline
    p2 = ( 90, 90, -70)   near sparger, below impeller
    p3 = ( 90, 90,  60)   just below free surface (z_surface = 70)

The reference lines follow plot_kla_fig25.py conventions:
    * M-Star Fig. 22 k_La = 4.1 /hr  (headline benchmark)
    * Experiment (Zakrzewski 2020) = 5.1 /hr
    * Experiment (M-Star docs) = 4.5 /hr
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from plot_kla import C_SAT_MOL_M3  # noqa: E402  same conversion as bubble stats


CSV_PATH = SCRIPT_DIR / "probes.csv"
OUT_PNG  = SCRIPT_DIR / "kla_probes.png"

# Reference k_La values (same convention as plot_kla_fig25.py)
MSTAR_KLA_FIG22_HR    = 4.1
EXPERIMENTAL_KLA_HR   = 4.5   # M-Star docs, unattributed
ZAKRZEWSKI_KLA_HR     = 5.1   # Thomas et al. 2021 citation

# The current LES run uses C_eq = 0.275 mol/m^3 (see kla_bioreactor.inp
# bubble.surface_C_eq_mol_m3).  C_SAT_MOL_M3 = 1.427 corresponds to the
# pure-O2 headspace (older primary run).  Fit each probe against the
# saturation value that matches the actual run configuration.
C_SAT_CURRENT = 0.275


# --------------------------------------------------------------------------
# Load probes.csv and dedup step-repeats from any restart-induced duplicate
# rows.  Same restart-safety pattern as plot_kla.parse_bubble_stats.  We
# use csv.DictReader because the CSV has a comment line (starting with #)
# right after the header that numpy.genfromtxt handles inconsistently.
# --------------------------------------------------------------------------
def parse_probes(path: Path) -> dict:
    """Return dict-of-arrays keyed by column name, sorted by step
    (last-write-wins on duplicates)."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            parsed = {}
            for k, v in row.items():
                if k is None or v is None or v == "":
                    continue
                try:
                    parsed[k] = float(v)
                except ValueError:
                    # Skip comment lines (e.g. "# probe positions ...") --
                    # DictReader gives them to us as rows whose fields all
                    # fail float() conversion.
                    continue
            if parsed:
                rows.append(parsed)

    if not rows:
        raise RuntimeError(f"no data rows found in {path}")

    # Dedup on step (last-write-wins), then sort ascending.
    by_step = {r["step"]: r for r in rows if "step" in r}
    rows = [by_step[s] for s in sorted(by_step)]

    cols = {}
    for r in rows:
        for k, v in r.items():
            cols.setdefault(k, []).append(v)

    out = {"step":        np.asarray(cols["step"], dtype=np.int64),
           "phys_time_s": np.asarray(cols["phys_time_s"], dtype=float)}
    for i in range(4):
        key = f"C_L_probe{i}_mol_m3"
        out[f"probe{i}_mol_m3"] = np.asarray(cols[key], dtype=float)
    return out


def fit_kla(t: np.ndarray, C: np.ndarray, C_sat: float,
            t_start: float = 0.2) -> float | None:
    """Constrained exponential fit; returns k_La in 1/s or None if too few data."""
    m = np.isfinite(t) & np.isfinite(C) & (t >= t_start) & (C > 0) & (C < C_sat)
    if m.sum() < 20:
        return None
    model = lambda tt, k: C_sat * (1.0 - np.exp(-k * tt))
    try:
        popt, _ = curve_fit(model, t[m], C[m],
                            p0=[1.0e-3], bounds=([1e-6], [1.0]),
                            maxfev=20000)
        return float(popt[0])
    except Exception:
        return None


def main() -> int:
    if not CSV_PATH.exists():
        print(f"ERROR: {CSV_PATH} not found", file=sys.stderr)
        return 1

    df = parse_probes(CSV_PATH)
    t = df["phys_time_s"]

    probe_labels = [
        "p0  (150, 90, -40)   impeller wake",
        "p1  ( 90, 90,  20)   mid-column",
        "p2  ( 90, 90, -70)   near sparger",
        "p3  ( 90, 90,  60)   near free surface",
    ]
    probe_colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]
    probe_cols   = [f"probe{i}_mol_m3" for i in range(4)]

    # Fit each probe against C_sat_current (matches the run's surface BC).
    fits_s = []
    for col in probe_cols:
        C = df[col]
        k = fit_kla(t, C, C_sat=C_SAT_CURRENT)
        fits_s.append(k)

    print(f"probes.csv   : {CSV_PATH}")
    print(f"n_samples    : {len(t)}   t_range: 0 -> {t[-1]:.2f} s")
    print(f"C_sat        : {C_SAT_CURRENT:.4f} mol/m^3  (matches surface_C_eq)")
    print("Per-probe fits (constrained exp against C_sat above).")
    print("NOTE: per-probe k_La is a LOCAL rate of approach to saturation, not")
    print("      the tank-average k_La that appears in kLa_fig25.png.  In a")
    print("      well-mixed tank the two agree; in a mixing-limited flow they")
    print("      diverge (probes near sparger/impeller saturate faster than")
    print("      the bulk average).")
    for lbl, k in zip(probe_labels, fits_s):
        if k is None:
            print(f"  {lbl:50s}  (insufficient data)")
        else:
            print(f"  {lbl:50s}  k_La = {k*3600:.3f} /hr")
    print(f"\nReferences:")
    print(f"  M-Star Fig.22 (probe)          : {MSTAR_KLA_FIG22_HR:.2f} /hr")
    print(f"  Experiment (M-Star docs)       : {EXPERIMENTAL_KLA_HR:.2f} /hr")
    print(f"  Experiment (Zakrzewski 2020)   : {ZAKRZEWSKI_KLA_HR:.2f} /hr")

    # ------------------------------------------------------------------
    # Two-panel plot: raw trace + fits (left) and log-linear view (right).
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6), constrained_layout=True)
    tmax = float(t[-1])
    tref = np.linspace(0.0, tmax * 1.05, 400)

    # Left panel: raw C_L(t) at each probe + per-probe fit + reference lines.
    ax = axes[0]
    for lbl, col, colour, k in zip(probe_labels, probe_cols, probe_colors, fits_s):
        C = df[col]
        ax.plot(t, C, "-", color=colour, lw=1.4, alpha=0.85, label=lbl)
        if k is not None:
            ax.plot(tref, C_SAT_CURRENT * (1.0 - np.exp(-k * tref)),
                    "--", color=colour, lw=1.4, alpha=0.6,
                    label=f"    fit  k_La = {k*3600:.2f} /hr")
    ax.axhline(C_SAT_CURRENT, color="k", ls=":", lw=1.0,
               label=f"C_sat (surface_bc) = {C_SAT_CURRENT} mol/m³")
    # Reference exponentials at experimental k_La — for visual anchor.
    ax.plot(tref, C_SAT_CURRENT * (1.0 - np.exp(-MSTAR_KLA_FIG22_HR/3600 * tref)),
            "-", color="black", lw=2.0, alpha=0.6,
            label=f"M-Star Fig.22:  k_La = {MSTAR_KLA_FIG22_HR:.2f} /hr")
    ax.plot(tref, C_SAT_CURRENT * (1.0 - np.exp(-EXPERIMENTAL_KLA_HR/3600 * tref)),
            "--", color="#7f7f7f", lw=1.6, alpha=0.7,
            label=f"experiment (M-Star docs):  k_La = {EXPERIMENTAL_KLA_HR:.2f} /hr")
    ax.plot(tref, C_SAT_CURRENT * (1.0 - np.exp(-ZAKRZEWSKI_KLA_HR/3600 * tref)),
            "--", color="#8c564b", lw=1.6, alpha=0.7,
            label=f"experiment (Zakrzewski 2020):  k_La = {ZAKRZEWSKI_KLA_HR:.2f} /hr")
    ax.set_xlabel("Time  (s)")
    ax.set_ylabel(r"$C_L$ at probe  (mol/m³)")
    ax.set_title("Per-probe dissolved-O₂ trace + constrained-exp fits", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.95)
    ax.set_xlim(0.0, tmax * 1.02)
    ax.set_ylim(bottom=0.0)

    # Right panel: ln(C_sat - C) vs t — should be a straight line with
    # slope = -k_La for a well-mixed exponential approach.  Deviations
    # from linearity reveal mixing-time transients per probe.
    ax = axes[1]
    for lbl, col, colour, k in zip(probe_labels, probe_cols, probe_colors, fits_s):
        C = df[col]
        m = np.isfinite(t) & np.isfinite(C) & (t > 0.2) & (C > 0) & (C < C_SAT_CURRENT)
        if m.sum() > 5:
            y = np.log(C_SAT_CURRENT - C[m])
            ax.plot(t[m], y, "-", color=colour, lw=1.4, alpha=0.85, label=lbl)
    # Reference slope lines at experimental k_La values.
    tref_r = np.linspace(0.2, tmax * 1.02, 200)
    ax.plot(tref_r, np.log(C_SAT_CURRENT) - MSTAR_KLA_FIG22_HR/3600 * tref_r,
            "-", color="black", lw=2.0, alpha=0.6,
            label=f"slope = -{MSTAR_KLA_FIG22_HR:.2f}/hr  (M-Star Fig.22)")
    ax.plot(tref_r, np.log(C_SAT_CURRENT) - ZAKRZEWSKI_KLA_HR/3600 * tref_r,
            "--", color="#8c564b", lw=1.6, alpha=0.7,
            label=f"slope = -{ZAKRZEWSKI_KLA_HR:.2f}/hr  (Zakrzewski 2020)")
    ax.set_xlabel("Time  (s)")
    ax.set_ylabel(r"$\ln(C_\mathrm{sat} - C_L)$  at probe")
    ax.set_title("Log-linear view — slope = −k_La if lumped-parameter model holds",
                 fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=7.5, framealpha=0.95)
    ax.set_xlim(0.0, tmax * 1.02)

    fig.suptitle(f"Fixed-location O₂ probes  (probes.csv, t_max = {tmax:.2f} s)",
                 fontsize=11)
    fig.savefig(OUT_PNG, dpi=140)
    print(f"\nsaved -> {OUT_PNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
