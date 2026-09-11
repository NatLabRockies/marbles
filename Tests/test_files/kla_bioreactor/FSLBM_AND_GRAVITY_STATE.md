# Current state — FSLBM + gravity forcing

This document describes the current implementation of the free-surface LBM
(FSLBM) and the macroscopic body-force (gravity + Lagrangian-bubble
back-coupling) in `Source/LBM.cpp`.  Line numbers refer to the current state
of the file at the time of writing — they may drift; use the function names
as the stable anchors.

It is a living document — please update it when the code below changes.

---

## 1. Step ordering (one LBM step on a level)

For each level `lev`, `LBM::advance_one_step` (≈ Source/LBM.cpp:800) executes:

```
 1. Moving-body sync           (reconstruct SDF, refill_and_spill)
 2. FSLBM advance              (fslbm_advance_surface)  ← replaces stream(m_f)
 3. Stream components          (stream(m_component_lattices[i]))
 4. Stream m_g                 (stream(m_g))
 5. Replenish m_g at interface (fslbm_replenish_g)      ← Change B (this rev.)
 6. Collide                    (f_to_macrodata → equilibrium → relax)
 7. Catalyst / reactions       (optional)
 8. Bubbles + macroscopic force (apply_macroscopic_forcing)  ← Change A
```

`apply_macroscopic_forcing` runs *after* `collide`.  This is the standard
exact-difference forcing order: relax first, then post-collision shift the
equilibrium by the velocity increment Δu = F·dt/ρ.  The next step's
`f_to_macrodata` automatically picks up the new velocity from the
post-shift `f`.

---

## 2. FSLBM — free-surface algorithm

### Cell taxonomy (Körner 2005)

For each cell:

- `CELL_SOLID`     — `IS_FLUID = 0`  (EB / impeller / baffle interior)
- `CELL_LIQUID`    — bulk fluid, φ = 1
- `CELL_GAS`       — empty headspace, φ = 0, populations zeroed
- `CELL_INTERFACE` — one-cell band straddling the surface,
                     φ ∈ (PHI_LO, PHI_HI) ≈ (1e-4, 1−1e-4)

`m_phi_fslbm` stores the fill fraction φ.  `m_cell_type` stores the integer
category.

### Per-step procedure (`fslbm_advance_surface`, ≈ Source/LBM.cpp:6383)

#### Body-motion sync (Pass 1 + Pass 2)

Reconcile cell types with the IS_FLUID mask written by the moving-body
refill:

- `LIQUID/INTERFACE` swept INTO by the body → `SOLID`, φ = 0.
- `SOLID` swept OUT of by the body → `LIQUID`, φ = 1 *unconditionally*.
  (The impeller is far below the free surface; there is no physical
  scenario where a body-vacated cell should be GAS.  An earlier
  `avg_phi`-heuristic produced spurious gas pockets and was removed.)

#### Step 0 — repair

Any `LIQUID/INTERFACE` cell with `rho < 0.01·rho_ref` is reseeded:

- `f` ← f_eq(ρ_ref, u=0, T_ref) using `pdiag_ref = (R/m_bar)·T_ref`
  (NOT `mesh_speed²`, which would seed at T = γ·T_ref).
- `g` ← g_eq(2ρe_ref, u=0, T_ref) via `set_extended_grad_expansion_generic`.

#### Step 0b — overshoot clamp

- `INTERFACE`: clamp `ρ > 5·ρ_ref` (ABB amplification risk at the surface).
- `LIQUID`: clamp NaN/Inf only; legitimate spill deposits are left to
  dissipate by streaming over a few steps.  (Clamping liquid at ρ_ref
  destroyed mass and broke the bulk flow.)

#### Step 1a — push streaming with bounce-back at solids

For each `LIQUID/INTERFACE` cell `iv`, push fq to neighbour `iv + e_q`:

- fluid neighbour → normal push.
- solid neighbour → bounce back (slot `f_star[iv][bq]`).
- gas neighbour   → leave the gas-facing slot zero
                    (filled in Step 1b in a separate pass to avoid the GPU
                    write race with Step 1b's ABB).

#### Pre-1b — interface normal & curvature

```
n̂ = ∇φ / ‖∇φ‖
κ = −∇·n̂      (κ > 0 ⇒ centre of curvature in the gas)
```

Wall correction with contact angle θ_W: at solid neighbours,

```
φ_ghost = φ + cos(θ_W) · ‖∇_∥ φ‖
```

θ = 90° (default) ⇒ cos = 0 ⇒ neutral wetting (original Körner).

#### Step 1b — Anti-Bounce-Back at gas neighbours (Source/LBM.cpp:6814–6885)

For each `INTERFACE` cell, for each direction q whose source is `GAS`:

```
f_star_q = f_eq_q(ρ_G, u, T) + f_eq_bq(ρ_G, u, T) − f_bq(x)
```

with gas-side density

```
ρ_G = ρ_ref + Δρ_Laplace
Δρ_Laplace = − 2σκ / (R_g T)
```

We use **ρ_ref**, NOT the local ρ.  Using local ρ creates a positive
feedback loop (elevated ρ → ABB targets elevated ρ → blow-up).
The result is `max(0, …)` clamped — ABB can produce f < 0 for large
non-equilibrium stress.

#### Step 2 — mass flux (pull view, pre-stream `f`)

```
Δm_iv = Σ_q  S_q · ( f_bq(iv + e_q)  −  f_q(iv) )
```

with

```
S_q = 1                         if neighbour is LIQUID
S_q = 0.5·(φ_iv + φ_neighbour)  if neighbour is INTERFACE
S_q = 0                         otherwise (GAS or SOLID)
```

#### Step 3 — commit

`m_f ← f_star`, `FillBoundary`.

#### Step 4 — phi update

```
φ^{n+1} = φ^n + Δm / ρ_post
```

`ρ_post` is computed directly from the *post-stream* `m_f`, NOT from
`m_macrodata` (which is one step stale and would be ≈ 0 for cells that
were just repaired / re-activated, sending φ to ∞).

#### Step 5 — cell conversion

- φ < `FSLBM_PHI_LO` ≈ 1e-4 → `GAS`; zero `f` and `g`; record
  excess = φ·ρ (≤ 0).
- φ > `FSLBM_PHI_HI` ≈ 1 − 1e-4 → `LIQUID`; record
  excess = (φ − 1)·ρ (≥ 0).

Component lattices in cells that just became `GAS` are **also zeroed**
in a follow-up kernel (Source/LBM.cpp:7099–7117) so dissolved scalars
don't persist in gas cells.

#### Step 5a — excess-mass redistribution (Pohl 2008 / Schwarzmeier 2023)

Convert-cell excess mass is distributed to neighbouring `INTERFACE`
cells weighted by `n̂ · ê_i` so that global mass is conserved.

#### Step 5b — interface spawning

Triggered only by conversions in Step 5 (Körner 2005 §3.3):

- A `LIQUID` cell with a neighbour that just became `GAS` is demoted
  to `INTERFACE` with φ = `FSLBM_PHI_HI` (PDFs already valid).
- A `GAS` cell with a neighbour that just became `LIQUID` is promoted
  to `INTERFACE` with φ = `FSLBM_PHI_LO`.  Its `f` and `g` are seeded
  from a fluid neighbour, normalised to `ρ_ref` so the cell enters the
  next step's streaming at the correct bulk density.  A spawn flag
  `flag[iv][2] = 1` is recorded.

#### Step 5c — component-lattice spawning

For every cell flagged in Step 5b, all component populations are
zeroed (Source/LBM.cpp:7267–7290).  The next step's `stream()` will
mass-conservingly refill them from neighbours; cloning a neighbour's
populations would duplicate mass because Eulerian components have no
volume-fraction weighting.

#### Step 6 — sync `IS_FLUID` and friends (`fslbm_sync_isfluid_markers`)

Derive `IS_FLUID_IDX` from updated `m_cell_type`, then recompute:

- `EB_BOUNDARY_IDX`
- `IS_FLUID_SIDE_IDX` (only adjacent to *moving* solid via `stat_mask == 1`,
  so impeller cells get the no-slip BC but baffle cells do not, and gas
  cells never get the impeller velocity)
- `IS_FLUID_SIDE_BOUNDARY_IDX`

---

## 3. Energy population `m_g` at the interface — Change B

After Step 4 above, the standard `stream(m_g)` runs on the
streamed-but-unconverted layout.  `INTERFACE` cells receive zero from
their gas-facing neighbours (gas cells carry no `g`).

### Closure (`fslbm_replenish_g`, Source/LBM.cpp:6349–6379) — symmetric bounce-back

```
for each CELL_INTERFACE iv:
    for each q with cell(iv + e_q) == CELL_GAS:
        g(iv, bounce_dirs[q]) ← g(iv, q)
```

### Why bounce-back

Donath (2011) gives no analytical g-closure — the dissertation only
treats hydrodynamic mass conservation.  This bounce-back closure is the
simplest physically meaningful condition:

- Σ_q g_q e_q · n̂ ≈ 0 in the gas-facing directions ⇒ approximately
  zero heat flux through the surface (adiabatic).
- Energy-conservative — no spurious thermal source on `INTERFACE` cells.
- Removes the previous T_ref Dirichlet injection that was depositing
  energy on every newly activated cell each step.
- Identical in form to `fslbm_replenish_components` — same logic that
  the surrounding code already uses for the component lattices.

---

## 4. Component lattices in the FSLBM pipeline

Components (dissolved scalars: O₂, products, etc.) are handled in five
places per step.  Important fact: `stream()` (Source/LBM.cpp:1076) keys
on **`IS_FLUID_IDX`**, which is set to 0 for *both* SOLID *and* GAS cells
by `fslbm_sync_isfluid_markers`.  So the standard streaming kernel's
bounce-back branch already provides a symmetric-bounce-back closure at
gas-facing slots, just like at solid walls.  The interface closure for
components exists — it just isn't in a function named after it.

| #   | Where                                          | What happens                                                                                              |
| --- | ---------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| 1   | refill_and_spill, Source/LBM.cpp:3049 (spill)  | Cells swept INTO by the body have their components zeroed and spilled to fluid neighbours.                |
| 1   | refill_and_spill, Source/LBM.cpp:3530 (refill) | Cells swept OUT of (newly fluid) get their `q=0` slot donor-filled from a normal-direction neighbour.     |
| 2   | stream(),         Source/LBM.cpp:840           | Push streaming + bounce-back at any non-fluid neighbour (gas OR solid) via `IS_FLUID == 0`.               |
| 3   | clamp,            Source/LBM.cpp:1125          | Optional `lbm.clamp_component_densities = 1`: zero all q if `Σ_q f_q < 0`.                                |
| 4   | collide,          Source/LBM.cpp:1618          | BGK relaxation toward `f^eq_q = ρ_comp · f̂^eq_q(u,T)` (unit-density shape scaled by local component ρ).   |
| 5a  | FSLBM Step 5,     Source/LBM.cpp:7100          | INTERFACE → GAS conversion: components zeroed in those cells (alongside f, g).                            |
| 5b  | FSLBM Step 5c,    Source/LBM.cpp:7268          | GAS → INTERFACE spawn: components in the spawned cell zeroed; next step's `stream()` mass-conservingly refills. |

### About `fslbm_replenish_components` (dead code)

A function `fslbm_replenish_components` is **defined**
(Source/LBM.cpp:6326) but **never invoked**.  It would *replace* the
bounce-back slot's value with the post-stream incoming-from-interior
population:

```
for each INTERFACE iv, for each q with cell(iv + e_q) == GAS:
    c[iv][ bounce_dirs[q] ] ← c[iv][ q ]      // post-stream values
```

That is **not the same** as the bounce-back that `stream()` already
applied (which uses pre-stream values).  Both give ≈ zero normal flux
on average, and the practical difference for dissolved-O₂ transport
is small.  Leaving it dead is defensible.

### Net assessment

Components are properly closed at the free surface, properly cleaned
at gas conversion, and properly seeded at interface spawning.  Nothing
is structurally missing.  `fslbm_replenish_components` is a legacy
alternative closure — keep it as documented dead code or delete it,
not a "fix that needs wiring in".

---

## 5. Gravity body force — Change A

Function: `apply_macroscopic_forcing(int lev, const amrex::MultiFab* force_mf)`
(Source/LBM.cpp:5819).

### 5.1 Units conversion

```
g_LB = g_phys · dt_phys² / dx_phys
```

Inputs from input file (defaults are 1.0, so unset ⇒ "force is already in
LB units"):

- `lbm.gravity   = gx gy gz` [m/s²]
- `lbm.dx_phys   = …`        [m]
- `lbm.dt_phys   = …`        [s]
- `lbm.dt_lev    = …`        backward-compat alias for `dt_phys`

For the kLa case:

```
g_LB,z = −9.81 · (1.19365e-5)² / 1.0e-3 ≈ −1.398 × 10⁻⁶  per step
```

(per-step velocity decrement in LB units).

### 5.2 Where it applies (Source/LBM.cpp:5870–5878)

- Skip cells with `IS_FLUID != 1` (solids, gas, masked).
- Skip `CELL_INTERFACE` (Donath 2011 p.122 — applying body force across
  the FSLBM ABB triggers acoustic shocks; the hydrostatic boundary is
  already enforced by Step 1b).
- Apply on `CELL_LIQUID` only — the bulk that needs the restoring buoyancy.

Body-force density per cell:

```
F = ρ · g_LB  +  ρ · clip(a_bubble, ±50·|g_LB|)
```

The bubble field is interpreted as an *acceleration* (LB velocity shift
per step) already normalised to the pure-liquid reference density;
multiplying by the local ρ converts it to a body-force density so the
1/ρ in the velocity shift below cancels.

### 5.3 Velocity shift

```
Δu = F · dt / ρ
```

### 5.4 Equilibrium update on `f` — exact-difference (Source/LBM.cpp:5985–5990)

```
f_q  +=  f_eq_q(ρ, u + Δu, Π_1)  −  f_eq_q(ρ, u, Π_0)
```

The diagonal stress entries include the *local* temperature
(NOT 1/3):

```
Π_aa = u_a² + r_temperature + dt · ω_corr · D_corr_a
r_temperature = R_g · T = P / ρ
ω_corr        = (2 − ω) / (2 ω ρ)
```

with `R_g = R_u / m_bar`, ω from the BGK relaxation, and the SGS
correction `D_corr` pulled from `m_derived`.  The mesh speed
(`m_dx_outer / m_dt_outer`) is the dimensional anchor of the equilibrium;
the equilibrium itself is a **product-form** equilibrium that uses
`R_g · T` as the diagonal stress entry — never 1/3.

### 5.5 Equilibrium update on `g` — thermodynamic consistency (Source/LBM.cpp:6005–6011)

The energy moment is

```
2ρe = ρ · (2 c_v T + |u|²)
```

so when u changes by Δu (with T held fixed), the energy moment also
changes by ρ · (|u + Δu|² − |u|²).  The Grad-expansion equilibrium
g_eq_q takes the heat-flux vector

```
q = 2 ρ u h
```

and the R-tensor

```
R_ab = 2 ρ u_a u_b (h + P/ρ)  +  2 P h δ_ab
h    = e + P/ρ        (functions of T and |u|²)
```

We compute both states (subscripts 0 and 1) and apply:

```
g_q  +=  g_eq_q( 2ρe_1, q_1, R_1 )  −  g_eq_q( 2ρe_0, q_0, R_0 )
```

Skipping this would inject a spurious internal-energy source
proportional to F · u, drifting T over time.

The argument `THETA0 = 1/3` passed to
`set_extended_grad_expansion_generic` is the **lattice temperature** of
the D3Q27 stencil (a property of the weights / abscissae used in the
Grad expansion) — it is NOT a flow speed of sound and is unrelated to
c_s².

### 5.6 Speed of sound — preserved everywhere

```
c_s = √( γ · R_g · T )      cell-local, never replaced by 1/√3.
```

The macroscopic equilibrium uses `R_g · T` as the diagonal-stress
entry; the acoustic c_s² = γ · R_g · T emerges from the moments of
that equilibrium.  Only `THETA0 = 1/3` appears in the code, exclusively
as the lattice temperature in the Grad expansion (mathematical, not
physical).

---

## 6. Files touched

- `Source/LBM.H`  — added `m_dx_phys`, `m_dt_phys`, `m_gravity`;
  renamed declaration `apply_bubble_body_force(…&)` →
  `apply_macroscopic_forcing(int, const amrex::MultiFab*)`.
- `Source/LBM.cpp`
  - parser block (extends `lbm.*`: `dx_phys`, `dt_phys`, `dt_lev`,
    `gravity`);
  - `apply_macroscopic_forcing` body (full rewrite, ~200 lines);
  - `advance_one_step` call site + gravity-only `else` branch;
  - `fslbm_replenish_g` rewritten to symmetric bounce-back;
  - pre-existing `m_bubble_params.dt_lev` → `dt_phys` bug fixed in 3 places.
- `Tests/test_files/kla_bioreactor/kla_bioreactor.inp` — added
  `lbm.gravity = 0.0 0.0 -9.81`.

---

## 7. Open questions / follow-ups

1. **`fslbm_replenish_components` is dead code** — see §4.  Action:
   either delete it, or document explicitly that the standard `stream()`
   already provides the analogous closure via the IS_FLUID mask.

2. **Smoke-test the kLa case** with the new gravity entry — confirm:
   - the gravity log line appears at startup;
   - no early divergence (first ~2000 steps);
   - the free surface starts settling toward hydrostatic equilibrium
     instead of producing floating clumps.

3. **Refined-level gravity branch.**  When bubbles are enabled, only
   `lev == 0` enters the bubble block, so `apply_macroscopic_forcing`
   on refined levels currently goes through the gravity-only `else`
   branch.  This is correct, but worth double-checking that the
   branch fires for every level when AMR is on.

---

*Last updated: 2026-06-11 — first draft after the GPU-port stabilisation
session that restored gravity forcing and replaced the m_g Dirichlet
closure with symmetric bounce-back.*
