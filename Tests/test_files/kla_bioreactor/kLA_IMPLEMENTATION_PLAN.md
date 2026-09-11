# Implementation Plan: Two-Phase kLa Mass Transfer Test Case

**Primary reference:** [Reference Solver CFD — Predicting Mass Transfer (kLa) in Gasified Bioreactors]((reference solver documentation))  
**Model details:** Thomas et al. (2021), *A mechanistic approach for predicting mass transfer in bioreactors*, CES 237, 116538. DOI: 10.1016/j.ces.2021.116538

**Status:** FSLBM free-surface fully implemented and stable. Latest commit `5d270e3` (branch `moving-body`): surface tension + contact angle. Lagrangian bubble / kLa modules not yet started.

---

## 1. Overview and Objective

Implement a sparged, stirred bioreactor test case that predicts the volumetric
oxygen mass transfer coefficient $k_L a$ from first principles. The simulation
couples:

1. A single-phase LBM fluid (already in codebase)
2. A **Lagrangian bubble phase** (new) — discrete bubbles tracked under
   Newton's laws with coalescence and breakup
3. A **dissolved oxygen scalar** (extends existing component lattice
   infrastructure) with a bubble-to-fluid source term
4. A free surface through which bubbles exit

Target validation metric: $k_L a \approx 4.1\ \text{hr}^{-1}$ matching the
Reference Solver benchmark at 400 RPM, 0.4 L/min gas flow.

---

## 2. Reference Case Parameters

| Quantity | Value |
|---|---|
| Tank diameter $T$ | 0.18 m |
| Impeller diameter $D_i$ | 0.058 m |
| Impeller speed $N$ | 400 RPM (tip speed 1.23 m/s) |
| Fluid density $\rho$ | 1000 kg/m³ |
| Kinematic viscosity $\nu$ | $10^{-6}$ m²/s |
| Fluid volume | 3.9 L (headspace ~0.02 m) |
| Lattice spacing $\Delta x$ | 1 mm → 180 cells across $T$ |
| Time step (Co = 0.1) | ~79 µs/step |
| Re (impeller) | 23 000 |
| Taylor length scale | 4.5 mm |
| Kolmogorov length scale | 90 µm |
| Sparger holes | 7 × 0.8 mm diameter, below impeller |
| Gas flow rate | 0.4 L/min |
| Initial bubble diameter | 1.1 mm |
| O₂ diffusion coefficient $D_{O_2}$ | $2 \times 10^{-9}$ m²/s (Sc ≈ 500) |
| O₂ solubility $S$ | 0.032 (dimensionless) |
| O₂ molar volume | 0.0224 m³/mol → $C_g = 0.0446$ mol/L |
| O₂ saturation concentration | $S \cdot C_g = 0.00144$ mol/L |
| Target $k_L a$ | ~4.1 hr⁻¹ (measured 4.5 hr⁻¹) |

---

## 3. Physics Modules Required

### 3.1  Lagrangian Bubble Module (new — largest work item)

Bubbles are discrete point particles tracked via Newton's second law
(Thomas et al. Eq. 9):

$$m_i \ddot{\mathbf{x}}_i = \mathbf{F}_{g} + \mathbf{F}_{a} + \mathbf{F}_{p} + \mathbf{F}_{D}$$

Forces to implement:

| Force | Expression | Notes |
|---|---|---|
| Net gravity + buoyancy | $(\rho_b - \rho_f) V_b \mathbf{g}$ | $\rho_b \ll \rho_f$ for gas |
| Drag | $\frac{1}{2} C_D \rho_f A_b |\mathbf{u}_{rel}| \mathbf{u}_{rel}$ | literature $C_D$ (Yang et al.) |
| Added mass | $C_{AM} \rho_f V_b \dot{\mathbf{u}}_{rel}$ | $C_{AM} = 0.5$ sphere |
| Pressure gradient | $-V_b \nabla p$ | from macrodata |

**Time integration**: velocity Verlet algorithm (Allen & Tildesley 2017), as used
in Thomas et al.

**Two-way coupling** (Thomas et al. Eq. 10): Only $\mathbf{F}_a$ (added mass) and
$\mathbf{F}_D$ (drag) are fed back as a body force onto the containing fluid
lattice cell — the pressure gradient force is **excluded** because it is already
incorporated in the Navier-Stokes momentum equation on the fluid side:

$$\mathbf{F}_{f,j} = -\sum_{i \in j} \left( \mathbf{F}_{i,a} + \mathbf{F}_{i,D} \right)$$

Fluid properties at the bubble location (velocity, pressure, $\epsilon$) obtained
by trilinear interpolation from the surrounding lattice nodes.

**Data structure**: GPU-resident array-of-structs — position, velocity, diameter,
O₂ moles. Must support dynamic add/remove (injection, exit, coalescence, breakup).
AMReX `ParticleContainer` is the natural choice.

### 3.2  Bubble Coalescence

Use critical-approach-Reynolds-number criterion (Boshenyatov 2012):

$$Re_a = \frac{U_a \, d_h}{\nu}, \quad d_h = \frac{2 d_1 d_2}{d_1 + d_2}$$

- $Re_a > 40$: coalesce — create new bubble conserving volume
- $Re_a \leq 40$: elastic bounce

Requires bubble neighbour distance calculation every $N_{coal}$ steps (suggested
$N_{coal} = 5$; start coalescence after 2 s physical to allow steady startup).

### 3.3  Bubble Breakup

Energy-dissipation-driven. Equilibrium diameter from Kolmogorov-surface-tension
balance (Hinze 1955; Thomas et al. Eq. 11):

$$D_e \approx \left(\frac{\sigma}{\rho}\right)^{3/5} \epsilon_i^{-2/5}$$

Breakup occurs when $D_e < D_{\text{bubble}}$ at the bubble's location.
Breakup conserves total bubble volume; two daughter bubbles produced per event.

**Daughter size distribution** — three models for sensitivity study
(Thomas et al. Fig. 10 shows predictions are insensitive across these at the
target operating condition):

1. **Triangle / 20-80 split** *(default, physically motivated)*: $f_v$ sampled
   from a triangle distribution with lower = 0, upper = 0.5, mode = 0.2
   (i.e. the most probable split is 20%/80% by volume — Sungkorn et al. 2012,
   Xing et al. 2015)
2. **Equal split**: $f_v = 0.5$
3. **Random**: $f_v = \text{rand}(0, 0.5)$

> Note: Thomas et al. also tested a U-shaped distribution but the triangle
> distribution with mode = 0.2 is the reported default. The model is
> insensitive to breakup kernel choice but sensitive to whether breakup is
> active at all (omitting breakup under-predicts $k_L a$ by ~2×).

### 3.3b  Bubble Size — Hydrostatic Pressure Correction

Although the effect is small for the 5 L vessel (fluid depth ≈ 0.12 m, pressure
change < 1%), Thomas et al. include Boyle's law diameter correction (Eq. 14)
for completeness:

$$\frac{d_i}{d_{i,\text{atm}}} = \left(1 + \frac{\rho g h}{P_{\text{atm}}}\right)^{-1/3}$$

where $h$ is bubble depth below the free surface. Apply at each timestep to
scale the registered bubble diameter before computing surface area / $k_L$.

---

### 3.4  Gas Sparger — Bubble Injection BC

7 injection points at fixed positions below the impeller. At each injection
event:

- Create a new bubble particle with $d = d_0 = 1.1$ mm, pure O₂ ($C_g^0 =
  0.0446$ mol/L)
- Injection rate: $\dot{n}_{inj} = Q_{gas} / V_{bubble,0}$ bubbles/s

Input parameters needed: `sparger_positions`, `sparger_flow_rate`,
`bubble_initial_diameter`, `bubble_O2_concentration`.

### 3.5  Free Surface Model

Thomas et al. model the top fluid surface using the **conservative phase-field
method** of Chiu & Lin (2011), *J. Comput. Phys.* 230, 185–204
(DOI: 10.1016/j.jcp.2010.09.021).

#### 3.5.1  Conservative Phase-Field Equation (Chiu & Lin)

A phase-field variable $\Phi \in [0, 1]$ (0 = gas headspace, 1 = liquid) marks
the gas-liquid interface. The governing equation in conservative form is:

$$\frac{\partial \Phi}{\partial t} + \nabla \cdot (\mathbf{u} \Phi) =
\gamma \nabla \cdot \left[ \nabla \Phi - \frac{\Phi(1-\Phi)}{\epsilon} \frac{\nabla \Phi}{|\nabla \Phi|} \right]$$

Key properties:
- Derived from the Allen-Cahn equation but written in conservative (divergence)
  form — guarantees exact mass conservation with conservative discretisation
- One-step method (simpler than two-step conservative level-set)
- Interface thickness controlled by parameter $\epsilon$ (set to $\sim 1.5 \Delta x$)
- $\gamma$ is an interface compression coefficient (set to $|\mathbf{u}|$)

**Boundary conditions**: no-flux for $\Phi$ at all walls; slip condition for
momentum at the top face (free surface).

**Mass-redistribution**: after each timestep, any undershoot/overshoot in $\Phi$
is clipped and the residual error uniformly redistributed among interface cells
(0.001 < $\Phi$ < 0.999) — ensures total liquid volume is exactly conserved.

#### 3.5.2  Coupling to Bubble Phase

- As sparger injects bubbles, the volume they displace causes the phase-field
  interface to rise (modelled implicitly through the momentum/continuity equations).
- When a bubble's centre-of-mass is above the interface ($\Phi < 0.5$ at bubble
  location), it is removed from the particle container; its remaining O₂ exits
  the system without transferring to the fluid.
- The dynamic free-surface height at any time is read as the isocontour
  $\Phi = 0.5$, integrated to give the current gas hold-up volume.

#### 3.5.3  Implementation Options

The Chiu & Lin model requires solving an additional PDE per timestep on the
full 3-D domain — non-trivial work. Two approaches:

| Option | Description | Cost |
|---|---|---|
| **A — Full phase-field** | Solve Chiu & Lin Eq. at each LBM step; dynamic free surface with proper mass conservation; gas hold-up tracked | +1–1.5 dev days; ~10–15% runtime overhead |
| **B — Z-threshold (simplified)** | Remove bubbles at fixed $z_{\text{top}}$; slip-wall fluid BC at top face; gas hold-up not tracked dynamically | 0.5 dev days; negligible overhead |

**Recommendation**: implement **Option A** to match the Thomas et al. model
faithfully, since the dynamic free surface affects the dissolved-O₂ saturation
boundary condition and the reported gas-hold-up check. Option B remains a
valid initial approximation if schedule is tight — the 5 L vessel has < 1%
gas hold-up and the free surface moves by only a few mm.

### 3.6  Local Energy Dissipation Rate Field

$k_L$ depends on local $\epsilon_i$ at each bubble. The LBM stress tensor
already provides the viscous dissipation. Full dissipation (resolved +
sub-grid) is (Thomas et al. Eq. 7):

$$\epsilon_j = \nu_{T,j} \, S_j^2$$

where $\nu_T = \nu + \nu_{\text{SGS}}$ is the total (molecular + eddy) viscosity
and $S_j = \sqrt{2 S_{ij} S_{ij}}$ is the resolved strain rate magnitude.

**The SGS closure is important** — Thomas et al. use Smagorinsky with $C_s = 0.1$
(calibrated from DNS, Yu et al. 2005) and validate that total dissipated power
$P_d = \int \epsilon_j \, dm$ matches shaft power $P_s = \tau \cdot 2\pi N$
to within ~5%. This energy balance is the primary sanity check that $\epsilon$
is correctly computed before enabling mass transfer:

$$\epsilon_{\text{SGS}} = (C_s \Delta x)^2 |S|^3, \quad C_s = 0.1$$

For the initial run the solver already has SGS via the standard LBM effective
viscosity (`nu_eff = nu + nu_sgs`). Store $\epsilon_j$ as a derived field in
`m_derived` per step; optionally time-average over a window for smoother
bubble interpolation.

### 3.7  Mass Transfer — Bubble-to-Fluid Source Term

Per-bubble transfer rate (Thomas et al. Eqs. 15–16):

$$\dot{n}_i = H_i \left( S_i C_{g,i} - C_{f,i} \right), \quad H_i = k_{L,i} A_i$$

Note the sign convention: $S_i$ multiplies only the gas-phase concentration
$C_{g,i}$, converting it to an equivalent dissolved-phase concentration, so the
**effective driving force** is $S_i C_{g,i} - C_{f,i}$ in units consistent
with the dissolved-phase concentration $C_f$.

Local mass transfer coefficient from Kawase et al. (1992) (cited as Thomas et al. Eq. 17),
with $C = 0.301$:

$$k_{L,i} = 0.301 \left(\epsilon_i \nu_i\right)^{1/4} Sc_i^{-1/2}, \quad Sc = \frac{\nu}{D_{O_2}}$$

Dimensional check: $[(\epsilon_i \nu_i)^{1/4}] = [(m^2 s^{-3} \cdot m^2 s^{-1})^{1/4}] = [m^4 s^{-4}]^{1/4} = m\,s^{-1}$ ✓

> Note: Thomas et al. write this as $C(\epsilon_i\nu_i)^a Sc_i^b$ and state $a=-\tfrac{1}{2}$,
> $b=\tfrac{1}{4}$, which appears inconsistent with both dimensional analysis and the
> Kawase et al. source. The formula above follows Kawase et al. (1992) directly and is
> dimensionally correct.

Source term deposited on fluid cell $j$ containing bubbles $i \in j$
(Thomas et al. Eqs. 20–21):

$$\dot{c}_j^{\text{src}} = \frac{1}{V_j} \sum_{i \in j} \dot{n}_i$$

This maps directly onto the existing `apply_reaction_source_terms()` pattern
— instead of a volumetric reaction, it is a bubble-localised source distributed
over the containing cell.

Bubble volume shrinks as O₂ transfers (Thomas et al. Eq. 18–19):

$$\dot{V}_i = -V_m \dot{n}_i, \quad V_m = \frac{R_u T}{P_g}$$

Using $V_m = 0.0224$ m³/mol at STP.

**Scalar transport**: apply a **Van Leer flux limiter** to the dissolved O₂
advection step to suppress numerical diffusion, as done in Thomas et al.
(the existing scalar lattice advection already uses the LBM streaming step,
but an additional Van Leer correction may be needed at high Sc).

### 3.8  k_La Calculation — Two Independent Methods

Thomas et al. use two methods as a cross-check for mass conservation (Eqs. 25–26):

1. **Probe-based** (Eq. 22): fit the time-evolution of $C_f(t)$ at a probe
   location to $dC_f/dt = k_L a \cdot (S C_b - C_f)$ to extract $k_L a$.

2. **Volume-averaged** (Eqs. 25–26): time-average the local field
   $k_L a_j = \frac{1}{\tau} \int_0^\tau k_{L,j}(t) A_{T,j}(t) / V_j \, dt$,
   then integrate over the vessel:
   $k_L a = \frac{1}{V_l} \iiint k_L a_j \, dV$.

Both methods should agree; discrepancy > 5% indicates a bug in species
conservation. Implement both in the post-processing notebook.

---

## 4. Changes to Existing Code

### 4.1  `Source/LBM.H` / `LBM.cpp`

| Change | Description |
|---|---|
| Add `BubbleContainer` member | AMReX `ParticleContainer`-derived class |
| `initialize_bubbles()` | Read sparger params, create initial empty container |
| `advance_bubbles()` | Integrate Newton's law, handle injection/exit |
| `compute_bubble_coalescence()` | Neighbour search + merge |
| `compute_bubble_breakup()` | Per-bubble $D_e$ check + split |
| `compute_bubble_fluid_coupling()` | Body force deposition to fluid cells |
| `apply_bubble_mass_transfer()` | Compute $k_L$, $\dot{n}_i$, update scalar + bubble diameter |
| `compute_energy_dissipation()` | Store $\epsilon$ field in `m_derived` |
| `write_bubble_stats()` | CSV output: mean diameter, total $k_L a$, bubble count  |

### 4.2  `Source/Constants.H`

Add physical constants: O₂ molar volume, solubility, surface tension for air-water.

### 4.3  `Source/Stencil.H` / `Stencil.cpp`

No changes expected.

### 4.4  New file: `Source/Bubble.H` / `Source/Bubble.cpp`

Encapsulate the `BubbleContainer` class, force calculation kernels, and
coalescence/breakup logic.

### 4.5  CMakeLists.txt

Add `Bubble.cpp` to the source list. Enable `AMReX_PARTICLES = ON` in
`CMake/set_amrex_options.cmake` (check if already on).

---

## 5. New Input File Parameters

```
# Two-phase bubble parameters
lbm.enable_bubbles        = 1

# Sparger
sparger.n_points          = 7              # number of injection holes
sparger.positions_file    = "sparger.csv"  # x y z of each hole
sparger.flow_rate         = 6.67e-6        # m³/s (= 0.4 L/min)
sparger.bubble_diameter   = 1.1e-3         # m

# Bubble physics
bubble.enable_coalescence      = 1
bubble.coalescence_Re_crit     = 40.0
bubble.coalescence_start_time  = 2.0       # s
bubble.coalescence_interval    = 5         # steps between neighbour checks
bubble.enable_breakup          = 1
bubble.breakup_model           = "triangle"  # "triangle" (default), "equal", "random"
bubble.surface_tension         = 0.072     # N/m (air-water)

# Dissolved oxygen scalar (component 0)
lbm.n_components               = 1
lbm.diffusivity_component_0    = <D_O2 in LB units>
bubble.O2_molar_volume         = 0.0224    # m³/mol
bubble.O2_solubility           = 0.032     # dimensionless
bubble.O2_initial_conc_bubble  = 44.6      # mol/m³ (= 0.0446 mol/L)
bubble.kL_coefficient          = 0.301     # C in k_L = C(εν)^(1/4) Sc^(-1/2)  [Thomas et al. Eq. 17]

# Free surface height (bubble exit threshold)
bubble.free_surface_z          = <fluid top z in LB cells>

# Output
amr.bubble_stats_int      = 200            # steps between bubble CSV writes
```

---

## 6. New Test Directory

`Tests/test_files/kla_bioreactor/`

Contents:
- `kla_bioreactor.inp` — input file (parameters above)
- `sparger.csv` — 7 sparger hole positions
- `pitched_impeller.stl` — reuse from `stirred_tank_reacting/`
- `baffled_wall.stl` — same
- Reference data: Reference Solver measured $k_L a = 4.5\ \text{hr}^{-1}$

---

## 7. New Post-Processing Notebook

`Tools/kla_plot.ipynb` — plots:
1. Time-evolution of dissolved O₂ at probe (fit to extract $k_L a$)
2. Bubble size distribution (histogram, vs. Reference Solver Fig. 21)
3. Total transfer rate $\dot{n}_T$ vs time (vs. Reference Solver Fig. 24)
4. Shaft power vs total dissipation (energy check, vs. Reference Solver Fig. 20)
5. Optional: spatial O₂ distribution (from plotfiles)

---

## 8. Validation Metrics

| Quantity | Target | Source |
|---|---|---|
| $k_L a$ at probe (exponential fit) | 4.1 hr⁻¹ | Reference Solver guide |
| $k_L a$ volume-averaged | ~4.1 hr⁻¹ | Thomas et al. (4.8 hr⁻¹ at similar conditions) |
| Probe vs. volume-average agreement | < 5% | mass conservation check (Thomas et al.) |
| Shaft power number $N_p$ | ~1.5 | Reference Solver guide |
| Shaft power = dissipated power $P_d$ | < 5% difference | Thomas et al. energy balance check |
| Mean bubble diameter above impeller | ~3 mm | Reference Solver guide |
| Gas hold-up | < 1% | Thomas et al. 5 L result |
| Power balance error | < 10% | energy conservation |

---

## 9. Implementation Sequence (suggested order)

1. **Enable `AMReX_PARTICLES`** and verify build — 0.5 day
2. **`Bubble.H/cpp` skeleton**: data structure, injection, advection (no forces yet) — 1 day
3. **Gravity + drag + added mass forces**; velocity Verlet integrator;
   verify terminal rise velocity of a single bubble — 0.5 day
4. **Two-way coupling**: deposit $\mathbf{F}_a + \mathbf{F}_D$ only back to
   fluid (not pressure gradient); verify momentum conservation — 0.5 day
5. **Trilinear interpolation** of fluid velocity, pressure, $\epsilon$ at bubble
   position — 0.5 day
6. **Energy dissipation field** `compute_energy_dissipation()`; validate
   $P_d = P_s$ to < 5% with single-phase impeller run — 0.5 day
7. **Dissolved O₂ scalar** with bubble source term and Van Leer limiter;
   Boyle's law diameter correction; single-bubble transfer test — 1 day
8. **Coalescence** (Re$_a > 40$ criterion, elastic bounce otherwise) — 1 day
9. **Breakup** (Hinze criterion; triangle daughter-size distribution,
   mode = 0.2) — 1 day
10. **Sparger injection BC** + free surface: implement Chiu & Lin
    conservative phase-field equation for gas-liquid interface tracking;
    validate with single-bubble rise reaching the surface and correct
    gas hold-up — 1–1.5 days
11. **Full bioreactor run** + notebook with both $k_L a$ methods (probe fit +
    volume-average); validation against Reference Solver targets — 1–2 days
12. **Sensitivity study**: breakup kernels (triangle / equal / random);
    confirm insensitivity consistent with Thomas et al. Fig. 10 — 0.5 day
13. **Optional**: CO₂ scrubbing (second scalar + bidirectional transfer) — 1 day

**Total estimate**: ~10–13 development days excluding HPC run time.

---

## 10. Key Open Questions / Decisions

| Question | Options | Recommendation |
|---|---|---|
| Particle container GPU/CPU? | CPU `ParticleContainer` vs. GPU-offloaded | Start CPU; move to GPU if bubble count > 10⁵ |
| SGS turbulence closure? | None (resolved only) vs. Smagorinsky $C_s = 0.1$ | **Include from the start** — Thomas et al. show it is required to correctly predict $\epsilon$ and hence $k_L$ |
| Free surface model? | **Full**: Chiu & Lin (2011) conservative phase-field PDE (Thomas et al. approach); **Simplified**: Z-threshold exit + slip-wall BC | ~~Full phase-field~~ **→ FSLBM (Körner 2005)** — see Section 11 below |

---

## 11. Free-Surface Model: FSLBM — As Implemented

**Decision (2026-04-27)**: Switched from `advance_phi()` Chiu & Lin PDE to
**Free-Surface LBM (FSLBM)** of Körner et al. (2005).  The Chiu & Lin
implementation was unstable at the triple-point (gas–liquid–solid contact line)
and could not be stabilized with Neumann BCs, γ = global |u|_max, or
interface-cell equilibrium reinitialization.

**Implementation status (2026-05-08)**: Fully implemented and stable.
Latest commit `5d270e3` on branch `moving-body`.
50-step test: Δm ≈ machine-epsilon, interface rho ∈ [0.9, 1.1].
18 000+ step test (kla_bioreactor): density and temperature stable.

Schwarzmeier et al. (2023, JCP 473, 111753) show FSLBM outperforms
conservative Allen–Cahn PFLBM (≡ Chiu & Lin LBM) in every gravity/inertia-
dominated test (rising bubble, Taylor bubble, drop impact) — exactly our
bioreactor regime.

---


### 11.1 Why FSLBM Wins for Bioreactors

| Test (Schwarzmeier 2023) | FSLBM | PFLBM (Chiu & Lin) |
|---|---|---|
| Gravity wave | ✅ converges | ✅ converges |
| Rising bubble (D≥16) | ✅ accurate | ❌ unstable at D<32 |
| Taylor bubble | ✅ Re agrees with exp. | ❌ unstable at D≤32 |
| Drop impact (crown) | ✅ captures ejected droplets | ❌ misses droplets |
| Memory per cell | ~1× (single PDF lattice) | ~2× (dual lattice) |
| Sensitivity to mobility M | None | High — no robust M |

**Root cause of PFLBM failure**: mobility M must be in a narrow window; below
it → instability; above it → non-physical film rupture.  No universal M exists.
FSLBM has no free parameter analogous to M.

---

### 11.2 Data Structures and Cell Types

#### New fields (as implemented in `Source/LBM.H`)

| Field | Type | Components | Ghost cells | Purpose |
|---|---|---|---|---|
| `m_cell_type[lev]` | `iMultiFab` | 1 | `m_f_nghost` | Integer cell classification |
| `m_phi_fslbm[lev]` | `MultiFab` | 1 | `m_f_nghost` | Fill level φ ∈ [0,1] |
| `m_free_surface_z` | `Real` | — | — | z-coordinate of initial free surface |
| `m_free_surface` | `bool` | — | — | Enable/disable FSLBM |
| `m_fslbm_sigma` | `Real` | — | — | Surface tension σ in LB units (default 0) |
| `m_fslbm_contact_angle_deg` | `Real` | — | — | Static contact angle θ in degrees (default 90) |

Input parameters:
```
lbm.fslbm_sigma         = <σ_phys * dt² / (ρ_phys * dx³)>   # LB surface tension
lbm.fslbm_contact_angle = 90.0                                # degrees; 90=neutral
```

Unit conversion: $\sigma_\text{LB} = \sigma_\text{phys} \cdot \Delta t^2 / (\rho_\text{phys} \cdot \Delta x^3)$.
For air–water ($\sigma=0.072$ N/m) at bioreactor scale ($\Delta t=1.19\times10^{-5}$ s, $\Delta x=10^{-3}$ m): $\sigma_\text{LB}\approx 1.026\times10^{-5}$.

#### Cell type encoding (`Source/Constants.H`)

```cpp
constexpr int CELL_SOLID     = 0;  // EB solid or impeller interior
constexpr int CELL_GAS       = 1;  // gas headspace — f=0, IS_FLUID=0
constexpr int CELL_INTERFACE = 2;  // φ ∈ (1e-4, 1−1e-4) — f valid, IS_FLUID=1
constexpr int CELL_LIQUID    = 3;  // φ = 1 — f valid, IS_FLUID=1

constexpr amrex::Real FSLBM_PHI_LO = 1.0e-4;     // INTERFACE→GAS threshold
constexpr amrex::Real FSLBM_PHI_HI = 1.0-1.0e-4; // INTERFACE→LIQUID threshold
```

| Type | φ | IS_FLUID_IDX | f present? | g present? |
|---|---|---|---|---|
| SOLID | 0 | 0 | No (zeroed) | No |
| GAS | 0 | 0 | No (zeroed) | No |
| INTERFACE | (PHI_LO, PHI_HI) | 1 | Yes | Yes |
| LIQUID | 1 | 1 | Yes | Yes |

---

### 11.3 Initialization: `fslbm_init_cell_type(lev)`

**Location**: `Source/LBM.cpp` line 5577

Sharp-interface initialization from `m_free_surface_z` (physical metres).
Uses a **two-cell-wide** interface band (±1·Δz) instead of ±0.5·Δz, so that
when the surface falls exactly on a cell face both adjacent cells become
INTERFACE with φ ∈ (PHI_LO, PHI_HI).  This prevents immediate conversion in
Step 5 on the first time step.

**Algorithm**:

```
For each cell at z_cell = plo[2] + (k + 0.5) * dz:
  if IS_FLUID_IDX == 0:                    → CELL_SOLID, φ = 0
  elif z_cell >= z_surf + dz:              → CELL_GAS,   φ = 0
  elif z_cell <= z_surf - dz:              → CELL_LIQUID, φ = 1
  else (|z_cell - z_surf| < dz):          → CELL_INTERFACE
    φ_lin = (z_surf - z_cell) / (2·dz) + 0.5    (linear ramp bottom→top: 0→1)
    φ = clamp(φ_lin, PHI_LO, PHI_HI)
```

After kernel: `m_phi_fslbm.FillBoundary()`, then call
`fslbm_sync_isfluid_markers(lev)` to derive IS_FLUID from cell type.

Key insight: IS_FLUID **must** be derived from `m_cell_type`, not from
`m_is_fluid_fraction` / SDF threshold (which uses a 0.5 threshold and would
reclassify FSLBM interface cells each step).

---

### 11.4 IS_FLUID Marker Synchronization: `fslbm_sync_isfluid_markers(lev)`

**Location**: `Source/LBM.cpp` line 5443

Derives all four `m_is_fluid` components from `m_cell_type`.  Called:
1. At the end of `fslbm_init_cell_type` (once, at startup)
2. At the start of `fslbm_advance_surface` (after body-motion sync)
3. At the end of `fslbm_advance_surface` (Step 6, after conversion)

**Step 1 — IS_FLUID_IDX** (over full ghost region, then `FillBoundary`):

```
IS_FLUID_IDX = (ct == CELL_LIQUID || ct == CELL_INTERFACE) ? 1 : 0
```

**Step 2 — EB_BOUNDARY_IDX**:

```
For cell iv:
  if all axis-aligned neighbors within nn=1 have IS_FLUID==0 → EB_BOUNDARY=0
  else if IS_FLUID(iv)==1 → EB_BOUNDARY=0
  else → EB_BOUNDARY=1
```

**Step 3 — IS_FLUID_SIDE_IDX**
(fluid cells adjacent to the *moving* impeller, used by `f_to_macrodata` velocity BC):

```
IS_FLUID_SIDE = 0  if IS_FLUID(iv)==0
else IS_FLUID_SIDE = 1  if any D3Q27 neighbor has ct==CELL_SOLID AND stationary_mask==1
                                                   (i.e., the moving impeller)
```

Key exclusions:
- **CELL_GAS** has IS_FLUID=0 from Step 1, so it can never receive IS_FLUID_SIDE=1.
  This prevents the impeller no-slip BC from being applied to gas headspace cells.
- **Baffle (stationary solid)**: `stationary_mask==0` → excluded from IS_FLUID_SIDE.

**Step 4 — IS_FLUID_SIDE_BOUNDARY_IDX**:

```
IS_FLUID_SIDE_BOUNDARY = 1 if IS_FLUID==1 AND IS_FLUID_SIDE==0
                            AND any D3Q27 neighbor has IS_FLUID_SIDE==1
```

---

### 11.5 Main Time Step: `fslbm_advance_surface(lev)`

**Location**: `Source/LBM.cpp` line 5657

Replaces `advance_phi(lev)` + `stream(lev, m_f)` in `advance()`.

The full `advance()` call order with FSLBM:

```
refill_and_spill(lev)          [if body_is_moving]
fslbm_advance_surface(lev)     [replaces advance_phi + stream(m_f)]
stream(lev, m_component_lattices[i])
stream(lev, m_g)
collide(lev)                    [f_to_macrodata → collide → relax]
```

#### Body-Motion Sync (pre-Step 0)

Runs immediately after `refill_and_spill`.  Reconciles `m_cell_type` (set by
`fslbm_sync_isfluid_markers` at the *previous* step end) against the IS_FLUID_IDX
just written by SDF reconstruction + `refill_and_spill`.

**Case A — body swept INTO a fluid cell** (`ct==LIQUID||INTERFACE`, `IS_FLUID==0`):

```
ct → CELL_SOLID
φ  → 0
```

Without this, Step 1 would still stream PDFs from a LIQUID cell that IS_FLUID=0,
and Step 5's threshold would re-activate it as INTERFACE next step — overriding
the body's IS_FLUID=0 each cycle.

**Case B — body swept OUT of a solid cell** (`ct==SOLID`, `IS_FLUID==1`):

Classify by height:

```
z_cell < z_surf → ct = CELL_LIQUID, φ = 1
z_cell ≥ z_surf → ct = CELL_GAS,   φ = 0
```

Then **seed both f AND g** to `feq(ρ=1, u=0, p_diag = mesh_speed²)`:

```cpp
const Real pdiag = mesh_speed * mesh_speed;   // cs² = T_ref
for (int q = 0; q < N_MICRO_STATES; ++q) {
    Real fval = set_extended_equilibrium_value(
        1.0, {0,0,0}, pdiag, pdiag, pdiag, mesh_speed, weights[q], evs[q]);
    f[iv][q] = fval;
    g[iv][q] = fval;   // MUST seed g; otherwise T=0 in collide → entropic blow-up
}
```

After body-motion sync: `m_cell_type.FillBoundary()`, then
`fslbm_sync_isfluid_markers(lev)`.

#### Step 0 — Repair Kernel

Finds any CELL_INTERFACE or CELL_LIQUID cell with Σf < 0.01 and seeds both
f AND g to `feq(ρ=1, u=0, pdiag=cs²)`.

**Why needed**: cells near the tank wall start as IS_FLUID=0 (SDF < 0.5) so
their f is zero-initialized.  `fslbm_sync_isfluid_markers` marks them IS_FLUID=1
without seeding f.  Within 4–5 steps these cells accumulate negative Δm → φ → −∞
and collapse the entire free surface.

After Step 0: `m_f.FillBoundary(); m_g.FillBoundary()`.

#### Step 1a — Push Streaming

Allocate `f_star` (MultiFab, N_MICRO_STATES comps, initialized to zero).

For each **LIQUID or INTERFACE** cell `iv`, for each direction `q`:

| Neighbor `ivn = iv + ev_q` | Action |
|---|---|
| CELL_LIQUID or CELL_INTERFACE | `f_star[ivn][q] = f_pre[iv][q]` (push) |
| CELL_GAS | leave `f_star[ivn][q] = 0` — Step 1b fills it via ABB |
| CELL_SOLID | bounce-back: `f_star[iv][bq] = f_pre[iv][q]` |

The gas slot is **intentionally left zero** in 1a to avoid a GPU write race:
both the push from a face-sharing interface neighbor AND the ABB fill of the
current cell want to write the same `f_star[iv][q]` slot.  Separating into
two passes eliminates the race.

#### Pre-Step-1b — Interface Normal n̂ and Curvature κ

Required for the Laplace pressure ABB (Step 1b) when `m_fslbm_sigma > 0`.
Allocated and computed every step; zero-cost when σ = 0 because Δρ = 0.

**Pass 1 — unit normal** $\hat{n} = \nabla\varphi/|\nabla\varphi|$ (`nhat_mf`, 3 comps, 1 ghost):

For every non-SOLID cell, central-difference gradient of φ with **contact-angle
wall correction** (Attar & Körner 2009 / Schwarzmeier 2023): for a SOLID neighbor
in direction $(\Delta i, \Delta j, \Delta k)$, substitute:

$$\varphi_\text{ghost} = \varphi_\text{fluid} + \cos(\theta_W)\,|\nabla_\perp\varphi|$$

where $\nabla_\perp$ is the gradient tangential to the wall-normal direction.
At $\theta_W = 90°$: $\cos = 0 \Rightarrow \varphi_\text{ghost} = \varphi_\text{fluid}$
— reproducing Körner neutral wetting.

Implementation (`Source/LBM.cpp`, inside `fslbm_advance_surface`):
```cpp
// 90°-base gradients (SOLID → phi of current cell)
auto pf = [&](ii,jj,kk) { return solid? phi : phi_arr[ii,jj,kk]; };
gx0 = 0.5*(pf(i+1,j,k) - pf(i-1,j,k));
gy0 = 0.5*(pf(i,j+1,k) - pf(i,j-1,k));
gz0 = 0.5*(pf(i,j,k+1) - pf(i,j,k-1));
// Contact-angle ghost for each solid-neighbor direction
auto phi_w = [&](ii,jj,kk) {
    if (!solid) return phi_arr[ii,jj,kk];
    if (di!=0) gt = sqrt(gy0²+gz0²);
    elif (dj!=0) gt = sqrt(gx0²+gz0²);
    else         gt = sqrt(gx0²+gy0²);
    return phi + l_cos_contact_angle * gt;
};
gpx = 0.5*(phi_w(i+1,j,k) - phi_w(i-1,j,k)); // same for y, z
nhat = grad_phi / max(|grad_phi|, 1e-8);
```

After kernel: `nhat_mf.FillBoundary()`.

**Pass 2 — curvature** $\kappa = -\nabla\cdot\hat{n}$ (`kappa_mf`, 1 comp, 1 ghost):

Only for CELL_INTERFACE cells. SOLID wall correction: for a SOLID neighbor,
use current cell's $\hat{n}$ component (zero normal-gradient of $\hat{n}$ at wall — Donath §6.3.5).

```cpp
kappa = -( 0.5*(nx_w(i+1,j,k)-nx_w(i-1,j,k))
         + 0.5*(ny_w(i,j+1,k)-ny_w(i,j-1,k))
         + 0.5*(nz_w(i,j,k+1)-nz_w(i,j,k-1)) );
```

In LB units $\Delta x=1$, no division by $\Delta x$ needed.  After kernel: `kappa_mf.FillBoundary()`.

Sign convention: $\kappa > 0$ when interface is concave toward gas (standard
bubble: gas inside → positive curvature → Laplace overpressure inside).

---

#### Step 1b — Anti-Bounce-Back (ABB) at Gas–Interface Boundary

For each **INTERFACE** cell `iv`, for each incoming direction `q` where the
source neighbor `src = iv + ev_bq` is CELL_GAS:

$$f^*_{\text{iv},q} = f^{eq}_{bq}(\rho_G, \mathbf{u}_\text{iv}, T_\text{iv})
                    + f^{eq}_{q}(\rho_G, \mathbf{u}_\text{iv}, T_\text{iv})
                    - f^\text{pre}_{\text{iv}, bq}$$

where $bq$ is the direction **outward** toward gas (opposite of $q$).

**Gas-side density** $\rho_G$ (Schwarzmeier 2023 Eq. 11–14, thermal extension):

$$\rho_G = \max\!\bigl(\rho_\text{iv} + \Delta\rho_\text{Laplace},\;10^{-3}\rho_\text{ref}\bigr), \quad
\rho_\text{iv} = \max\!\bigl(\textstyle\sum_q f^\text{pre}_{\text{iv},q},\;\rho_\text{ref}\bigr)$$

$$\Delta\rho_\text{Laplace} = \frac{-2\sigma\kappa}{R_g T_\text{iv}}, \qquad R_g = \frac{R_u}{m_\text{bar}}$$

where $\kappa$ is from Pre-Step-1b (wall-corrected), $\sigma$ = `m_fslbm_sigma`,
and $\rho_\text{ref}$ = `ic_constant.density` (read at startup).

**When σ = 0**: Δρ = 0 → ρ_G = ρ_iv = exact original Körner behaviour,
regardless of local temperature T_iv.  This is the key stability property:
no interaction between the Laplace extension and the thermal degrees of freedom
when surface tension is disabled.

**Rationale for ρ_iv ≥ ρ_ref clamp**: for `ρ_iv > ρ_ref` (impeller driving
ρ above reference), `f_pre[bq] ≈ feq_bq(ρ_actual)` — without the clamp
`f_star[q]` would go negative and compound through collide.

**Clamp to zero**: `f_star[iv][q] = max(0, ...)`.  Non-equilibrium stress near
the impeller can still push the result negative.

$T_\text{iv}$ and $\mathbf{u}_\text{iv}$ from `m_macrodata`; $T_\text{iv}$ clamped to `≥ 1e-10`.

Pressure tensor diagonal:

$$p_{xx} = u_x^2 + \frac{R_u}{m_\text{bar}} T_\text{iv}, \quad p_{yy} = u_y^2 + \frac{R_u}{m_\text{bar}} T_\text{iv}, \quad p_{zz} = u_z^2 + \frac{R_u}{m_\text{bar}} T_\text{iv}$$

#### Step 2 — Mass Flux Accumulation

Allocate `mass_flux` (MultiFab, 1 comp, initialized to zero).

For each CELL_INTERFACE cell `iv` (pull scheme, pre-streaming `f_pre`):

$$\Delta m = \sum_{q=0}^{N-1} S_q \bigl[f^\text{pre}_{(\mathbf{x}+\mathbf{c}_q),\bar{q}} - f^\text{pre}_{\mathbf{x},q}\bigr]$$

$$S_q = \begin{cases}
1 & \text{ct}(\text{ivn}) = \text{CELL\_LIQUID} \\
\tfrac{1}{2}(\varphi_\text{iv} + \varphi_\text{ivn}) & \text{ct}(\text{ivn}) = \text{CELL\_INTERFACE} \\
0 & \text{ct}(\text{ivn}) = \text{CELL\_GAS or CELL\_SOLID}
\end{cases}$$

where `ivn = iv + ev_q`.

**Diagnostic** (printed every step to stdout):
```
FSLBM rho step=NNN liq=[rho_min, rho_max] ifc=[rho_min, rho_max]
```
If `rho_max > 2.0` for interface cells, the cell index of the maximum is also printed.

#### Step 3 — Commit Streamed PDFs

```
MultiFab::Copy(m_f[lev], f_star, 0, 0, N_MICRO_STATES, nGrowVect)
m_f[lev].FillBoundary()
```

#### Step 4 — Fill Level Update

For each CELL_INTERFACE cell:

$$\varphi^{n+1} = \varphi^n + \frac{\Delta m}{\rho_\text{post}}, \quad
\rho_\text{post} = \max\!\left(\sum_q f^{n+1}_{\text{iv},q},\ 10^{-4}\right)$$

where $f^{n+1}$ is the **current post-streaming `m_f`** (written in Step 3).

**Critical**: do NOT use `m_macrodata[RHO_IDX]`.  `m_macrodata` is updated by
`collide()` which runs **after** `fslbm_advance_surface`.  Cells that were
IS_FLUID=0 last step and just re-activated have `m_macrodata[RHO]=0`, so any
finite Δm → φ = +∞, collapsing the surface within 2–3 steps.
Reading ρ directly from `m_f` eliminates this (see Section 11.7, Bug 3).

Clamp: `φ = clamp(φ, 0, 1)`.

#### Step 5a — Cell Conversion

`mass_flux` MultiFab is reset to zero and reused as a **conversion flag** scratch
(`+1` = just converted to GAS, `−1` = just converted to LIQUID).

For each CELL_INTERFACE cell:

```
if φ < FSLBM_PHI_LO:    ct → CELL_GAS,   φ → 0, f → 0, flag → +1
elif φ > FSLBM_PHI_HI:  ct → CELL_LIQUID,          flag → −1
```

Using a flag array ensures spawn (Step 5b) runs **only for cells that converted
this step**, not for every stable interface cell adjacent to liquid/gas.

#### Step 5b — Spawn New Interface Cells

After `mass_flux.FillBoundary(); m_cell_type.FillBoundary(); m_f.FillBoundary(); m_g.FillBoundary()`:

For each CELL_LIQUID or CELL_GAS cell `iv`, check face-connected neighbors
(±x, ±y, ±z only):

| Cell type | Trigger | Action |
|---|---|---|
| CELL_LIQUID | Any face-neighbor flag > +0.5 (neighbor → GAS this step) | Demote: `ct → CELL_INTERFACE`, `φ → PHI_HI`.  PDFs already valid. |
| CELL_GAS | Any face-neighbor flag < −0.5 (neighbor → LIQUID this step) | Promote: `ct → CELL_INTERFACE`, `φ → PHI_LO`.  Seed f AND g from nearest LIQUID/INTERFACE neighbor, normalized to ρ=1. |

**Normalization on spawn** (GAS → INTERFACE):

```cpp
Real rho_ivn = max(Σ_q f[ivn][q], 1e-10);
for q: f[iv][q] = f[ivn][q] / rho_ivn;
       g[iv][q] = g[ivn][q] / rho_ivn;
```

Normalizing to ρ=1 prevents inheriting an elevated density from a high-ρ LIQUID
neighbor — otherwise the spawned cell enters push streaming with ρ > 1,
compounding each conversion cycle until blow-up.

#### Step 6 — Finalize

```
fslbm_sync_isfluid_markers(lev)
m_cell_type[lev].FillBoundary()
m_phi_fslbm[lev].FillBoundary()
m_f        [lev].FillBoundary()
m_g        [lev].FillBoundary()
```

---

### 11.6 Interaction with `refill_and_spill` (as Implemented)

`refill_and_spill()` runs **before** `fslbm_advance_surface` in `advance()`.
Before the FSLBM guards were added, three categories of incorrect behavior
occurred because `refill_and_spill` uses IS_FLUID_IDX from the SDF
reconstruction, which disagrees with the FSLBM cell-type classification.

#### IS_FLUID oscillation — mass theft (Root Cause 1)

CELL_GAS cells in the headspace are tagged IS_FLUID=0 by `fslbm_sync_isfluid_markers`.
But `reconstruct_body_sdf` writes IS_FLUID=1 for headspace cells (SDF > 0.5
because the headspace is not inside the impeller).  This causes every headspace
cell to appear as *`newly_fluid`* every step.

`refill_and_spill` Steps 5 and 6 then transfer f[0] from LIQUID donors to these
"newly_fluid" gas cells, and Step 7 zeroes the donor f — draining mass from the
liquid phase continuously.

**Fix**: guard refill Steps 5 and 6:

```cpp
if (is_free_surface && (ct == CELL_GAS || ct == CELL_INTERFACE
                        || ct == CELL_LIQUID)) continue;
```

#### CELL_INTERFACE f zeroed by Step 8

`refill_and_spill` Step 8 zeroes f for cells with IS_FLUID=0.  After an
impeller step, some CELL_INTERFACE cells near the tank wall briefly have
IS_FLUID=0 (from SDF, before `fslbm_sync_isfluid_markers` corrects it).
Step 8 zeros their f → negative-ρ cells → cascade collapse.

**Fix**: skip Step 8 for CELL_INTERFACE cells:

```cpp
if (is_free_surface && ct == CELL_INTERFACE) continue;
```

#### CELL_INTERFACE f mass spilled by Step 4b

`refill_and_spill` Step 4b spills excess f from newly-solid cells to neighbors.
Some CELL_INTERFACE cells' f was incorrectly identified as "excess" and spilled
into gas cells.

**Fix**: skip Step 4b for FSLBM cells:

```cpp
if (is_free_surface && (ct == CELL_INTERFACE || ct == CELL_LIQUID
                        || ct == CELL_GAS)) continue;
```

#### Summary of all `refill_and_spill` FSLBM guards

| Step | Skip condition |
|---|---|
| Step 4b — spill to neighbors | `is_free_surface && (ct==INTERFACE\|\|LIQUID\|\|GAS)` |
| Step 5 — refill first pass (donor count) | `is_free_surface && (ct==GAS\|\|INTERFACE\|\|LIQUID)` |
| Step 6 — refill second pass (transfer) | same as Step 5 |
| Step 8 — zero solid f | `is_free_surface && ct==INTERFACE` |

---

### 11.7 Root-Cause Bug Summary

Three independent bugs all caused the free surface to blow up within 5–10 steps.
All three were fixed in commit `f1ac8de`.

#### Bug 1: `refill_and_spill` draining liquid mass into gas cells

**Mechanism**: IS_FLUID oscillation (SDF overrides `fslbm_sync`) causes gas cells
to appear `newly_fluid` → refill Steps 5&6 steal f[0] from LIQUID/INTERFACE donors
→ Step 7 zeroes donors → mass lost every step.

**Symptom**: ρ in liquid/interface cells decreases monotonically; after ~10 steps
Σf < 0 → `feq(ρ<0)` overflow.

**Fix**: Guard refill Steps 5 and 6 (Section 11.6).

#### Bug 2: `m_g` never seeded when `m_f` initialized

**Mechanism**: Case B (body exposed a new cell), Step 0 repair, and Step 5b spawn
all seeded `m_f` but left `m_g = 0`.  BGK collision computes T = Σg/(2ρcᵥ) = 0;
the entropic H-function $H \propto \ln(f/f^{eq})$ diverges with T=0.
Result: ρ = −1.64×10⁹ after one collide step.

**Symptom**: ρ in newly activated cells goes to −10⁹ on the first collide call;
surface collapses immediately.

**Fix**: Wherever f is seeded, also seed g identically:

```cpp
g[iv][q] = set_extended_equilibrium_value(1.0, {0,0,0}, pdiag, pdiag, pdiag, ...);
```

Applied in: Case B body-motion sync, Step 0 repair kernel, Step 5b GAS-promoted spawn.
(CELL_LIQUID demote in 5b keeps existing valid PDFs — no seeding needed.)

#### Bug 3: φ update used stale `m_macrodata` ρ

**Mechanism**: Step 4 originally used `m_macrodata[iv][RHO_IDX]` (from previous
step's collide).  Cells just activated with IS_FLUID=0 last step have
`m_macrodata[RHO] ≈ 0`.  Any nonzero Δm → φ = +∞ → immediate CELL_LIQUID
conversion → cascade.

**Symptom**: Surface collapses on the first step the impeller moves (Case B activation).

**Fix**: Compute ρ from the post-streaming `m_f` directly:

```cpp
Real rho_post = 0;
for (int q = 0; q < N_MICRO_STATES; ++q) rho_post += f_post[iv][q];
rho_post = max(rho_post, 1e-4);
phi += dm / rho_post;
```

---

### 11.8 Treatment of `m_g` at the Free Surface

`m_g` encodes energy: $\sum_q g_q = 2\rho c_v T$.

**Implemented**: bounce-back for `g` at gas cells (IS_FLUID_IDX==0) via the
existing `stream(lev, m_g)` call.  Gives a zero-gradient (Neumann) temperature
BC — physically correct for an adiabatic headspace.

**Optional Phase 2** (not yet implemented): proper ABB for `g`:

$$g_{\bar q}(x) = g_q^{eq}(2\rho_0 c_v T_x, \mathbf{u}_x) + g_{\bar q}^{eq}(2\rho_0 c_v T_x, \mathbf{u}_x) - g_q(x)$$

Implement only if thermal gradients at the free surface become important.

---

### 11.9 Treatment of Component Lattices at Free Surface

#### Streaming — zero-flux BC by construction

`stream(lev, m_component_lattices[i])` is called unchanged after `fslbm_advance_surface`.
The `stream()` function gates all streaming on `IS_FLUID_IDX`:

```cpp
if (is_fluid_arrs[nbx](iv, IS_FLUID_IDX) == 1) {
    if (is_fluid_arrs[nbx](ivn, IS_FLUID_IDX) != 0)
        f_star[ivn][q] = f[iv][q];          // push to fluid neighbor
    else
        f_star[iv][bounce_dirs[q]] = f[iv][q];  // bounce-back at solid/gas
}
```

Since `fslbm_sync_isfluid_markers` sets `IS_FLUID_IDX = 0` for all CELL_GAS cells,
the gas headspace automatically acts as a **solid wall** for component streaming.
Result: bounce-back → zero-flux (Neumann) BC at the free surface.

No FSLBM-specific code is needed in `stream()` or in any component handling.

#### Collision — gated on IS_FLUID_IDX

`collide()` for component lattices also gates on `IS_FLUID_IDX`:

```cpp
if (is_fluid_arrs[nbx](iv, IS_FLUID_IDX) == 1) {
    // BGK: f_eq_comp[q] = rho_comp * eq_unit[q]
    //      f += omega_comp * (f_eq_comp - f)
}
```

Gas cells (IS_FLUID=0) are skipped entirely — their component PDFs remain at zero.

`omega_comp = 1.0 / (D / (R/m_bar * T * dt) + 0.5)` uses `T` from `m_macrodata`.
For newly spawned interface cells, T is valid because Step 0 and Case B seeding
now initialize `g` correctly (fixing T) before `collide()` runs each step.

#### Cell transitions and component PDFs

| Transition | Component PDF behavior |
|---|---|
| GAS → INTERFACE (spawn, Step 5b) | **Left at zero** — cell was gas, had no dissolved O₂. Correct: net zero concentration at a newly created surface cell; O₂ diffuses in from adjacent liquid over subsequent steps. |
| LIQUID → INTERFACE (demotion, Step 5b) | **Left unchanged** — cell's dissolved O₂ is preserved. Correct: the cell was already carrying O₂ as a liquid cell. |
| INTERFACE → GAS (conversion, Step 5a) | Component PDFs zeroed together with `f` (all zeroed in the `for q` loop). Any O₂ in the cell is lost — negligible in practice since PHI<1e-4 means the cell barely contributed. |
| INTERFACE → LIQUID | No change to component PDFs. |
| SOLID → LIQUID/GAS (Case B, body sweep out) | Component lattice **not seeded** — starts at zero. Correct: the cell was inside the impeller, carried no dissolved O₂. |

#### Note on hard-coded value in entropic solve (`1.0`)

`collide()` builds a unit-density equilibrium shape for the entropy reference:

```cpp
eq_unit_arr(iv, q+NQ) = set_extended_equilibrium_value(
    1.0, zero_vel, p_by_rho, p_by_rho, p_by_rho, ...);
```

This `1.0` is **correct and intentional** — it is the normalisation density for the
H-function `H = Σ f_q ln(f_q / f_ref_q)`.  Using ρ=1 makes the entropy reference
independent of the actual (variable) `rho_comp`, which is required for the
Ansumali–Karlin entropic alpha solve.  This is not a physical density assumption.

#### Future: Dirichlet BC at free surface (headspace desorption)

If desorption into the headspace becomes important, replace bounce-back with ABB:

$$f^*_{\text{iv},q} = f^{eq}_{bq}(\rho_{\text{comp},G}, \mathbf{u}_\text{iv}, T_\text{iv})
+ f^{eq}_{q}(\rho_{\text{comp},G}, \mathbf{u}_\text{iv}, T_\text{iv})
- f^\text{pre}_{\text{iv},bq}, \quad \rho_{\text{comp},G} = 0$$

This gives a Dirichlet (zero O₂ concentration) BC at the gas–liquid interface.
Not needed for the 10 s kLa benchmark where surface desorption is negligible.

---

### 11.10 Code File Map (as Implemented)

| File | Location | Content |
|---|---|---|
| `Source/LBM.H` | — | `m_phi_fslbm`, `m_cell_type`, `m_free_surface`, `m_free_surface_z`, `m_fslbm_sigma`, `m_fslbm_contact_angle_deg`; prototypes for `fslbm_advance_surface`, `fslbm_init_cell_type`, `fslbm_sync_isfluid_markers` |
| `Source/Constants.H` | — | `CELL_SOLID/GAS/INTERFACE/LIQUID`; `FSLBM_PHI_LO`; `FSLBM_PHI_HI` |
| `Source/LBM.cpp` | line 5443 | `LBM::fslbm_sync_isfluid_markers(lev)` |
| `Source/LBM.cpp` | line 5577 | `LBM::fslbm_init_cell_type(lev)` |
| `Source/LBM.cpp` | line 5657 | `LBM::fslbm_advance_surface(lev)` — incl. Pre-Step-1b (nhat/kappa) |
| `Source/LBM.cpp` | `advance()` | call order: refill_and_spill → fslbm_advance_surface → stream(components) → stream(m_g) → collide |
| `Source/LBM.cpp` | `refill_and_spill()` | FSLBM guards in Steps 4b, 5, 6, 8 |

---

### 11.11 Implementation Phases Status

**Phase 1 — Minimum viable** ✅ Complete (commit `f1ac8de`)

- [x] `m_cell_type` field, constants, allocation
- [x] `fslbm_init_cell_type`: sharp-interface initialization from `m_free_surface_z`
- [x] `fslbm_sync_isfluid_markers`: IS_FLUID/EB_BOUNDARY/IS_FLUID_SIDE from cell type
- [x] `fslbm_advance_surface`: body-motion sync, Step 0 repair, Steps 1a/1b/2/3/4/5a/5b/6
- [x] ABB at gas→interface with local ρ (not hardcoded ρ=1) and result clamped ≥ 0
- [x] φ update from Σf_post (not stale macrodata)
- [x] m_g seeded alongside m_f everywhere f is initialized
- [x] `refill_and_spill` FSLBM guards (Steps 4b, 5, 6, 8)
- [x] 50-step stability test: Δm ≈ machine epsilon, interface ρ ∈ [0.9, 1.1]

**Phase 1b — Surface tension + contact angle** ✅ Complete (commit `5d270e3`)

- [x] `m_fslbm_sigma`: input `lbm.fslbm_sigma`; LB surface tension
- [x] `m_fslbm_contact_angle_deg`: input `lbm.fslbm_contact_angle`; default 90°
- [x] Pre-Step-1b: `nhat_mf` (Pass 1) — $\hat{n}=\nabla\varphi/|\nabla\varphi|$ with contact-angle ghost-φ correction
- [x] Pre-Step-1b: `kappa_mf` (Pass 2) — $\kappa=-\nabla\cdot\hat{n}$ with SOLID wall correction
- [x] Step 1b ABB: additive Laplace correction $\Delta\rho = -2\sigma\kappa/(R_g T)$
- [x] Stability property: when $\sigma=0$, $\Delta\rho=0$ → exact Körner behaviour regardless of $T$
- [x] 18 000+ step test (kla_bioreactor, $\sigma_\text{LB}=1.026\times10^{-5}$): stable

**Phase 2 — Robust surface (not yet started)**

- [ ] Excess-mass redistribution (`fslbm_distribute_excess_mass`)
- [ ] Push/pull forced-fill/empty for isolated interface cells (Thürey 2007)
- [ ] Bubble exit detection via `cell_type == CELL_GAS`
- [ ] 500-step test with impeller spinning up

**Phase 3 — Optional**

- [ ] ABB for `m_g` at free surface (Section 11.8)
- [ ] ABB for component lattices with `ρ_comp = 0` (head-space desorption)

---

### 11.12 Key Equations Reference

**Interface normal** (Pre-Step-1b, Pass 1; Attar & Körner 2009):

$$\hat{n} = \frac{\nabla_\text{CA}\varphi}{|\nabla_\text{CA}\varphi|}, \qquad
(\nabla_\text{CA}\varphi)_x = \tfrac{1}{2}\bigl(\varphi_\text{ghost}(i{+}1,j,k)-\varphi_\text{ghost}(i{-}1,j,k)\bigr)$$

$$\varphi_\text{ghost}(ii,jj,kk) = \begin{cases}
\varphi + \cos(\theta_W)\,|\nabla_\perp\varphi| & \text{SOLID neighbor} \\
\varphi(ii,jj,kk) & \text{otherwise}
\end{cases}$$

**Curvature** (Pre-Step-1b, Pass 2):

$$\kappa = -\nabla_\text{wall}\cdot\hat{n}, \quad
(\nabla_\text{wall})_x \hat{n}_x = \tfrac{1}{2}\bigl[n_{x,w}(i{+}1,j,k)-n_{x,w}(i{-}1,j,k)\bigr]$$

$$n_{x,w}(ii,jj,kk) = \begin{cases} \hat{n}_x(i,j,k) & \text{SOLID neighbor} \\ \hat{n}_x(ii,jj,kk) & \text{otherwise}\end{cases}$$

**ABB for f with Laplace pressure** (Körner 2005 / Schwarzmeier 2023 Eq. 11-14, thermal):

$$f^*_{\text{iv},q} = f^{eq}_{bq}(\rho_G, \mathbf{u}_\text{iv}, T_\text{iv})
                    + f^{eq}_{q}(\rho_G, \mathbf{u}_\text{iv}, T_\text{iv})
                    - f^\text{pre}_{\text{iv},bq}$$

$$\rho_G = \max\!\bigl(\rho_\text{iv} + \Delta\rho_\text{Laplace},\; 10^{-3}\rho_\text{ref}\bigr), \quad
\rho_\text{iv} = \max\!\bigl(\textstyle\sum_q f^\text{pre}_{\text{iv},q},\; \rho_\text{ref}\bigr)$$

$$\Delta\rho_\text{Laplace} = \frac{-2\sigma\kappa}{R_g T_\text{iv}}, \quad R_g = \frac{R_u}{m_\text{bar}}$$

When $\sigma=0$: $\Delta\rho=0 \Rightarrow \rho_G=\rho_\text{iv}$ (exact original Körner behaviour).

**Mass flux** (Step 2):

$$\Delta m = \sum_q S_q\bigl[f^\text{pre}_{(\mathbf{x}+\mathbf{c}_q),\bar{q}} - f^\text{pre}_{\mathbf{x},q}\bigr],\quad
S_q = \begin{cases}1 & \text{LIQUID nbr}\\\tfrac{\varphi_\mathbf{x}+\varphi_{\mathbf{x}+\mathbf{c}_q}}{2} & \text{INTERFACE nbr}\\0 & \text{GAS/SOLID}\end{cases}$$

**Fill level update** (Step 4, as implemented):

$$\varphi^{n+1} = \varphi^n + \frac{\Delta m}{\max\!\left(\sum_q f^{n+1}_{\text{iv},q},\ 10^{-4}\right)}$$

(Körner 2005 Eq. 6 writes $(\varphi^n \rho^\text{pre} + \Delta m)/\rho^\text{post}$;
the implemented simplified form drops the $\varphi^n \rho^\text{pre}$ numerator
correction — valid for small per-step ρ variation in the bioreactor regime.)

**Conversion**: $\varphi < 10^{-4} \Rightarrow \text{CELL\_GAS}$;
$\varphi > 1-10^{-4} \Rightarrow \text{CELL\_LIQUID}$

**References**:
- Körner et al. (2005) J. Stat. Phys. 121: original FSLBM
- Schwarzmeier et al. (2023) JCP 473, 111753: FSLBM vs PFLBM benchmark, surface tension ABB Eq. 11–14
- Thürey (2007) PhD thesis, Univ. Erlangen: excess mass, push/pull rules
- Attar & Körner (2009) J. Colloid Interface Sci. 335: contact angle via ghost-φ
- Donath (2011) PhD thesis, Univ. Erlangen: obstacle normal, wall correction for κ (§6.3.5)
- waLBerla: `src/lbm/free_surface/dynamics/functionality/` (ABB, AdvectMass, CellConversion)

---

## 12. Open Decisions (kLa modules — not yet started)

| Question | Options | Recommendation |
|---|---|---|
| Particle container GPU/CPU? | CPU `ParticleContainer` vs GPU-offloaded | Start CPU; move to GPU if bubble count > 10⁵ |
| SGS turbulence closure? | None vs. Smagorinsky $C_s=0.1$ | **Include from start** — required for correct ε and hence $k_L$ |
| Bubble-fluid force interpolation? | Nearest-cell vs. trilinear | **Trilinear** (Thomas et al. explicit) |
| O₂ diffusivity in LB units | Sc=500 → `diffusivity≈1.6e-9` LB | Use entropic component solver (already in codebase) |
| Daughter bubble distribution? | Equal, random, triangle(mode=0.2) | **Triangle(0, 0.5, mode=0.2)** — Thomas et al.; others for sensitivity |
| Energy balance validation? | Optional vs. mandatory | **Mandatory** before enabling mass transfer — $P_s = P_d$ to <5% |
