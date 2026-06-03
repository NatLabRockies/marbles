# Thermal FSLBM Component Boundaries & Outgassing Plan

This document outlines the proposed theoretical framework and implementation for managing `m_component_lattices` (dissolved scalars, e.g., oxygen) at the dynamic Free-Surface LBM (FSLBM) boundary. 

## Completed: Step A (The Impermeable Lid / Bounce-Back)
**Status: Implemented and Re-Validated**
*(Update: Standard `stream()` inherently bounces populations from exact `CELL_GAS` domains since they are marked `IS_FLUID=0`. Any custom bounce-back injection actually overrides and duplicates this natural conservation, causing extreme scalar imbalances. The custom boundary wrapper has been carefully stripped out entirely. Standard advective bounce-back preserves the 0-gradient condition completely naturally and keeps mass conserved.)*

## Completed: Step B (Spawning Components & Absolute Safeties)
**Status: Implemented and Stabilized**

**The Mass Cloning Pitfall:**
When a cell volume expands and converts from `GAS` to `CELL_INTERFACE`, it previously triggered a search for fluid neighbors. Making exact fraction clones of adjacent variables (`c_arrs(iv) = c_arrs(ivn)`) or filling it with equilibrium $1.0^{LB}$ weights (to prevent zero-sinks) directly *duplicates and injects mass out of nowhere*, causing explosive non-physical macroscopic buildup ($Y_0 \rightarrow 10^{18}$). Since `m_component_lattices` operates as an Eulerian density variable independently tracking mass $\sum c_{arr}$, standard PUSH streaming natively streams mass into the newly formed empty space. 

**Corrected Implementation:**
1. **Absolute Vacuums allowed safely:** When `CELL_GAS` converts to `CELL_INTERFACE`, the scalar components are uniformly initialized to exactly **$0.0$**. The subsequent `stream()` step naturally pushes fluid into the empty cell from interior liquid neighbors, achieving an organic, mass-conserving interpolation gradient without inventing mass. 
2. **Collision Temperature Caps:** The reason a `0.0` distribution originally crashed the solver was because sub-grid droplet anomalies drove the *main free-surface fluid* temperature up to infinity (`2.1685e+293`), scrambling `eq_ref_q`. Both `macrodata_to_equilibrium` and `relax_f_to_equilibrium` entropic functions now deploy an absolute $T \leq 0.5$ thermodynamic safeguard threshold against unrecoverable acoustic surges: `T_safe = amrex::min(temperature, 0.5 / specific_gas_constant)`. The BGK loop naturally skips sparse states $rho\_comp \le 1.0e-30$ without issues now that the primary temperature is rigorously capped.

## Next Action: Step C (Controlled FSLBM Outgassing via kLa)
**Status: Pending**

Once the surface is impermeable via Step A and B, mass only leaves the solver by explicit localized physical source terms. We compute the outgassing delta identically to the bubble module.
Since our LB density unit for component $O_2$ is already defined as $1.0^{LB} = C_{ref}$, the outgassing source term simplifies.

```cpp
// Evaluated only for CELL_INTERFACE facing GAS
const amrex::Real d_rho_LB = (k_La_surface * dt_phys) * (1.0 - rho_comp_LB);

// Inject this delta along equilibrium weights (borrowing l_T_ref assumption)
for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
    const amrex::Real f_eq_unit = set_extended_equilibrium_value(...);
    fO2_arrs[nbx](iv, q) += d_rho_LB * f_eq_unit; 
}
```
