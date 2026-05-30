# Thermal FSLBM Component Boundaries & Outgassing Plan

This document outlines the proposed theoretical framework and implementation for managing `m_component_lattices` (dissolved scalars, e.g., oxygen) at the dynamic Free-Surface LBM (FSLBM) boundary. 

Currently, the component lattices simply stream $0.0$ from gas neighbors. This acts as an artificial vacuum sink. To transition the simulation from an artificial vacuum to a controlled kLa interfacial mass-transfer model (matching the bubble solvers), we intend to perform the following steps:

## Step A: The Impermeable Lid (Zero-Gradient / Bounce-Back)
We must first prevent simple outgassing due to lattice streaming. When missing directions are advected from the `CELL_GAS` vacuum into the `CELL_INTERFACE` components, they should perform a localized bounce-back. This creates a zero-mass-flux (impermeable lid) boundary. 

```cpp
// This would be wrapped in a function `fslbm_replenish_components(int lev)`
// In a loop over `amrex::ParallelFor(m_component_lattices[c][lev], ...)`
if (ct_arrs[nbx](src, 0) == CELL_GAS) {
    // Missing population due to gas boundary -> purely bounce back outgoing 
    // population to enforce 0-gradient impermeability
    f_comp_arrs[nbx](iv, q) = f_comp_arrs[nbx](iv, bounce_dirs[q]);
}
```

*Note: In the current AMD/LBM branch, basic stream() operations actually contain logic that bounces populations back from solid/gas cells. This function would guarantee complete impermeability against gas specifically if standard boundary conditions drop it.*

## Step B: Spawning Components from Gas $\rightarrow$ Interface
When the total volume expands and a cell converts from `GAS` to `CELL_INTERFACE`, the new interface cell presently has $0.0$ oxygen concentration. This behaves as an artificial dilution of the scalar fields.

Whenever $f$ and $g$ populations are spawned from a fluid neighbor in Step 5 (Mass Conversion), the component lattices should do a pure non-equilibrium 1:1 copy of the exact same neighbor.

```cpp
// Within fslbm_advance_surface() after checking if a neighbor converted to GAS:
// ...
if (ctn == CELL_LIQUID || ctn == CELL_INTERFACE) {
    // [Existing logic copies and normalizes f and g]
    const amrex::Real inv_rho = l_fslbm_rho_ref / rho_ivn;
    for (int q = 0; q < N_MICRO_STATES; ++q) {
        f_arrs[nbx](iv, q) = f_arrs[nbx](ivn, q) * inv_rho;
        g_arrs[nbx](iv, q) = g_arrs[nbx](ivn, q) * inv_rho;
    }
    // [New Logic: Mark this cell as needing component copies]
    flag_arrs[nbx](iv, 0) = amrex::Real(-2.0);
    return;
}
//...
// Then in a subsequent component specific loop:
if (flag_ro[nbx](iv, 0) < amrex::Real(-1.5)) {
    // Find valid neighbor and directly copy to preserve scalar C:
    for (int q = 0; q < N_MICRO_STATES; ++q) {
        c_arrs[nbx](iv, q) = c_arrs[nbx](ivn, q);
    }
}
```

## Step C: Controlled Outgassing via kLa
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
