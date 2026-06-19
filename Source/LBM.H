#ifndef LBM_H
#define LBM_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BCUtil.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_ErrorList.H>
#include "EB.H"
#include "BC.H"
#include "Constants.H"
#include "FillPatchOps.H"
#include "IC.H"
#include "Utilities.H"
#include "Stencil.H"
#include "Bubble.H"

namespace lbm {

class LBM : public amrex::AmrCore
{
public:
    LBM();
    ~LBM() override;

    // initializes multilevel data
    void init_data();

    // advance solution to final time
    void evolve();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping
    // and fill with existing fine and coarse data. overrides the pure virtual
    // function in AmrCore
    void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    void ClearLevel(int lev) override;

    // Make a new level from scratch using provided BoxArray and
    // DistributionMapping. Only used during initialization. overrides the pure
    // virtual function in AmrCore
    void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

    // compute dt
    static amrex::Real est_time_step(const int lev);

    // convert macrodata to equilibrium
    void macrodata_to_equilibrium(const int lev);

    // calculate the macro fluid properties from the distributions and store
    // part of (4.45) of http://dx.doi.org/10.3929/ethz-b-000607045
    void f_to_macrodata(const int lev);

    // Stream the information to the neighbor particles
    void stream(const int lev, amrex::Vector<amrex::MultiFab>& fs);

    // Collide the particles
    void collide(const int lev);

    // Relax the particles toward the equilibrium state
    void relax_f_to_equilibrium(const int lev);

    // Compute derived quantities
    void compute_derived(const int lev);

    // Compute corrections
    void compute_q_corrections(const int lev);

    // Compute forces on EB
    void compute_eb_forces();

    // Initialize microdata
    void initialize_f(const int lev);

    // Initialize is_fluid
    void initialize_is_fluid(const int lev);

    // Initialize moving body shape from file
    void initialize_moving_body_shape(int lev);

    // Initialize stationary body from file (STL/CSV)
    void init_stationary_body(int lev);

    // Update integer is_fluid mask from the fractional field and recompute
    // boundary markers (steps 1 and 2 of refill_and_spill preparation).
    // If threshold is omitted, class default m_is_fluid_fraction_threshold is used.
    void update_is_fluid_from_fraction_and_mark(const int lev,
                                               amrex::Real threshold = -1.0);

    // Refill and spill algorithm (steps 3 and 4). This does the following
    // (per-level):
    //  - snapshot previous is_fluid (component 0), threshold fractional
    //    field into integer mask, recompute boundary markers (EB and side
    //    boundaries), then
    //  - SPILL: for cells that became solid, distribute their m_f/m_g
    //    populations to neighboring secondary-boundary (component 3)
    //    cells into matching population components using stencil weights
    //    (fallback to component 2 if needed). Donors are zeroed.
    //  - REFILL: for cells that became fluid, collect immediate boundary
    //    neighbors (component 2) and compute weighted sums to fill each
    //    population component; subtract the same weighted amounts from
    //    donor neighbors so overall conservation holds.
    // If threshold is negative, the class default is used.
    void refill_and_spill(const int lev, amrex::Real threshold = -1.0);

    // Reconstruct fractional mask from body signed distance function.
    // For rigid body motion (translation + rotation), this updates
    // m_is_fluid_fraction by evaluating the body geometry in its current
    // position/orientation without numerical diffusion.
    void reconstruct_body_sdf(const int lev, amrex::Real time);

    // Update m_body_angular_velocity from the target value and the ramp
    // schedule.  No-op when m_body_angular_velocity_ramp_steps == 0
    // (legacy instant-on behaviour).  Called at the start of advance()
    // before any kernel captures m_body_angular_velocity.
    void update_body_angular_velocity_for_ramp();

    // Advance the conservative phase-field (Chiu & Lin 2011, Eq. 18) for
    // the free liquid-gas surface by one LBM time step.
    // Evolves m_is_fluid_fraction in place; calls refill_and_spill internally.
    void advance_phi(const int lev);

    // FSLBM (Körner 2005): initialize cell type + fill level from z_surf.
    void fslbm_init_cell_type(const int lev);

    // FSLBM: one free-surface sweep (streaming with ABB + φ update + conversions).
    // Replaces both advance_phi(lev) and stream(lev, m_f).
    void fslbm_advance_surface(const int lev);

    // FSLBM: replenish m_g populations in INTERFACE cells facing gas.
    // Prevents thermal energy drain at the free surface.
    void fslbm_replenish_g(const int lev);
    void fslbm_replenish_components(const int lev);

    // FSLBM: derive IS_FLUID_IDX (and boundary markers) from m_cell_type.
    // Replaces update_is_fluid_from_fraction_and_mark for FSLBM to avoid
    // conflict with the moving-body SDF that also lives in m_is_fluid_fraction.
    void fslbm_sync_isfluid_markers(const int lev);

    // Initialize masks
    void initialize_mask(const int lev);

    // Fill microdata inside the EB
    void fill_f_inside_eb(const int lev);

    // get a field based on a variable name
    std::unique_ptr<amrex::MultiFab>
    get_field(const std::string& name, const int lev, const int ngrow);

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> plot_file_mf();

    // clamp negative component densities to zero
    void clamp_negative_component_densities(const int lev);

    // Apply operator-split reaction source terms for the two-step catalytic
    // mechanism  S+C <-> I -> P+C  to all component lattices on one level.
    // Source increments are distributed uniformly across all N_MICRO_STATES
    // populations so that no net momentum is added by the chemistry.
    void apply_reaction_source_terms(int lev);

    // Inject the catalyst (component 1) into a box region at a prescribed
    // simulation step.  Executed only once (guarded by m_cat_inject_done).
    void apply_timed_catalyst_injection(int lev);

    // Open/close the per-step species statistics CSV file.
    void open_species_stats_file(bool initialize);
    void close_species_stats_file();

    // Append one row (step, LB-time, mean rho per species, conservation
    // residuals) to the species statistics file.
    void write_species_stats();

private:
    // read parameters
    void read_parameters();

    // read tagging parameters
    void read_tagging_parameters();

    // a wrapper for EstTimeStep(0
    void compute_dt();

    // advance a level by dt
    void time_step(const int lev, const amrex::Real time, const int iteration);

    // advance a single level for a single time step, updates flux registers
    void advance(
        const int lev,
        const amrex::Real time,
        const amrex::Real dt_lev,
        const int iteration,
        const int ncycle);

    // perform work after a time step
    void post_time_step();

    // Check if a field exists
    bool check_field_existence(const std::string& name);

    // Get field component
    static int get_field_component(
        const std::string& name, const amrex::Vector<std::string>& varnames);

    // set covered coarse cells to be the average of overlying fine cells
    void average_down(amrex::IntVect crse_ng);

    // more flexible version of AverageDown() that lets you average down across
    // multiple levels
    void average_down_to(int crse_lev, amrex::IntVect crse_ng);

    // Check distribution
    void sanity_check_f(const int lev);

    // Set the user defined BC functions
    void set_bcs();

    // Set the user defined IC functions
    void set_ics();

    // get plotfile name
    std::string plot_file_name(const int step) const;

    // get chkfile name
    std::string chk_file_name(const int step) const;

    // set plotfile variables names
    amrex::Vector<std::string> plot_file_var_names() const;

    // write plotfile to disk
    void write_plot_file();

public:
    // write checkpoint file to disk
    void write_checkpoint_file() const;

    // read checkpoint file from disk
    void read_checkpoint_file();
private:

    // utility to skip to next line in Header
    static void goto_next_line(std::istream& is);

    // open the forces file
    void open_forces_file(const bool initialize);

    // close the forces file
    void close_forces_file();

    // write the forces to file
    void output_forces_file(const amrex::Vector<amrex::Real>& forces);

    amrex::Vector<int> m_isteps;    // which step?
    amrex::Vector<int> m_nsubsteps; // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> m_ts_new;
    amrex::Vector<amrex::Real> m_ts_old;
    amrex::Vector<amrex::Real> m_dts;

    // maximum number of steps and stop time
    int m_max_step = std::numeric_limits<int>::max();
    amrex::Real m_stop_time = std::numeric_limits<amrex::Real>::max();

    // variable names for output
    amrex::Vector<std::string> m_lbm_varnames;
    amrex::Vector<std::string> m_macrodata_varnames;
    amrex::Vector<std::string> m_microdata_varnames;
    amrex::Vector<std::string> m_microdata_g_varnames;
    amrex::Vector<std::string> m_deriveddata_varnames;
    amrex::Vector<std::string> m_idata_varnames;
    amrex::Vector<std::string> m_fracdata_varnames;

    // array of multifabs to store the solution at each level of refinement
    amrex::Vector<amrex::MultiFab> m_macrodata;
    amrex::Vector<amrex::MultiFab> m_f;
    amrex::Vector<amrex::Vector<amrex::MultiFab>> m_component_lattices;
    amrex::Vector<amrex::MultiFab> m_g;
    amrex::Vector<amrex::MultiFab> m_eq;
    amrex::Vector<amrex::MultiFab> m_eq_g;
    amrex::Vector<amrex::MultiFab> m_derived;
    amrex::Vector<amrex::iMultiFab> m_is_fluid;
    amrex::Vector<amrex::MultiFab> m_plt_mf;
    amrex::Vector<amrex::iMultiFab> m_mask;
    
    //non-integer mask to identify fluid cells. Will be moved with advection equation for moving solid body
    amrex::Vector<amrex::MultiFab> m_is_fluid_fraction;

    // Moving body reference data
    amrex::Gpu::DeviceVector<uint16_t> m_body_voxel_data;
    amrex::IntVect m_body_voxel_dims{0};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_body_voxel_origin{0.0};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_body_voxel_dx{0.0};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_body_initial_center{0.0};
    bool m_using_voxel_body = false;

    // Stationary body reference data
    amrex::Vector<amrex::iMultiFab> m_stationary_mask;
    bool m_has_stationary_body = false;

    // FSLBM cell type per cell: CELL_SOLID(0), CELL_GAS(1), CELL_INTERFACE(2), CELL_LIQUID(3)
    amrex::Vector<amrex::iMultiFab> m_cell_type;

    // FSLBM fill level φ ∈ [0,1] — separate from m_is_fluid_fraction which is
    // used by the moving-body SDF / refill_and_spill and must not be modified
    // by FSLBM code.
    amrex::Vector<amrex::MultiFab> m_phi_fslbm;
 
    const int m_macrodata_nghost = 1;
    const int m_f_nghost = 3;
    const int m_eq_nghost = 0;
    const int m_derived_nghost = 0;

    amrex::Vector<int> m_bc_lo{0, 0, 0};
    amrex::Vector<int> m_bc_hi{0, 0, 0};
    amrex::GpuArray<int, AMREX_SPACEDIM * 2> m_bc_type = {0};
    std::string m_velocity_bc_type{"noop"};

    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<amrex::BCRec> m_bcs;

    std::unique_ptr<FillPatchOpsBase> m_fillpatch_op;
    amrex::Vector<std::unique_ptr<FillPatchOpsBase>> m_component_fillpatch_ops;
    std::unique_ptr<FillPatchOpsBase> m_fillpatch_g_op;

    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>>
        m_factory;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int m_regrid_int = 2;

    // if >= 0 we restart from a checkpoint
    std::string m_restart_chkfile;

    // plotfile prefix and frequency
    std::string m_plot_file{"plt"};

    int m_plot_int = -1;

    // screen output frequency (steps between step-header prints, default = plot_int)
    int m_print_int = 1000;

    // checkpoint prefix and frequency
    std::string m_chk_file{"chk"};
    int m_chk_int = -1;

    // digits used in chk and plt file names
    int m_file_name_digits{5};

    // error tags
    amrex::Vector<amrex::AMRErrorTag> m_err_tags;

    // User-input quantities (will be overridden by ParmParse)
    amrex::Real m_dx_outer = 1.0;
    amrex::Real m_dt_outer = 1.0;
    amrex::Real m_nu = 1.0;
    amrex::Real m_alpha = 1.0;
    int m_n_components = 0;
    amrex::Vector<amrex::Real> m_component_diffusivities;

    // Physical (SI) cell size and time step.  Default 1.0 means "already in LB
    // units"; override via lbm.dx_phys [m] and lbm.dt_phys [s] when the input
    // file contains a body force in physical units (e.g. lbm.gravity in m/s^2).
    // These also feed BubbleParams (m_bubble_params.dx_phys, .dt_phys).
    amrex::Real m_dx_phys = 1.0;
    amrex::Real m_dt_phys = 1.0;

    // External body force per unit mass in physical units [m/s^2].
    // Converted to LB acceleration  g_LB = g_phys * dt_phys^2 / dx_phys  at the
    // use site (apply_macroscopic_forcing).  Default = no gravity.
    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, 0.0};

    // save streaming data to plot file
    bool m_save_streaming = true;

    // save derived data to plot file
    bool m_save_derived = true;

    // compute and output forces
    bool m_compute_forces = false;

    // clamp negative component densities before collision
    bool m_clamp_component_densities = false;

    // Mass-fraction window [Y_min, Y_max] for the component-lattice clamp.
    // Cells with Y_c = rho_comp / rho_main outside this window are
    // *rescaled* (or zeroed if rho_main is too small) to drag the zeroth
    // moment back into the physical band [0, 1] with a numerical safety
    // margin.  Cells INSIDE the window are left untouched — small negative
    // excursions produced by the entropic Newton solve / streaming
    // round-off during start-up transients must be preserved (zeroing
    // them was found to break mass conservation in the single-phase
    // scalar transport validation tests).
    // Defaults [-0.1, 1.1] give 10× headroom on each side of the
    // physically valid range Y ∈ [0, 1].
    // Set via lbm.component_y_min / lbm.component_y_max.
    amrex::Real m_clamp_component_y_min = -0.1;
    amrex::Real m_clamp_component_y_max =  1.1;

    // Carrier-fluid density floor below which the component populations
    // in the cell are zeroed unconditionally — "mass fraction" is
    // meaningless when there is no carrier mass (gas / solid cells).
    // Default 1e-6.  Set via lbm.component_rho_main_floor.
    amrex::Real m_clamp_component_rho_main_floor = 1.0e-6;

    // use entropic alpha stabilisation for component lattices
    bool m_use_entropic_components = false;

    // use entropic alpha stabilisation for the main f distribution
    bool m_use_entropic_f = false;

    // ---------------------------------------------------------------
    // Lagrangian bubble mass transfer physics
    // ---------------------------------------------------------------

public:
    /// Apply macroscopic body force to the fluid distributions m_f and m_g
    /// using the exact-difference forcing scheme:
    ///     delta_f_q = f_eq(rho, u + F*dt/rho, T) - f_eq(rho, u, T)
    ///     delta_g_q = g_eq(2*rho*e1, ...)        - g_eq(2*rho*e0, ...)
    /// where the energy moments are recomputed with the shifted velocity.
    /// Skips CELL_INTERFACE cells (Donath 2011, p122) to prevent acoustic
    /// shocks at the FSLBM boundary.  Always adds gravity m_gravity (in LB
    /// units g_LB = g_phys * dt_phys^2 / dx_phys); if force_mf != nullptr
    /// also adds the bubble-induced acceleration (capped at 50x gravity).
    void apply_macroscopic_forcing(int lev, const amrex::MultiFab* force_mf);

    /// Apply bubble O2 source to the dissolved-O2 component lattice.
    /// o2_src_mf: 1-component source rate in mol/(m³·s).
    void apply_bubble_o2_source(int lev, const amrex::MultiFab& o2_src_mf);

private:

    // file for forces
    std::string m_forces_file{"forces.txt"};

    // file handle for forces
    std::fstream m_forces_stream;

    // Derived values (will be updated after reading in user inputs)
    amrex::Real m_mesh_speed = 1.0;
    amrex::Real m_cs = 1.0;
    amrex::Real m_cs_2 = 1.0;

    std::unique_ptr<ic::InitializerBase> m_ic_op;
    amrex::Vector<std::unique_ptr<ic::InitializerBase>> m_component_ic_ops;
    std::string m_ic_type;

    // LBM Scheme selector. "energyD3Q27" activates product equilibrium,
    // energy equation etc.
    std::string m_model_type = "energyD3Q27";

    // Initial temperature.
    amrex::Real m_initialTemperature = 1.0 / 3.0;

    // threshold for converting fractional mask into integer is_fluid
    amrex::Real m_is_fluid_fraction_threshold = 0.5;

    // ---------------------------------------------------------------
    // Free-surface (conservative phase-field, Chiu & Lin 2011)
    // ---------------------------------------------------------------
    // Enable by setting  lbm.free_surface = 1  in the input file.
    // m_is_fluid_fraction plays the role of Phi in [0,1]:
    //   Phi = 1  liquid cells  (LBM is active)
    //   Phi = 0  gas headspace (LBM is inactive — like a moving solid)
    // The initial flat interface at  z = m_free_surface_z  (LB cells)
    // is smeared over ~1.5 cells with a tanh profile.
    bool        m_free_surface       = false;
    amrex::Real m_free_surface_z     = 0.0;   // interface z, in LB cell index units
    amrex::Real m_free_surface_eps   = 0.0;   // interface half-width  (~1.5*dx, set by init)
    amrex::Real m_phi_M0             = 0.0;   // initial total Phi mass for mass redistribution
    amrex::Real m_phi_gamma_coeff    = 1.0;   // multiplier on gamma = gamma_coeff * |u|
    // FSLBM reference density: read from ic_constant.density at startup.
    // Used for ABB ρ₀, repair/spawn seeding (feq target), and repair threshold.
    // Defaults to 1.0 (the conventional LBM reference) but scales correctly when
    // the initial bulk density differs from 1.
    amrex::Real m_fslbm_rho_ref      = 1.0;
    // Surface tension in LB units for the ABB gas-pressure BC.
    // p_G = p_atm - 2*sigma*kappa;  rho_G = p_G / (Rg*T_interface).
    // Default 0 = flat interface (no Laplace pressure).
    // Set via lbm.fslbm_sigma in the input file.
    amrex::Real m_fslbm_sigma        = 0.0;
    // Static contact angle for all solid walls [degrees].
    // 90 = neutral wetting (default, same as original Körner model).
    // < 90 = hydrophilic  ;  > 90 = hydrophobic.
    // Set via lbm.fslbm_contact_angle in the input file.
    amrex::Real m_fslbm_contact_angle_deg = 90.0;

    // Stranded-cell sweep search radius (LB cells).
    // FSLBM Step 0_pre converts any CELL_INTERFACE cell that has no
    // CELL_LIQUID in a (2R+1)^3 box of neighbours to CELL_GAS
    // (zeroes φ, f, g).  This is a purely topological criterion:
    // "cannot drain to bulk fluid via streaming within R steps".
    // Default R=2 (5x5x5 box = 124 neighbours) catches isolated
    // droplets and detached splash sheets two cells from any liquid.
    // R=1 reduces to a strict D3Q27 isolation test (very conservative).
    // R must not exceed m_f_nghost (=3) so reads stay inside the ghost
    // layer.  Set to 0 to disable the sweep entirely.
    // Set via lbm.fslbm_strand_search_radius in the input file.
    int m_fslbm_strand_search_radius = 2;

    // Dirichlet thermal boundary condition at FSLBM interface cells:
    //   T(CELL_INTERFACE) = T_ref  (= m_initialTemperature)
    // Physical justification: the gas headspace (and any free-surface
    // interface cell facing it) is in thermal equilibrium with the
    // atmospheric reservoir.  Surface tension / vapour transport
    // homogenises temperature across the interface much faster than
    // the lattice can resolve.
    //
    // Implementation: after the normal collision step, for each
    // CELL_INTERFACE cell, REBUILD g_eq using the cell's actual
    // (rho, u) and T = T_ref, then overwrite g_arr.  f is left
    // untouched, so mass and momentum are conserved.  Any T anomaly
    // injected by the ABB / replenish_g / streaming asymmetries at
    // the interface is healed within one step instead of accumulating
    // over thousands of steps until a single cell blows up.
    //
    // Default 0 (legacy: interface T evolves freely from g moments).
    // Set lbm.fslbm_interface_isothermal = 1 to activate.
    bool m_fslbm_interface_isothermal = false;

    // Local-density blending for the Step 1b ABB gas-side target ρ_G.
    //
    // Compressible-FSLBM mass leak (June 2026): with ρ_G = ρ_ref and
    // local cells driven to ρ_iv > ρ_ref by the impeller-driven flow,
    // each INTERFACE-GAS link injects ≈ 2·w_q·(ρ_G − ρ_iv) of mass,
    // which is systematically negative in the entire surface band.
    // Summed over ~88k cells × ~10 gas links per step ≈ −2-3 mass/step.
    //
    // The Donath/Körner FSLBM hides this because they assume
    // incompressibility (ρ_iv ≈ ρ_ref).  In the compressible regime,
    // we blend the ABB target toward the local cell density:
    //
    //   ρ_G = (1 − β) · ρ_ref + β · ρ_iv_smoothed   (then add Δρ_laplace)
    //
    // where ρ_iv_smoothed is a 6-face-neighbour-average of CELL_LIQUID
    // and CELL_INTERFACE densities (skipping GAS/SOLID).  β=0 is the
    // legacy ρ_G = ρ_ref.  β=1 is fully local — the existing code
    // comment in Step 1b warns this creates a positive feedback loop;
    // empirically β ≤ 0.3 is safe.  Recommended for compressible runs:
    // β = 0.2.  Set lbm.fslbm_abb_local_rho_blend in the input file.
    amrex::Real m_fslbm_abb_local_rho_blend = 0.0;

    // Per-cell phi-correction (Variant D, June 2026; corrected formula).
    //
    // Closes the FSLBM compressible-regime mass leak by routing only
    // the GAS-link ABB-induced f-mass change at INTERFACE cells back
    // into the phi field.  Per-cell algorithm at each FSLBM step:
    //   m_pre[iv]   = Σ_q f[iv][q]                       (Step 0c snapshot, comp 0)
    //   Δm_full[iv] = Σ_{q : nbr(q) ∈ LIQ/IFC}            (Step 0c, comp 1)
    //                    ( f_pre[nbr(q)][bq] − f_pre[iv][q] )
    //   ... step 1a / 1b / 2 / 3 / 4 run ...
    //   m_post[iv]  = Σ_q f[iv][q]                       (after Step 3)
    //   Δm_ABB      = (m_post − m_pre) − Δm_full
    //   phi[iv] -= Δm_ABB / max(rho_post, ε·ρ_ref)
    //
    // Δm_full is the un-weighted (S_q ≡ 1) bulk-fluid streaming flux —
    // i.e. what the post-stream m_post would be if there were no ABB.
    // Subtracting it from (m_post − m_pre) leaves precisely the GAS-link
    // ABB residual, which is the actual leak.  Solid bounce-back at iv
    // contributes 0 to (m_post − m_pre) so it needs no separate term.
    //
    // f and g moments are NOT modified, so velocity, temperature, and
    // stress are untouched.  Only phi (the volume-fraction ledger) is
    // adjusted.  Schwarzmeier 2023 §3.3 proves no ABB variant balances
    // forces at moving interfaces, so the leak is fundamental; this
    // closes the bookkeeping at the post-processing level.
    //
    // The previous (mass_flux[0]-based) formulation was incorrect: it
    // used Step 2's Körner-weighted Δm_intended in place of Δm_full,
    // which leaks (1 − S_q)·flux at every IFC-IFC link into "unintended"
    // and destabilises the interface band (run 14397288: 6× worse leak
    // than the no-correction baseline 14348303).
    //
    // Default 0 (legacy: leak observable as M_tot drift).  Set
    // lbm.fslbm_abb_mass_correction = 1 to activate.
    bool m_fslbm_abb_mass_correction = false;

    // Per-cell pre-FSLBM-step buffer used by the phi correction above.
    // Allocated 2 components, no ghost.  Filled at the very top of
    // fslbm_advance_surface for INTERFACE cells; read by the correction
    // kernel just after Step 4.
    //   comp 0: m_pre   = Σ_q f_pre[iv][q]
    //   comp 1: Δm_full = un-weighted Σ over LIQ/IFC links
    amrex::Vector<amrex::MultiFab> m_pre_fslbm_mass;

    // Variant E: global f₀ mass clamp (June 2026).
    //
    // At the end of each fslbm_advance_surface call, compute the current
    // total fluid mass
    //   M_current = Σ_{LIQUID} ρ_iv + Σ_{INTERFACE} φ_iv · ρ_iv
    // (same formula the [mass_diag] line prints as M_tot), compare with
    // the target captured at the first call, and inject the per-cell
    // deficit
    //   ε_cell = (M_target − M_current) / N_ifc
    // into f[iv][0] at every CELL_INTERFACE cell.  Because c₀ = 0, the
    // injection preserves momentum, stress, heat-flux contribution, and
    // every other moment of f — only ρ shifts.  Per-pass mass restoration
    // is ε_cell · Σ_{IFC} φ_iv ≈ 0.5·deficit (under-relaxed by <φ>_ifc),
    // so 99% of any transient gap closes within ~10 steps; the slight
    // under-relaxation is intentionally stabilising.
    //
    // Default 0.  Set lbm.fslbm_global_mass_clamp = 1 to activate.
    // Mutually compatible with Variant D, though enabling both is
    // pointless (E does the bookkeeping correction globally).
    bool m_fslbm_global_mass_clamp = false;

    // Interval (in steps) between clamp passes.  Each pass costs ~one
    // 3-component MultiFab::sum() + one small ParallelFor over
    // CELL_INTERFACE cells + one FillBoundary on m_f (~150 µs/step on
    // 180³ H100, or ~3% of step time; worse on multi-GPU where the
    // reduction crosses the network).  Default 10 keeps the deficit
    // bounded at ~10× per-step leak (~10 mass units = 3×10⁻⁶ relative
    // to M_target ~3.4e6) which is still well below noise.  Set to 1
    // to clamp every step.  Ignored when m_fslbm_global_mass_clamp = 0.
    int m_fslbm_global_mass_clamp_interval = 10;

    // Target mass captured at first call to fslbm_advance_surface.  Stays
    // < 0 until then; on first call M_target = M_current_at_step_0.
    // NOTE: not currently checkpointed — restarts re-capture from the
    // first post-restart step (so M_target locks in any leak that
    // occurred up to the checkpoint).  Acceptable for cold-start runs;
    // restart support is a few-line addition in read_checkpoint_file
    // if/when needed.
    amrex::Real m_fslbm_mass_target = amrex::Real(-1.0);

    // Moving body parameters
    bool m_body_is_moving = false;
    amrex::GpuArray<amrex::Real, 3> m_body_velocity = {0.0, 0.0, 0.0};
    // m_body_angular_velocity is the CURRENT (possibly ramped) angular
    // velocity used by every body kernel.  m_body_angular_velocity_target
    // is the steady-state value read from the input file
    // (body.angular_velocity).  m_body_angular_velocity_ramp_steps
    // controls a linear ramp at startup (step 0 -> 0; ramp_steps -> target):
    //   0  : no ramp; impeller spins at the target rate from step 0 (legacy)
    //   N>0: omega(step) = (step / N) * target   for step <= N
    //        omega(step) =                target for step >  N
    // A long ramp lets the wall boundary layers develop without seeing the
    // full Ma=0.067 jet at startup, which produced CELL_LIQUID rho
    // excursions to ~300 at the baffle wall in run 14310435 (from-zero).
    amrex::GpuArray<amrex::Real, 3> m_body_angular_velocity = {0.0, 0.0, 0.0};
    amrex::GpuArray<amrex::Real, 3> m_body_angular_velocity_target = {0.0, 0.0, 0.0};
    int m_body_angular_velocity_ramp_steps = 0;
    amrex::GpuArray<amrex::Real, 3> m_body_center = {0.0, 0.0, 0.0};
    amrex::Real m_body_rotation_angle = 0.0;  // accumulated rotation
    std::string m_body_geom_type = "";

    // Body is isothermal.
    bool m_bodyIsIsothermal = false;

    // Entire fluid is isothermal (force T = body_temperature everywhere).
    bool m_fluidIsIsothermal = false;

    // Body temperature.
    amrex::Real m_bodyTemperature = 1.0 / 3.0;

    // Adiabatic exponent.
    amrex::Real m_adiabaticExponent = 5.0 / 3.0;

    // Universal Gas constant, LB unit. Never change.
    const amrex::Real m_R_u = constants::AIR_MOLAR_MASS;

    // Mean molecular mass, Safe default air. Unit kg/kmol.
    amrex::Real m_m_bar = constants::AIR_MOLAR_MASS;

    // Set the actual reference speed of sound
    amrex::Real m_speedOfSound_Ref = std::sqrt(
        m_adiabaticExponent * (m_R_u / m_m_bar) * m_initialTemperature);

    // ---------------------------------------------------------------
    // Two-step catalytic reaction parameters
    // Mechanism:  S + C <--k_f/k_r--> I --k_p--> P + C
    //   component 0 = S  (substrate)
    //   component 1 = C  (catalyst)
    //   component 2 = I  (intermediate)
    //   component 3 = P  (product)
    // ---------------------------------------------------------------
    bool        m_enable_reactions    = false;
    amrex::Real m_rxn_k_forward       = 0.0;  // bimolecular  [LB_rho^-1 step^-1]
    amrex::Real m_rxn_k_reverse       = 0.0;  // unimolecular [step^-1]
    amrex::Real m_rxn_k_product       = 0.0;  // unimolecular [step^-1]

    // Timed catalyst injection
    bool        m_cat_inject_done     = false;
    int         m_cat_inject_step     = -1;   // LB step for catalyst release
    amrex::Real m_cat_inject_density  = 0.0;  // LB rho assigned in box
    amrex::RealVect m_cat_inject_box_lo{AMREX_D_DECL(0.0, 0.0, 0.0)};
    amrex::RealVect m_cat_inject_box_hi{AMREX_D_DECL(0.0, 0.0, 0.0)};

    // Species statistics output
    std::string   m_species_stats_file{"species_stats.csv"};
    std::fstream  m_species_stats_stream;

    // ---------------------------------------------------------------
    // Lagrangian bubble physics (two-phase kLa mass transfer)
    // ---------------------------------------------------------------
    bool             m_enable_bubbles      = false;
    BubbleManager    m_bubbles;
    BubbleParams     m_bubble_params;

    /// Concentration reference: 1 LB_rho unit ≡ m_bubble_o2_C_ref [mol/m³]
    /// Used to convert k_La O2 source from mol/(m³·s) → LB_rho/step.
    amrex::Real      m_bubble_o2_C_ref     = 1.0;

    /// Internal step counter for coalescence and stats intervals
    int              m_bubble_step_counter = 0;
};
} // namespace lbm
#endif
