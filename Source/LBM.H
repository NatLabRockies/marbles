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

    // write checkpoint file to disk
    void write_checkpoint_file() const;

    // read checkpoint file from disk
    void read_checkpoint_file();

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

    // save streaming data to plot file
    bool m_save_streaming = true;

    // save derived data to plot file
    bool m_save_derived = true;

    // compute and output forces
    bool m_compute_forces = false;

    // clamp negative component densities before collision
    bool m_clamp_component_densities = false;

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

    // Moving body parameters
    bool m_body_is_moving = false;
    amrex::GpuArray<amrex::Real, 3> m_body_velocity = {0.0, 0.0, 0.0};
    amrex::GpuArray<amrex::Real, 3> m_body_angular_velocity = {0.0, 0.0, 0.0};
    amrex::GpuArray<amrex::Real, 3> m_body_center = {0.0, 0.0, 0.0};
    amrex::Real m_body_rotation_angle = 0.0;  // accumulated rotation
    std::string m_body_geom_type = "";

    // Body is isothermal.
    bool m_bodyIsIsothermal = false;

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
};
} // namespace lbm
#endif
