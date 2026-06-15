#ifndef BUBBLE_H
#define BUBBLE_H

#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <fstream>
#include <string>

namespace lbm {

namespace BubbleIdx {
    enum Real {
        VX       = 0,
        VY       = 1,
        VZ       = 2,
        DIAMETER = 3,
        N_O2     = 4,
        AX       = 5,
        AY       = 6,
        AZ       = 7,
        BREAKUP_COOLDOWN = 8,
        DN_I     = 9,
        EPS_CACHED = 10,
        NReal    = 11
    };
}

using BubbleContainer = amrex::ParticleContainer<BubbleIdx::NReal, 0, 0, 0, amrex::ManagedArenaAllocator>;
using BubbleParticle  = typename BubbleContainer::ParticleType;

struct BubbleParams
{
    amrex::Real dx_phys  = 1.0;
    amrex::Real dt_phys  = 1.0;
    amrex::Real nu_lb    = 0.01;

    int n_sparger_points = 7;
    amrex::Vector<amrex::Real> sparger_x;
    amrex::Vector<amrex::Real> sparger_y;
    amrex::Vector<amrex::Real> sparger_z;
    amrex::Real flow_rate       = 6.67e-6;
    amrex::Real d0              = 1.1e-3;

    amrex::Real O2_molar_volume = 0.0224;
    amrex::Real O2_solubility   = 0.032;
    amrex::Real O2_init_conc    = 44.6;
    amrex::Real kL_coeff        = 0.301;
    amrex::Real D_O2            = 2.0e-9;

    amrex::Real rho_fluid       = 1000.0;
    amrex::Real nu_fluid        = 1.0e-6;
    amrex::Real surf_tension    = 0.072;
    amrex::Real P_atm           = 101325.;
    amrex::Real g_grav          = 9.81;

    amrex::Real free_surface_z  = 1.0e30;

    bool        enable_coalescence     = true;
    amrex::Real Re_a_crit             = 40.0;
    amrex::Real coal_start_time       = 2.0;
    int         coal_interval         = 5;

    bool        enable_breakup        = true;
    int         breakup_model         = 0;
    amrex::Real min_diameter          = 1.0e-4;
    amrex::Real eps_max               = 5000.0;

    amrex::Real C_ref               = 100.0;

    int         stats_int             = 200;
    std::string stats_file            = "bubble_stats.csv";
};

class BubbleManager
{
public:
    BubbleManager() = default;
    ~BubbleManager() { close_stats_file(); }

    static void read_params(BubbleParams& p);

    void initialize(
        const amrex::Geometry&            geom,
        const amrex::BoxArray&            ba,
        const amrex::DistributionMapping& dm,
        const BubbleParams&               params,
        bool                              append_stats = false);

    void advance(
        amrex::Real                dt,
        const amrex::MultiFab&     macrodata,
        const amrex::MultiFab&     derived,
        const amrex::MultiFab&     o2_conc_mf,
        const amrex::Geometry&     geom,
        amrex::MultiFab&           fluid_force_mf,
        amrex::MultiFab&           o2_src_mf,
        amrex::Real                phys_time,
        const amrex::MultiFab*     phi_mf       = nullptr,
        // Cell-type field used for solid-wall blocking.  We accept the FSLBM
        // m_cell_type iMultiFab (values CELL_LIQUID / CELL_INTERFACE / CELL_GAS
        // / CELL_SOLID) rather than the raw m_is_fluid mask: the raw mask sets
        // both GAS and SOLID to 0, so testing IS_FLUID==0 inside advance()
        // wrongly treated bubbles rising into the gas headspace as colliding
        // with a wall.  Testing == CELL_SOLID isolates real wall contact.
        const amrex::iMultiFab*    cell_type_mf = nullptr);

    void inject_bubbles(amrex::Real dt_phys);
    void do_coalescence(amrex::Real phys_time);
    void write_stats(int step, amrex::Real phys_time);

    BubbleContainer&       container()       { return m_container; }
    const BubbleContainer& container() const { return m_container; }

    bool is_initialized() const { return m_initialized; }

    void Checkpoint(const std::string& dir, const std::string& name) const;
    void Restart(const std::string& dir, const std::string& name);

private:
    void do_breakup(
        const amrex::MultiFab& derived,
        const amrex::Geometry& geom);

    void open_stats_file(bool append = false);
    void close_stats_file();

    BubbleParams m_params;
    BubbleContainer m_container;
    bool m_initialized = false;
    bool m_particles_ever_injected = false;
    amrex::Vector<amrex::Real> m_injection_residuals;
    std::ofstream m_stats_stream;
};

} // namespace lbm

#endif
