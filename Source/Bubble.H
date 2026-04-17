/**
 * @file Bubble.H
 * @brief Lagrangian bubble tracking for kLa two-phase mass transfer.
 *
 * Implements Newton's law bubble dynamics (Thomas et al. 2021, CES 237,
 * 116538) with:
 *   - Velocity Verlet time integration
 *   - Forces: gravity/buoyancy, Schiller-Naumann drag, added mass,
 *     pressure gradient
 *   - Two-way coupling (F_a + F_D fed back, not pressure gradient)
 *   - Trilinear interpolation of fluid velocity, pressure, epsilon
 *   - Hinze breakup (triangle daughter distribution, mode = 0.2)
 *   - Boshenyatov coalescence (Re_a > 40 criterion)
 *   - Kawase k_L mass transfer (Eq. 17 of Thomas et al.)
 *   - Sparger injection BC
 *   - Boyle's law diameter correction
 *
 * Coordinate convention
 * ---------------------
 * Bubble POSITIONS are stored in the AMReX geometry coordinate frame
 * (lattice-cell units, matching prob_lo/prob_hi from the input file).
 * Bubble VELOCITIES and all force/transport calculations are in SI
 * (m/s, N, mol/s etc.) and converted using dx_phys [m/cell] and
 * dt_phys [s/step] supplied in BubbleParams.
 */
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

// ---------------------------------------------------------------------------
// Particle real-data component indices (stored in AoS struct)
// ---------------------------------------------------------------------------
namespace BubbleIdx {
    enum Real {
        VX       = 0,  ///< x-velocity in SI [m/s]
        VY       = 1,  ///< y-velocity in SI [m/s]
        VZ       = 2,  ///< z-velocity in SI [m/s]
        DIAMETER = 3,  ///< bubble diameter in SI [m]
        N_O2     = 4,  ///< moles of O2 in bubble [mol]
        AX       = 5,  ///< x-acceleration from previous step [m/s²] — for vel. Verlet
        AY       = 6,  ///< y-acceleration from previous step [m/s²]
        AZ       = 7,  ///< z-acceleration from previous step [m/s²]
        NReal    = 8
    };
} // namespace BubbleIdx

// Traditional AoS ParticleContainer:
//   NStructReal = 8  (vx, vy, vz, diameter, n_O2, ax, ay, az)
//   NStructInt  = 0
using BubbleContainer = amrex::ParticleContainer<BubbleIdx::NReal, 0>;
using BubbleParticle  = typename BubbleContainer::ParticleType;

// ---------------------------------------------------------------------------
// Runtime parameters — set from ParmParse sections "sparger" and "bubble"
// ---------------------------------------------------------------------------
struct BubbleParams
{
    // Unit conversions (set by LBM after reading its own params)
    amrex::Real dx_phys  = 1.0;   ///< metres per lattice cell
    amrex::Real dt_phys  = 1.0;   ///< seconds per lattice step
    amrex::Real nu_lb    = 0.01;  ///< LB kinematic viscosity (dx²/step)

    // Sparger geometry
    int n_sparger_points = 7;
    amrex::Vector<amrex::Real> sparger_x; ///< positions in LB cells
    amrex::Vector<amrex::Real> sparger_y;
    amrex::Vector<amrex::Real> sparger_z;
    amrex::Real flow_rate       = 6.67e-6; ///< total gas flow [m³/s]
    amrex::Real d0              = 1.1e-3;  ///< initial bubble diameter [m]

    // O2 / mass transfer
    amrex::Real O2_molar_volume = 0.0224;  ///< m³/mol at STP
    amrex::Real O2_solubility   = 0.032;   ///< dimensionless Bunsen coefficient
    amrex::Real O2_init_conc    = 44.6;    ///< mol/m³ — pure O2 at STP
    amrex::Real kL_coeff        = 0.301;   ///< C in Kawase k_L formula
    amrex::Real D_O2            = 2.0e-9;  ///< O2 diffusivity in fluid [m²/s]

    // Fluid properties
    amrex::Real rho_fluid       = 1000.0;  ///< kg/m³
    amrex::Real nu_fluid        = 1.0e-6;  ///< m²/s (SI kinematic viscosity)
    amrex::Real surf_tension    = 0.072;   ///< N/m air-water
    amrex::Real P_atm           = 101325.; ///< Pa
    amrex::Real g_grav          = 9.81;    ///< m/s²

    // Free surface
    amrex::Real free_surface_z  = 1.0e30;  ///< LB cells — bubbles above exit

    // Coalescence
    bool        enable_coalescence     = true;
    amrex::Real Re_a_crit             = 40.0;
    amrex::Real coal_start_time       = 2.0;  ///< s — delay after startup
    int         coal_interval         = 5;    ///< steps between checks

    // Breakup
    bool        enable_breakup        = true;
    int         breakup_model         = 0;  ///< 0=triangle, 1=equal, 2=random

    // Concentration reference: 1 LB_rho ≡ C_ref mol/m³ dissolved O2
    // Set by LBM after reading lbm.O2_concentration_reference.
    amrex::Real C_ref               = 100.0; ///< mol/m³ per LB_rho unit

    // Output
    int         stats_int             = 200;
    std::string stats_file            = "bubble_stats.csv";
};

// ---------------------------------------------------------------------------
// BubbleManager: owns BubbleContainer + all bubble physics
// ---------------------------------------------------------------------------
class BubbleManager
{
public:
    BubbleManager() = default;
    ~BubbleManager() { close_stats_file(); }

    /// Read all "sparger.*" and "bubble.*" ParmParse entries into p.
    static void read_params(BubbleParams& p);

    /// Initialise particle container; call once after AMReX geometry is known.
    void initialize(
        const amrex::Geometry&            geom,
        const amrex::BoxArray&            ba,
        const amrex::DistributionMapping& dm,
        const BubbleParams&               params);

    /**
     * @brief Main per-step update.
     *
     * Performs (in order):
     *  1. Boyle's law diameter correction
     *  2. Force computation (gravity/buoyancy, drag, added mass, pressure
     *     gradient) by trilinear interpolation from macrodata/derived
     *  3. Velocity Verlet position + velocity update
     *  4. Remove bubbles that have crossed the free surface
     *  5. Do breakup check at new positions
     *  6. Deposit two-way coupling force (F_a + F_D) onto fluid_force_mf
     *  7. Compute k_L, Dn_i and deposit O2 source onto o2_src_mf
     *  8. Shrink bubble volumes by transferred O2
     *
     * @param dt          LB time step (steps — use 1 per level advance)
     * @param macrodata   LBM macrodata MultiFab (rho, u, P_ij, ...)
     * @param derived     LBM derived MultiFab (vorticity, epsilon)
     * @param o2_conc_mf  Dissolved O2 concentration MultiFab (1 comp, LB rho units)
     * @param geom        AMReX Geometry (LB coordinate frame)
     * @param fluid_force_mf (output) 3-comp force density [LB force/volume]
     * @param o2_src_mf      (output) O2 source [mol/(m³·s)] per cell
     * @param phys_time   Physical simulation time in seconds
     * @param phi_mf      Phase-field Φ ∈ [0,1] (liquid=1, gas=0); when provided,
     *                    bubbles are removed when Φ < 0.5 at their location
     *                    (Chiu & Lin free surface).  Pass nullptr to fall back
     *                    to the fixed free_surface_z threshold.
     */
    void advance(
        amrex::Real                dt,
        const amrex::MultiFab&     macrodata,
        const amrex::MultiFab&     derived,
        const amrex::MultiFab&     o2_conc_mf,
        const amrex::Geometry&     geom,
        amrex::MultiFab&           fluid_force_mf,
        amrex::MultiFab&           o2_src_mf,
        amrex::Real                phys_time,
        const amrex::MultiFab*     phi_mf = nullptr);

    /// Sparger injection: create new bubbles at a rate matching flow_rate.
    /// Call every step before advance().
    void inject_bubbles(amrex::Real dt_phys);

    /// Bubble–bubble coalescence check (call every coal_interval steps).
    void do_coalescence(amrex::Real phys_time);

    /// Write per-step bubble statistics CSV row.
    void write_stats(int step, amrex::Real phys_time);

    BubbleContainer&       container()       { return m_container; }
    const BubbleContainer& container() const { return m_container; }

    bool is_initialized() const { return m_initialized; }

private:
    // -----------------------------------------------------------------------
    // Helper methods called from advance()
    // -----------------------------------------------------------------------

    /// Apply Boyle's law diameter scaling (depth-pressure correction).
    void apply_boyles_law(const amrex::Geometry& geom);

    /**
     * Compute all forces on each bubble by trilinear interpolation from
     * macrodata and derived fields.  Results stored in m_forces.
     */
    void compute_forces(
        const amrex::MultiFab& macrodata,
        const amrex::MultiFab& derived,
        const amrex::Geometry& geom);

    /**
     * Detect and execute bubble breakup.
     * Bubbles with D > De (Hinze scale) are split into two daughters.
     * New bubbles are appended to the container.
     * @param derived  LBM derived MultiFab containing epsilon at EPSILON_IDX.
     * @param geom     AMReX Geometry for trilinear interpolation.
     */
    void do_breakup(
        const amrex::MultiFab& derived,
        const amrex::Geometry& geom);

    /**
     * Remove bubbles whose centre has risen above the free surface.
     */
    /// @param phi_mf  If non-null, use Φ < 0.5 criterion; else use free_surface_z.
    void remove_exited_bubbles(const amrex::Geometry& geom,
                               const amrex::MultiFab* phi_mf = nullptr);

    /**
     * Deposit two-way coupling reaction force (F_a + F_D, SI N) converted
     * to LB force density onto fluid_force_mf.
     */
    void deposit_fluid_forces(
        amrex::MultiFab&       fluid_force_mf,
        const amrex::Geometry& geom);

    /**
     * Compute per-bubble mass transfer rates dn_i [mol/s],
     * deposit as O2 source density [mol/(m³·s)] onto o2_src_mf,
     * and shrink bubble volumes accordingly.
     *
     * @param o2_src_mf   (output) O2 source rate per cell [mol/(m³·s)]
     * @param o2_conc_mf  Dissolved O2 in LB density units (comp 0)
     * @param derived     LBM derived field (epsilon at EPSILON_IDX)
     * @param geom        Geometry for cell volume
     */
    void deposit_o2_sources(
        amrex::MultiFab&       o2_src_mf,
        const amrex::MultiFab& o2_conc_mf,
        const amrex::MultiFab& derived,
        const amrex::Geometry& geom);

    // -----------------------------------------------------------------------
    // Interpolation helpers
    // -----------------------------------------------------------------------

    /**
     * Trilinear interpolation of a scalar field at position (px, py, pz)
     * in LB cell coordinates.
     */
    static amrex::Real trilinear_interp(
        const amrex::MultiFab& mf,
        int                    comp,
        const amrex::Geometry& geom,
        amrex::Real px, amrex::Real py, amrex::Real pz);

    // -----------------------------------------------------------------------
    // Statistics file
    // -----------------------------------------------------------------------
    void open_stats_file();
    void close_stats_file();

    // -----------------------------------------------------------------------
    // Data members
    // -----------------------------------------------------------------------
    BubbleContainer m_container;
    BubbleParams    m_params;

    // Per-bubble force cache.  One entry per particle in traversal order.
    // Filled each step by compute_forces(), consumed by advance().
    struct BubbleForces {
        amrex::Real fx_a, fy_a, fz_a;  ///< added-mass force [N, SI]
        amrex::Real fx_D, fy_D, fz_D;  ///< drag force [N, SI]
        amrex::Real dn_i;               ///< O2 transfer rate [mol/s]
    };
    amrex::Vector<BubbleForces> m_forces;

    // Fractional injection counter (accumulated between integer events)
    amrex::Real m_injection_residual = 0.0;

    bool       m_initialized   = false;
    std::fstream m_stats_stream;
};

} // namespace lbm
#endif // BUBBLE_H
