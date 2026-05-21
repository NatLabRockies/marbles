/**
 * @file Bubble.cpp
 * @brief Implementation of Lagrangian bubble tracking for kLa mass transfer.
 *
 * Reference:
 *   Thomas et al. (2021), CES 237, 116538. DOI: 10.1016/j.ces.2021.116538
 *   Rettinger & Rüde (2018), Comput Fluids 172 — point-particle LBM–DEM coupling (cited at Eq. 9 of Thomas et al.)
 *   Tenneti, Garg & Subramaniam (2011), Int J Multiph Flow 37(9), 1072–1092 — drag law (Eq. 18)
 *   Kawase et al. (1992), Biotechnol. Bioeng. 39(11) — k_L formula
 *   Boshenyatov (2012) — coalescence criterion Re_a > 40
 *   Hinze (1955) — breakup equilibrium diameter
 */
#include "Bubble.H"
#include "Constants.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <cmath>
#include <algorithm>
#include <random>
#include <sstream>
#include <iomanip>

namespace lbm {

// ============================================================================
// Drag coefficient for a sphere in a suspension.
//
// Implements the Tenneti, Garg & Subramaniam (2011) drag law, which is the
// correlation used by Rettinger & Rüde (2018) (their Table 1, ref [38]),
// as cited by Thomas et al. (2021) Eq. 9.
//
// Reference:
//   Tenneti, Garg & Subramaniam (2011), Int J Multiph Flow 37(9), 1072–1092.
//   DOI: 10.1016/j.ijmultiphaseflow.2011.05.010
//
// The force on particle i is:
//
//   F_d = 3π μ d (u_f - u_p) * eps_f * F(Re_p, eps_p)
//
// where F is the dimensionless drag (Tenneti et al. Eq. 18):
//
//   F = F_0 + eps_p * F_1 + F_2
//
//   F_0  =  (1 + 0.15*Re^0.687) / eps_f^3
//
//   F_1  =  10 / eps_f^2           (dense-packing correction)
//
//   F_2  =  (eps_p / eps_f^3) * 0.413*Re/24
//           * (eps_f^{-1} + 3*eps_p*eps_f + 8.4*Re^{-0.343})
//           / (1 + 10^{3*eps_p} * Re^{-(1+4*eps_p)/2})
//
// This function returns C_D such that:
//   F_d  =  0.5 * C_D * rho_f * A_b * |u_rel|^2
//        =  0.5 * C_D * rho_f * (π r²) * |u_rel|^2
//
// Conversion from Tenneti's F to C_D:
//   F_d(Tenneti) = 3π μ d |u_rel| eps_f F
//              = 0.5 * (24/Re) * (1/eps_f^3) * rho_f * A_b * |u_rel|^2 * eps_f * eps_f^3 * F/F_SN
//   →  C_D = (24 / Re) * eps_f^3 * F          [exactly, by construction]
//       where F is the full Tenneti expression above.
//
// Dilute limit (eps_p → 0, eps_f → 1):
//   F_0 → 1 + 0.15*Re^0.687,  F_1=F_2 → 0
//   C_D → (24/Re)*(1+0.15*Re^0.687)  — recovers Schiller-Naumann exactly.
//
// Parameters:
//   Re    — bubble Reynolds number = |u_rel| * d / nu_f
//   eps_p — local solid (bubble) volume fraction in the cell [0,1)
// ============================================================================
static amrex::Real bubble_drag_cd(amrex::Real Re, amrex::Real eps_p)
{
    if (Re < 1.0e-10) { return 0.0; }

    // Guard: clamp eps_p to physically meaningful range
    eps_p = std::max(0.0, std::min(eps_p, 0.99));
    const amrex::Real eps_f = 1.0 - eps_p;
    const amrex::Real eps_f3 = eps_f * eps_f * eps_f;

    // F_0: Schiller-Naumann corrected by 1/eps_f^3
    const amrex::Real F0 = (1.0 + 0.15 * std::pow(Re, 0.687)) / eps_f3;

    // F_1: dense-packing term
    const amrex::Real F1 = 10.0 * eps_p / (eps_f * eps_f);

    // F_2: cross term (Re × eps_p), Tenneti et al. Eq. 18
    amrex::Real F2 = 0.0;
    if (eps_p > 1.0e-8) {
        const amrex::Real numer = eps_f / eps_f3 + 3.0 * eps_p * eps_f + 8.4 * std::pow(Re, -0.343);
        const amrex::Real denom = 1.0 + std::pow(10.0, 3.0 * eps_p) *
                                        std::pow(Re, -(1.0 + 4.0 * eps_p) / 2.0);
        F2 = (eps_p / eps_f3) * (0.413 * Re / 24.0) * (numer / denom);
    }

    // F = F0 + F1 + F2  (Tenneti Eq. 18)
    const amrex::Real F = F0 + F1 + F2;

    // Convert to C_D convention: F_d = 0.5 * C_D * rho_f * A_b * |u_rel|^2
    // F_d(Tenneti) = 3π μ d |u_rel| eps_f F = 0.5 * C_D * rho_f * π*(d/2)^2 * |u_rel|^2
    // → C_D = 24*eps_f*F / Re
    const amrex::Real Cd = 24.0 * eps_f * F / Re;

    // High-Re cap (Tenneti correlation is validated for Re < ~300)
    return (Re < 1000.0) ? Cd : std::min(Cd, 0.44);
}

// ============================================================================
// Utility: trilinear interpolation in a MultiFab at LB-cell position (px,py,pz)
// ============================================================================
amrex::Real BubbleManager::trilinear_interp(
    const amrex::MultiFab& mf,
    int                    comp,
    const amrex::Geometry& geom,
    amrex::Real px, amrex::Real py, amrex::Real pz)
{
    // Cell-centre coordinates in LB cells: cell (i,j,k) has centre at (i+0.5, j+0.5, k+0.5)
    const amrex::Real* plo    = geom.ProbLo();
    const amrex::Real* dx_arr = geom.CellSize();

    // Fractional index
    amrex::Real fi = (px - plo[0]) / dx_arr[0] - 0.5;
    amrex::Real fj = (py - plo[1]) / dx_arr[1] - 0.5;
    amrex::Real fk = (pz - plo[2]) / dx_arr[2] - 0.5;

    int i0 = static_cast<int>(std::floor(fi));
    int j0 = static_cast<int>(std::floor(fj));
    int k0 = static_cast<int>(std::floor(fk));
    int i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1;

    amrex::Real ti = fi - i0;   // [0,1] weight toward i1
    amrex::Real tj = fj - j0;
    amrex::Real tk = fk - k0;

    const amrex::Box& domain = geom.Domain();
    // Clamp to valid domain
    auto clamp = [](int v, int lo, int hi) { return std::max(lo, std::min(hi, v)); };
    i0 = clamp(i0, domain.smallEnd(0), domain.bigEnd(0));
    i1 = clamp(i1, domain.smallEnd(0), domain.bigEnd(0));
    j0 = clamp(j0, domain.smallEnd(1), domain.bigEnd(1));
    j1 = clamp(j1, domain.smallEnd(1), domain.bigEnd(1));
    k0 = clamp(k0, domain.smallEnd(2), domain.bigEnd(2));
    k1 = clamp(k1, domain.smallEnd(2), domain.bigEnd(2));

    // Find which FArrayBox owns cell (i0,j0,k0)
    // For simplicity, iterate over tiles on the MultiFab and find the box.
    amrex::Real val = 0.0;
    bool found = false;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        if (!bx.contains(amrex::IntVect(i0, j0, k0))) { continue; }
        const auto& arr = mf.const_array(mfi);

        auto get = [&](int ii, int jj, int kk) {
            amrex::Real v = arr(amrex::IntVect(clamp(ii, bx.smallEnd(0), bx.bigEnd(0)),
                                     clamp(jj, bx.smallEnd(1), bx.bigEnd(1)),
                                     clamp(kk, bx.smallEnd(2), bx.bigEnd(2))), comp);
            // EB/solid cells may contain NaN; treat as zero for interpolation purposes
            return std::isfinite(v) ? v : amrex::Real(0.0);
        };

        val =
            (1-ti)*(1-tj)*(1-tk)*get(i0,j0,k0) +
               ti *(1-tj)*(1-tk)*get(i1,j0,k0) +
            (1-ti)*   tj *(1-tk)*get(i0,j1,k0) +
               ti *   tj *(1-tk)*get(i1,j1,k0) +
            (1-ti)*(1-tj)*   tk *get(i0,j0,k1) +
               ti *(1-tj)*   tk *get(i1,j0,k1) +
            (1-ti)*   tj *   tk *get(i0,j1,k1) +
               ti *   tj *   tk *get(i1,j1,k1);
        found = true;
        break;
    }
    if (!found) {
        // Particle is outside valid domain — return zero
        val = 0.0;
    }
    return val;
}

// ============================================================================
// read_params
// ============================================================================
void BubbleManager::read_params(BubbleParams& p)
{
    {
        amrex::ParmParse pp("sparger");
        pp.query("n_points",   p.n_sparger_points);
        pp.query("flow_rate",  p.flow_rate);
        pp.query("bubble_diameter", p.d0);

        // Sparger positions — expected as flat list x0 y0 z0 x1 y1 z1 ...
        amrex::Vector<amrex::Real> pos_flat;
        if (pp.queryarr("positions", pos_flat)) {
            p.sparger_x.resize(p.n_sparger_points);
            p.sparger_y.resize(p.n_sparger_points);
            p.sparger_z.resize(p.n_sparger_points);
            for (int i = 0; i < p.n_sparger_points; ++i) {
                p.sparger_x[i] = pos_flat[3*i + 0];
                p.sparger_y[i] = pos_flat[3*i + 1];
                p.sparger_z[i] = pos_flat[3*i + 2];
            }
        } else {
            // Fallback: sparger positions file not provided — will be set to default
            p.sparger_x.assign(p.n_sparger_points, 0.0);
            p.sparger_y.assign(p.n_sparger_points, 0.0);
            p.sparger_z.assign(p.n_sparger_points, 0.0);
        }
    }

    {
        amrex::ParmParse pp("bubble");
        pp.query("O2_molar_volume",  p.O2_molar_volume);
        pp.query("O2_solubility",    p.O2_solubility);
        pp.query("O2_initial_conc",  p.O2_init_conc);
        pp.query("kL_coefficient",   p.kL_coeff);
        pp.query("D_O2",             p.D_O2);
        pp.query("rho_fluid",        p.rho_fluid);
        pp.query("nu_fluid",         p.nu_fluid);
        pp.query("surface_tension",  p.surf_tension);
        pp.query("P_atm",            p.P_atm);
        pp.query("g_grav",           p.g_grav);
        pp.query("free_surface_z",   p.free_surface_z);

        // Coalescence
        int ienb = p.enable_coalescence ? 1 : 0;
        pp.query("enable_coalescence",   ienb); p.enable_coalescence = ienb;
        pp.query("coalescence_Re_crit",  p.Re_a_crit);
        pp.query("coalescence_start_time", p.coal_start_time);
        pp.query("coalescence_interval", p.coal_interval);

        // Breakup
        int ienbk = p.enable_breakup ? 1 : 0;
        pp.query("enable_breakup",       ienbk); p.enable_breakup = ienbk;
        std::string bmodel = "triangle";
        pp.query("breakup_model",        bmodel);
        if      (bmodel == "equal" )   { p.breakup_model = 1; }
        else if (bmodel == "random")   { p.breakup_model = 2; }
        else                           { p.breakup_model = 0; } // triangle
        pp.query("min_diameter", p.min_diameter);
        pp.query("eps_max",      p.eps_max);

        pp.query("stats_int",  p.stats_int);
        pp.query("stats_file", p.stats_file);
    }
}

// ============================================================================
// initialize
// ============================================================================
void BubbleManager::initialize(
    const amrex::Geometry&            geom,
    const amrex::BoxArray&            ba,
    const amrex::DistributionMapping& dm,
    const BubbleParams&               params)
{
    m_params = params;
    m_container.Define(geom, dm, ba);
    m_container.resizeData();        // resizes m_particles to finestLevel()+1 levels
                                     // (Define() only sets GDB; constructors call
                                     //  resizeData() automatically but Define() does not)
    // Per-hole residuals, pre-seeded to 1.0 so the first call always injects.
    m_injection_residuals.assign(m_params.n_sparger_points, 1.0);
    m_initialized = true;
    open_stats_file();
    amrex::Print() << "[BubbleManager] Initialized. n_sparger_points = "
                   << m_params.n_sparger_points << "\n";
}

// ============================================================================
// inject_bubbles
//   Inject bubbles at sparger locations.
//   dt_phys: physical seconds in this step.
//   All injection is staged on rank 0; Redistribute() scatters to correct ranks.
// ============================================================================
void BubbleManager::inject_bubbles(amrex::Real dt_phys)
{
    if (!m_initialized) { return; }
    if (m_params.sparger_x.empty()) { return; }

    const amrex::Real V0 = (amrex::Math::pi<amrex::Real>() / 6.0) *
                           std::pow(m_params.d0, 3.0); // m³ per initial bubble
    // Bubbles per second per hole
    const amrex::Real rate_per_hole =
        (m_params.flow_rate / V0) / m_params.n_sparger_points;

    // Each hole independently accumulates its fractional bubble count so that
    // at low injection rates every hole fires equally rather than only hole 0.
    bool any_inject = false;
    for (int ih = 0; ih < m_params.n_sparger_points; ++ih) {
        m_injection_residuals[ih] += rate_per_hole * dt_phys;
        if (static_cast<int>(m_injection_residuals[ih]) > 0) any_inject = true;
    }

    if (!any_inject) {
        // Still call Redistribute to ensure tile map is populated even on steps
        // where no bubbles are injected (e.g. first step when residual < 1).
        if (!m_particles_ever_injected) { return; }  // m_particles not yet sized
        m_container.Redistribute();
        return;
    }

    // All injection on rank 0; Redistribute() scatters to correct ranks
    if (amrex::ParallelDescriptor::IOProcessor()) {
        auto& ptile = m_container.DefineAndReturnParticleTile(0, 0, 0);
        static std::mt19937 rng(42);
        std::uniform_real_distribution<amrex::Real> jitter(-0.3, 0.3);
        const amrex::Geometry& geom = m_container.Geom(0);
        const amrex::Real* plo  = geom.ProbLo();
        const amrex::Real* phi  = geom.ProbHi();

        for (int ih = 0; ih < m_params.n_sparger_points; ++ih) {
            int n_here = static_cast<int>(m_injection_residuals[ih]);
            m_injection_residuals[ih] -= static_cast<amrex::Real>(n_here);
            for (int ib = 0; ib < n_here; ++ib) {
                BubbleParticle p;
                p.id()  = BubbleParticle::NextID();
                p.cpu() = amrex::ParallelDescriptor::MyProc();
                p.pos(0) = m_params.sparger_x[ih] + jitter(rng);
                p.pos(1) = m_params.sparger_y[ih] + jitter(rng);
                p.pos(2) = m_params.sparger_z[ih] + jitter(rng);
                // Clamp to valid geometry
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    p.pos(d) = std::max(plo[d] + 0.5,
                                        std::min(phi[d] - 0.5, p.pos(d)));
                }
                p.rdata(BubbleIdx::VX)       = 0.0;
                p.rdata(BubbleIdx::VY)       = 0.0;
                p.rdata(BubbleIdx::VZ)       = 0.0;
                p.rdata(BubbleIdx::DIAMETER) = m_params.d0;
                p.rdata(BubbleIdx::N_O2)     = m_params.O2_init_conc * V0;
                p.rdata(BubbleIdx::AX)       = 0.0;
                p.rdata(BubbleIdx::AY)       = 0.0;
                p.rdata(BubbleIdx::AZ)       = 0.0;
                ptile.push_back(p);
            }
        }
    }
    m_container.Redistribute();
    m_particles_ever_injected = true;
}

// ============================================================================
// apply_boyles_law
//   Scale bubble diameter based on hydrostatic pressure at depth h below surface.
//   d / d_atm = (1 + rho*g*h / P_atm)^(-1/3)
// ============================================================================
void BubbleManager::apply_boyles_law(const amrex::Geometry& geom)
{
    const amrex::Real dx      = m_params.dx_phys;
    const amrex::Real rho_f   = m_params.rho_fluid;
    const amrex::Real g       = m_params.g_grav;
    const amrex::Real P_atm   = m_params.P_atm;
    const amrex::Real z_surf  = m_params.free_surface_z;  // LB cells

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& pbox = kv.second;
            auto& aos  = pbox.GetArrayOfStructs();
            for (auto& p : aos()) {
                if (!p.id().is_valid()) { continue; }
                // Depth below free surface in metres
                amrex::Real z_pos = p.pos(2);
                amrex::Real h = (z_surf - z_pos) * dx;  // metres; positive below surface
                if (h <= 0.0) { h = 0.0; }
                amrex::Real d_atm = p.rdata(BubbleIdx::DIAMETER);
                amrex::Real d_new = d_atm * std::pow(1.0 + rho_f * g * h / P_atm, -1.0/3.0);
                p.rdata(BubbleIdx::DIAMETER) = d_new;
            }
        }
    }
}

// ============================================================================
// compute_forces
//   Net buoyancy + drag + added mass forces per bubble.
//   Results go into m_forces (indexed in same traversal order as main loop).
// ============================================================================
void BubbleManager::compute_forces(
    const amrex::MultiFab& macrodata,
    const amrex::MultiFab& derived,
    const amrex::Geometry& geom)
{
    m_forces.clear();

    const amrex::Real dx     = m_params.dx_phys;
    const amrex::Real dt     = m_params.dt_phys;
    const amrex::Real rho_f  = m_params.rho_fluid;
    const amrex::Real nu_f   = m_params.nu_fluid;
    const amrex::Real g      = m_params.g_grav;
    const amrex::Real C_AM   = constants::C_ADDED_MASS;

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& pbox = kv.second;
            auto& aos  = pbox.GetArrayOfStructs();
            for (auto& p : aos()) {
                if (!p.id().is_valid()) {
                    m_forces.push_back({0,0,0, 0,0,0, 0});
                    continue;
                }

                const amrex::Real d  = p.rdata(BubbleIdx::DIAMETER);  // m
                if (d <= 1.0e-5) {
                    // Tiny bubble slipped through — skip forces, mark invalid
                    p.id() = -1;
                    m_forces.push_back({0,0,0, 0,0,0, 0});
                    continue;
                }
                const amrex::Real r  = 0.5 * d;
                const amrex::Real Vb = (4.0/3.0) * amrex::Math::pi<amrex::Real>() * r*r*r;
                const amrex::Real Ab = amrex::Math::pi<amrex::Real>() * r*r;   // cross-section

                // Bubble velocity in SI
                const amrex::Real vbx = p.rdata(BubbleIdx::VX);
                const amrex::Real vby = p.rdata(BubbleIdx::VY);
                const amrex::Real vbz = p.rdata(BubbleIdx::VZ);

                // Interpolate fluid velocity at bubble location (LB units → SI)
                // macrodata comps: RHO=0, VELX=1, VELY=2, VELZ=3
                const amrex::Real ufx_lb = trilinear_interp(macrodata, constants::VELX_IDX,
                                                             geom, p.pos(0), p.pos(1), p.pos(2));
                const amrex::Real ufy_lb = trilinear_interp(macrodata, constants::VELY_IDX,
                                                             geom, p.pos(0), p.pos(1), p.pos(2));
                const amrex::Real ufz_lb = trilinear_interp(macrodata, constants::VELZ_IDX,
                                                             geom, p.pos(0), p.pos(1), p.pos(2));
                const amrex::Real ufx = ufx_lb * dx / dt;  // m/s
                const amrex::Real ufy = ufy_lb * dx / dt;
                const amrex::Real ufz = ufz_lb * dx / dt;

                // Relative velocity: u_rel = v_b - u_f
                const amrex::Real urx = vbx - ufx;
                const amrex::Real ury = vby - ufy;
                const amrex::Real urz = vbz - ufz;
                const amrex::Real ur_mag = std::sqrt(urx*urx + ury*ury + urz*urz);

                // Bubble Reynolds number
                const amrex::Real Re_b = (ur_mag > 0.0) ? (ur_mag * d / nu_f) : 0.0;

                // Local solid volume fraction: bubble volume / cell volume
                // (single bubble per cell approximation; conservative lower bound)
                const amrex::Real dx_phys = m_params.dx_phys;
                const amrex::Real V_cell  = dx_phys * dx_phys * dx_phys;
                const amrex::Real eps_p_local = std::min(Vb / V_cell, 0.99);

                const amrex::Real Cd   = bubble_drag_cd(Re_b, eps_p_local);

                // Drag force (SI, N) — in direction opposing relative motion
                amrex::Real fd_factor = -0.5 * Cd * rho_f * Ab * ur_mag;
                const amrex::Real FDx = fd_factor * urx;
                const amrex::Real FDy = fd_factor * ury;
                const amrex::Real FDz = fd_factor * urz;

                // Net buoyancy: F = -rho_f * Vb * g in z (upward = +z)
                // (rho_b ≈ 0, so F_g + F_pressure ≈ -rho_f * Vb * g_z)
                // g_grav points downward; buoyancy is upward (+z here)
                const amrex::Real FBx = 0.0;
                const amrex::Real FBy = 0.0;
                const amrex::Real FBz = rho_f * Vb * g;  // upward [N]

                // Effective mass (added mass dominates for gas bubbles)
                const amrex::Real m_eff = C_AM * rho_f * Vb;

                // Acceleration from buoyancy + drag
                const amrex::Real ax_new = (FBx + FDx) / m_eff;
                const amrex::Real ay_new = (FBy + FDy) / m_eff;
                const amrex::Real az_new = (FBz + FDz) / m_eff;

                // Added mass force (Newton's 3rd law back to fluid)
                const amrex::Real FAx = C_AM * rho_f * Vb * ax_new;
                const amrex::Real FAy = C_AM * rho_f * Vb * ay_new;
                const amrex::Real FAz = C_AM * rho_f * Vb * az_new;

                BubbleForces bf;
                bf.fx_a = FAx; bf.fy_a = FAy; bf.fz_a = FAz;
                bf.fx_D = FDx; bf.fy_D = FDy; bf.fz_D = FDz;
                bf.dn_i = 0.0;  // filled later in deposit_o2_sources
                m_forces.push_back(bf);

                // Store new acceleration in particle for next step's velocity Verlet
                p.rdata(BubbleIdx::AX) = ax_new;
                p.rdata(BubbleIdx::AY) = ay_new;
                p.rdata(BubbleIdx::AZ) = az_new;
            }
        }
    }
}

// ============================================================================
// do_breakup
//   Hinze criterion: break if D > D_e = (sigma/rho_f)^0.6 * epsilon^(-0.4)
//   Epsilon is interpolated from the derived MultiFab (EPSILON_IDX),
//   which is computed in compute_derived() using the Smagorinsky SGS model.
//   LB units → SI: eps_SI = eps_LB * dx²/dt³.
// ============================================================================
void BubbleManager::do_breakup(
    const amrex::MultiFab& derived,
    const amrex::Geometry& geom)
{
    if (!m_params.enable_breakup) { return; }

    static std::mt19937 rng_brk(12345);
    std::uniform_real_distribution<amrex::Real> uniform_dist(0.0, 1.0);

    const amrex::Real sigma  = m_params.surf_tension;
    const amrex::Real rho_f  = m_params.rho_fluid;
    const amrex::Real dx     = m_params.dx_phys;
    const amrex::Real dt     = m_params.dt_phys;
    // Conversion factor: eps_LB → eps_SI [m²/s³]
    const amrex::Real eps_conv = dx * dx / (dt * dt * dt);

    amrex::Vector<BubbleParticle> new_bubbles;

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& pbox = kv.second;
            auto& aos  = pbox.GetArrayOfStructs();
            for (auto& p : aos()) {
                if (!p.id().is_valid()) { continue; }

                // Interpolate turbulent dissipation rate from derived field [LB → SI]
                const amrex::Real eps_lb = trilinear_interp(
                    derived, constants::EPSILON_IDX,
                    geom, p.pos(0), p.pos(1), p.pos(2));
                const amrex::Real eps_SI = std::min(
                    std::max(eps_lb * eps_conv, 1.0e-10),
                    m_params.eps_max);  // cap: prevent boundary-layer over-fragmentation

                const amrex::Real d  = p.rdata(BubbleIdx::DIAMETER);

                // Skip if bubble is already at or below minimum allowed size
                if (d <= m_params.min_diameter) { continue; }

                const amrex::Real D_e = std::pow(sigma / rho_f, 0.6) *
                                        std::pow(eps_SI, -0.4);

                if (d <= D_e) { continue; }  // no breakup

                // Compute daughter sizes BEFORE invalidating the parent so we
                // can skip gracefully if either daughter falls below min_diameter.
                amrex::Real fv = 0.5;
                if      (m_params.breakup_model == 0) {
                    // Triangle distribution: lower=0, upper=0.5, mode=0.2
                    // Inverse CDF for triangle distribution
                    const amrex::Real u   = uniform_dist(rng_brk);
                    const amrex::Real low = 0.0, high = 0.5, mode = 0.2;
                    const amrex::Real Fc  = (mode - low) / (high - low);
                    if (u < Fc) {
                        fv = low + std::sqrt(u * (high - low) * (mode - low));
                    } else {
                        fv = high - std::sqrt((1.0 - u) * (high - low) * (high - mode));
                    }
                } else if (m_params.breakup_model == 2) {
                    fv = 0.5 * uniform_dist(rng_brk);
                } else {
                    fv = 0.5; // equal split
                }

                const amrex::Real V_total = (amrex::Math::pi<amrex::Real>() / 6.0) *
                                            std::pow(d, 3.0);
                const amrex::Real V1 = fv * V_total;
                const amrex::Real V2 = (1.0 - fv) * V_total;
                const amrex::Real d1 = std::cbrt(6.0 * V1 / amrex::Math::pi<amrex::Real>());
                const amrex::Real d2 = std::cbrt(6.0 * V2 / amrex::Math::pi<amrex::Real>());
                const amrex::Real n1 = p.rdata(BubbleIdx::N_O2) * fv;
                const amrex::Real n2 = p.rdata(BubbleIdx::N_O2) * (1.0 - fv);

                // Skip breakup if either daughter would fall below min_diameter.
                // Prevents runaway fragmentation from near-wall epsilon spikes.
                if (d1 < m_params.min_diameter || d2 < m_params.min_diameter) {
                    continue;
                }

                // All checks passed — invalidate the parent and create two daughters
                p.id() = -1;

                // Create daughter bubbles (add to new_bubbles list)
                for (int id = 0; id < 2; ++id) {
                    BubbleParticle daughter;
                    daughter.id()  = BubbleParticle::NextID();
                    daughter.cpu() = amrex::ParallelDescriptor::MyProc();
                    daughter.pos(0) = p.pos(0);
                    daughter.pos(1) = p.pos(1);
                    daughter.pos(2) = p.pos(2);
                    daughter.rdata(BubbleIdx::VX) = p.rdata(BubbleIdx::VX);
                    daughter.rdata(BubbleIdx::VY) = p.rdata(BubbleIdx::VY);
                    daughter.rdata(BubbleIdx::VZ) = p.rdata(BubbleIdx::VZ);
                    daughter.rdata(BubbleIdx::AX) = p.rdata(BubbleIdx::AX);
                    daughter.rdata(BubbleIdx::AY) = p.rdata(BubbleIdx::AY);
                    daughter.rdata(BubbleIdx::AZ) = p.rdata(BubbleIdx::AZ);
                    daughter.rdata(BubbleIdx::DIAMETER) = (id == 0) ? d1 : d2;
                    daughter.rdata(BubbleIdx::N_O2)     = (id == 0) ? n1 : n2;
                    new_bubbles.push_back(daughter);
                }
            }
        }
    }

    // Add new daughters to the container.
    // Each rank stages its daughters in its own tile (0,0) then Redistribute
    // scatters them to the correct process/tile based on position.
    if (!new_bubbles.empty()) {
        auto& ptile = m_container.DefineAndReturnParticleTile(0, 0, 0);
        for (auto& nb : new_bubbles) {
            ptile.push_back(nb);
        }
    }

    // Redistribute: removes particles with id=-1 and scatters new daughters
    m_container.Redistribute();
}

// ============================================================================
// do_coalescence
//   Boshenyatov criterion: Re_a > 40 → coalesce; else elastic bounce.
//   O(N²) particle-pair check — adequate for N < 10^4.
// ============================================================================
void BubbleManager::do_coalescence(amrex::Real phys_time)
{
    if (!m_params.enable_coalescence) { return; }
    if (phys_time < m_params.coal_start_time) { return; }

    const amrex::Real nu_f  = m_params.nu_fluid;
    const amrex::Real dx    = m_params.dx_phys;

    // Collect all valid bubble data CPU-side
    struct BData {
        int lev; amrex::Long idx_in_box;
        amrex::Real x, y, z, vx, vy, vz, d, n_O2;
        bool valid = true;
    };
    amrex::Vector<BData> bvec;

    // Map from global index to (lev, grid, tile, local_index) for merging
    struct BRef { int lev; int grid; int tile; int local_idx; };
    amrex::Vector<BRef> brefs;

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& pbox  = kv.second;
            auto& aos   = pbox.GetArrayOfStructs();
            int grid    = kv.first.first;
            int tile    = kv.first.second;
            for (int li = 0; li < static_cast<int>(aos().size()); ++li) {
                const auto& p = aos()[li];
                if (!p.id().is_valid()) { continue; }
                BData bd;
                bd.x   = p.pos(0) * dx;  // convert to metres
                bd.y   = p.pos(1) * dx;
                bd.z   = p.pos(2) * dx;
                bd.vx  = p.rdata(BubbleIdx::VX);
                bd.vy  = p.rdata(BubbleIdx::VY);
                bd.vz  = p.rdata(BubbleIdx::VZ);
                bd.d   = p.rdata(BubbleIdx::DIAMETER);
                bd.n_O2 = p.rdata(BubbleIdx::N_O2);
                bd.valid = true;
                bvec.push_back(bd);
                BRef br{lev, grid, tile, li};
                brefs.push_back(br);
            }
        }
    }

    const int N = static_cast<int>(bvec.size());
    for (int i = 0; i < N; ++i) {
        if (!bvec[i].valid) { continue; }
        for (int j = i+1; j < N; ++j) {
            if (!bvec[j].valid) { continue; }

            const amrex::Real dxi = bvec[i].d;
            const amrex::Real dxj = bvec[j].d;
            const amrex::Real d_h = 2.0 * dxi * dxj / (dxi + dxj);  // harmonic diameter

            // Centre-to-centre distance
            const amrex::Real rx = bvec[i].x - bvec[j].x;
            const amrex::Real ry = bvec[i].y - bvec[j].y;
            const amrex::Real rz = bvec[i].z - bvec[j].z;
            const amrex::Real dist = std::sqrt(rx*rx + ry*ry + rz*rz);

            // Only consider overlapping or very close bubbles
            if (dist > 0.5 * (dxi + dxj) * 1.5) { continue; }

            // Approach velocity (component along line of centres)
            const amrex::Real dvx = bvec[i].vx - bvec[j].vx;
            const amrex::Real dvy = bvec[i].vy - bvec[j].vy;
            const amrex::Real dvz = bvec[i].vz - bvec[j].vz;
            const amrex::Real dv_mag = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

            const amrex::Real Re_a = dv_mag * d_h / nu_f;

            if (Re_a > m_params.Re_a_crit) {
                // Coalescence: conserve volume and momentum
                const amrex::Real V_i = (amrex::Math::pi<amrex::Real>() / 6.0) * dxi*dxi*dxi;
                const amrex::Real V_j = (amrex::Math::pi<amrex::Real>() / 6.0) * dxj*dxj*dxj;
                const amrex::Real V_new = V_i + V_j;
                const amrex::Real d_new = std::cbrt(6.0 * V_new / amrex::Math::pi<amrex::Real>());
                // Momentum conservation (using m ≈ 0, so mass average is unweighted)
                const amrex::Real vx_new = 0.5 * (bvec[i].vx + bvec[j].vx);
                const amrex::Real vy_new = 0.5 * (bvec[i].vy + bvec[j].vy);
                const amrex::Real vz_new = 0.5 * (bvec[i].vz + bvec[j].vz);
                const amrex::Real n_new  = bvec[i].n_O2 + bvec[j].n_O2;
                // Position: centre-of-volume
                const amrex::Real x_new = 0.5 * (bvec[i].x + bvec[j].x);
                const amrex::Real y_new = 0.5 * (bvec[i].y + bvec[j].y);
                const amrex::Real z_new = 0.5 * (bvec[i].z + bvec[j].z);

                // Update particle i with coalesced data
                auto& br_i = brefs[i];
                auto& aos_i = m_container.GetParticles(br_i.lev)[{br_i.grid, br_i.tile}]
                              .GetArrayOfStructs();
                auto& pi = aos_i()[br_i.local_idx];
                pi.pos(0) = x_new / dx;
                pi.pos(1) = y_new / dx;
                pi.pos(2) = z_new / dx;
                pi.rdata(BubbleIdx::VX)       = vx_new;
                pi.rdata(BubbleIdx::VY)       = vy_new;
                pi.rdata(BubbleIdx::VZ)       = vz_new;
                pi.rdata(BubbleIdx::DIAMETER) = d_new;
                pi.rdata(BubbleIdx::N_O2)     = n_new;

                // Invalidate particle j
                auto& br_j = brefs[j];
                auto& aos_j = m_container.GetParticles(br_j.lev)[{br_j.grid, br_j.tile}]
                              .GetArrayOfStructs();
                aos_j()[br_j.local_idx].id() = -1;

                bvec[i].d   = d_new;
                bvec[i].vx  = vx_new; bvec[i].vy = vy_new; bvec[i].vz = vz_new;
                bvec[i].n_O2 = n_new;
                bvec[i].x   = x_new; bvec[i].y = y_new; bvec[i].z = z_new;
                bvec[j].valid = false;
            }
            // else: elastic bounce — for now do nothing (simple model)
        }
    }
    m_container.Redistribute();
}

// ============================================================================
// remove_exited_bubbles
//   Remove bubbles that have risen above the free surface.
// ============================================================================
void BubbleManager::remove_exited_bubbles(const amrex::Geometry& geom,
                                           const amrex::MultiFab* phi_mf)
{
    const amrex::Real z_surf = m_params.free_surface_z;

    if (phi_mf != nullptr) {
        // Phase-field mode: remove bubble when Φ < 0.5 at its cell.
        // Uses nearest-cell lookup (floor of position in LB coordinates).
        const auto& prob_lo = geom.ProbLoArray();
        const auto& dx_arr  = geom.CellSizeArray();
        const amrex::Box domain = geom.Domain();

        for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
            for (auto& kv : m_container.GetParticles(lev)) {
                auto& aos = kv.second.GetArrayOfStructs();
                for (auto& p : aos()) {
                    if (!p.id().is_valid()) { continue; }

                    // Cell index containing bubble centre
                    const int ci = static_cast<int>(std::floor(
                        (p.pos(0) - prob_lo[0]) / dx_arr[0]));
                    const int cj = static_cast<int>(std::floor(
                        (p.pos(1) - prob_lo[1]) / dx_arr[1]));
                    const int ck = static_cast<int>(std::floor(
                        (p.pos(2) - prob_lo[2]) / dx_arr[2]));

                    const amrex::IntVect iv(AMREX_D_DECL(ci, cj, ck));
                    if (!domain.contains(iv)) {
                        // Outside domain entirely — always remove
                        p.id() = -1;
                        continue;
                    }

                    // Find the MFIter tile that owns this cell
                    for (amrex::MFIter mfi(*phi_mf); mfi.isValid(); ++mfi) {
                        if (mfi.validbox().contains(iv)) {
                            const amrex::Real phi_val =
                                (*phi_mf).array(mfi)(iv, 0);
                            if (phi_val < 0.5) {
                                p.id() = -1;  // in gas region — exit
                            }
                            break;
                        }
                    }
                }
            }
        }
    } else {
        // Fallback: fixed z-threshold
        for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
            for (auto& kv : m_container.GetParticles(lev)) {
                auto& aos = kv.second.GetArrayOfStructs();
                for (auto& p : aos()) {
                    if (!p.id().is_valid()) { continue; }
                    if (p.pos(2) >= z_surf) {
                        p.id() = -1;
                    }
                }
            }
        }
    }
    m_container.Redistribute();
}

// ============================================================================
// deposit_fluid_forces
//   Convert force (SI, N) to LB force density and add into fluid_force_mf.
//   Uses a two-pass approach: collect (cell, value) pairs, then iterate
//   over MultiFab tiles and add contributions (safe, no GPU atomics needed).
// ============================================================================
void BubbleManager::deposit_fluid_forces(
    amrex::MultiFab&       fluid_force_mf,
    const amrex::Geometry& geom)
{
    const amrex::Real dx        = m_params.dx_phys;
    const amrex::Real dt        = m_params.dt_phys;
    const amrex::Real rho_ref   = m_params.rho_fluid;
    // F_LB_vol = F_SI [N] * dt^2 / (rho_ref [kg/m³] * dx^4 [m^4])
    const amrex::Real conv = dt * dt / (rho_ref * std::pow(dx, 4.0));

    // Pass 1: collect all force deposits
    struct ForceDep { int ci, cj, ck; amrex::Real fx, fy, fz; };
    amrex::Vector<ForceDep> deps;
    deps.reserve(static_cast<int>(m_forces.size()));

    const amrex::Real* plo    = geom.ProbLo();
    const amrex::Real* dx_arr = geom.CellSize();
    const amrex::Box& dom     = geom.Domain();

    {
        int fi = 0;
        for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
            for (const auto& kv : m_container.GetParticles(lev)) {
                const auto& aos = kv.second.GetArrayOfStructs();
                for (const auto& p : aos()) {
                    if (!p.id().is_valid()) { ++fi; continue; }

                    const auto& bf = m_forces[fi];
                    const amrex::Real Ffx = -(bf.fx_a + bf.fx_D) * conv;
                    const amrex::Real Ffy = -(bf.fy_a + bf.fy_D) * conv;
                    const amrex::Real Ffz = -(bf.fz_a + bf.fz_D) * conv;

                    const int ci = static_cast<int>(std::floor((p.pos(0) - plo[0]) / dx_arr[0]));
                    const int cj = static_cast<int>(std::floor((p.pos(1) - plo[1]) / dx_arr[1]));
                    const int ck = static_cast<int>(std::floor((p.pos(2) - plo[2]) / dx_arr[2]));

                    if (dom.contains(amrex::IntVect(AMREX_D_DECL(ci, cj, ck)))) {
                        deps.push_back({ci, cj, ck, Ffx, Ffy, Ffz});
                    }
                    ++fi;
                }
            }
        }
    }

    // Pass 2: deposit into MultiFab (MFIter outer loop)
    for (amrex::MFIter mfi(fluid_force_mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = fluid_force_mf.array(mfi);
        for (const auto& d : deps) {
            if (bx.contains(amrex::IntVect(AMREX_D_DECL(d.ci, d.cj, d.ck)))) {
                arr(d.ci, d.cj, d.ck, 0) += d.fx;
                arr(d.ci, d.cj, d.ck, 1) += d.fy;
                arr(d.ci, d.cj, d.ck, 2) += d.fz;
            }
        }
    }
}

// ============================================================================
// deposit_o2_sources
//   Kawase k_L, Eq. 17 of Thomas et al.
//   k_L = 0.301 * (epsilon * nu)^(1/4) * Sc^(-1/2)
//
//   Mass transfer rate per bubble [mol/s]:
//   dn_i/dt = k_L * A_i * (S_i * C_g,i - C_f,i)
//
//   Source deposited on cell j [mol/(m³·s)]:
//   cdot_j = (1/Vcell) * sum_{i in j} dn_i
// ============================================================================
void BubbleManager::deposit_o2_sources(
    amrex::MultiFab&       o2_src_mf,
    const amrex::MultiFab& o2_conc_mf,
    const amrex::MultiFab& derived,
    const amrex::Geometry& geom)
{
    const amrex::Real dx     = m_params.dx_phys;
    const amrex::Real dt     = m_params.dt_phys;
    const amrex::Real nu_f   = m_params.nu_fluid;
    const amrex::Real D_O2   = m_params.D_O2;
    const amrex::Real C_kL   = m_params.kL_coeff;
    const amrex::Real Sc     = nu_f / D_O2;            // Schmidt number
    const amrex::Real S_i    = m_params.O2_solubility; // dimensionless Bunsen. 
    const amrex::Real dx3    = dx * dx * dx;            // cell volume [m³]
    // epsilon LB → SI: eps_SI [m²/s³] = eps_LB * dx²/dt³
    const amrex::Real eps_conv = dx * dx / (dt * dt * dt);

    // Pass 1: compute per-bubble mass transfer and collect deposits
    struct SrcDep { int ci, cj, ck; amrex::Real src; };
    amrex::Vector<SrcDep> deps;
    deps.reserve(static_cast<int>(m_forces.size()));

    const amrex::Real* plo    = geom.ProbLo();
    const amrex::Real* dx_arr = geom.CellSize();
    const amrex::Box& dom     = geom.Domain();

    {
        int fi = 0;
        for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
            for (auto& kv : m_container.GetParticles(lev)) {
                auto& pbox = kv.second;
                auto& aos  = pbox.GetArrayOfStructs();
                for (auto& p : aos()) {
                    if (!p.id().is_valid()) { ++fi; continue; }

                    const amrex::Real d   = p.rdata(BubbleIdx::DIAMETER);
                    const amrex::Real r   = 0.5 * d;
                    const amrex::Real Vb  = (4.0/3.0) * amrex::Math::pi<amrex::Real>() * r*r*r;
                    const amrex::Real Ab  = amrex::Math::pi<amrex::Real>() * d * d;

                    // Gas-phase O2 concentration [mol/m³]
                    const amrex::Real n_O2  = p.rdata(BubbleIdx::N_O2);
                    const amrex::Real C_g_i = (Vb > 0.0) ? (n_O2 / Vb) : 0.0;

                    // Fluid O2 concentration at bubble location [LB rho → mol/m³]
                    const amrex::Real C_f_lb_raw = trilinear_interp(o2_conc_mf, 0,
                                                                  geom,
                                                                  p.pos(0), p.pos(1), p.pos(2));

                    // Clamp negative concentrations to zero: negative values are
                    // numerical artifacts from early-time oscillations and have no
                    // physical meaning.  Treating them as zero keeps the driving
                    // force positive (bubble→liquid) which is stabilizing.
                    const amrex::Real C_f_lb = std::max(C_f_lb_raw, amrex::Real(0.0));

                    // Convert from LB_rho to physical units: 1 LB_rho = C_ref mol/m³
                    const amrex::Real C_f_i = C_f_lb * m_params.C_ref;

                    // Local epsilon [LB → SI]
                    const amrex::Real eps_lb = trilinear_interp(derived, constants::EPSILON_IDX,
                                                                  geom,
                                                                  p.pos(0), p.pos(1), p.pos(2));
                    const amrex::Real eps_SI = std::max(eps_lb * eps_conv, 1.0e-10);

                    // Kawase k_L [m/s]  (Eq. 17, Thomas et al.)
                    const amrex::Real k_L = C_kL * std::pow(eps_SI * nu_f, 0.25) *
                                            std::pow(Sc, -0.5);

                    // Driving force [mol/m³]: S_i * C_g,i - C_f,i
                    const amrex::Real driving = S_i * C_g_i - C_f_i;

                    // Mass transfer rate [mol/s]
                    const amrex::Real dn_i = k_L * Ab * driving;

                    // Cache for budget output
                    if (fi < static_cast<int>(m_forces.size())) {
                        m_forces[fi].dn_i = dn_i;
                    }

                    // Shrink bubble O2 content
                    amrex::Real n_O2_new = n_O2 - dn_i * dt;
                    if (n_O2_new < 0.0) { n_O2_new = 0.0; }
                    p.rdata(BubbleIdx::N_O2) = n_O2_new;

                    // Update bubble diameter from new mole count
                    const amrex::Real Vb_new = n_O2_new * m_params.O2_molar_volume;
                    // Minimum bubble diameter (10 µm): below this the effective
                    // mass → 0 causing acceleration → Inf in compute_forces.
                    const amrex::Real d_min = 1.0e-5;  // metres
                    if (Vb_new > 0.0) {
                        const amrex::Real d_new =
                            std::cbrt(6.0 * Vb_new / amrex::Math::pi<amrex::Real>());
                        if (d_new < d_min) {
                            p.id() = -1;  // Too small — remove
                        } else {
                            p.rdata(BubbleIdx::DIAMETER) = d_new;
                        }
                    } else {
                        p.id() = -1;  // Bubble fully dissolved — mark for removal
                    }

                    // Source deposit [mol/(m³·s)]
                    const amrex::Real src_rate = dn_i / dx3;

                    const int ci = static_cast<int>(std::floor((p.pos(0) - plo[0]) / dx_arr[0]));
                    const int cj = static_cast<int>(std::floor((p.pos(1) - plo[1]) / dx_arr[1]));
                    const int ck = static_cast<int>(std::floor((p.pos(2) - plo[2]) / dx_arr[2]));
                    if (dom.contains(amrex::IntVect(AMREX_D_DECL(ci, cj, ck)))) {
                        deps.push_back({ci, cj, ck, src_rate});
                    }
                    ++fi;
                }
            }
        }
    }

    // Pass 2: deposit into o2_src_mf (MFIter outer loop)
    amrex::Real dep_max = 0.0;
    for (amrex::MFIter mfi(o2_src_mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto arr = o2_src_mf.array(mfi);
        for (const auto& d : deps) {
            if (bx.contains(amrex::IntVect(AMREX_D_DECL(d.ci, d.cj, d.ck)))) {
                arr(d.ci, d.cj, d.ck, 0) += d.src;
                dep_max = std::max(dep_max, std::abs(d.src));
            }
        }
    }
    // Diagnostic: print deposit summary (controlled by a local print interval)
    {
        static int dbg_ctr = 0;
        if (++dbg_ctr % 800 == 1) {
            amrex::Print() << "[O2_deposit] n_deps=" << deps.size()
                           << "  max_src_rate=" << dep_max << " mol/(m3*s)\n";
        }
    }

    // Remove bubbles that fully dissolved (id = -1)
    m_container.Redistribute();
}

// ============================================================================
// advance  — main per-step driver
// ============================================================================
void BubbleManager::advance(
    amrex::Real                dt,
    const amrex::MultiFab&     macrodata,
    const amrex::MultiFab&     derived,
    const amrex::MultiFab&     o2_conc_mf,
    const amrex::Geometry&     geom,
    amrex::MultiFab&           fluid_force_mf,
    amrex::MultiFab&           o2_src_mf,
    amrex::Real                phys_time,
    const amrex::MultiFab*     phi_mf)
{
    if (!m_initialized) { return; }
    if (!m_particles_ever_injected) { return; }  // m_particles not yet sized
    const amrex::Real dt_phys = m_params.dt_phys;
    const amrex::Real dx_phys = m_params.dx_phys;

    // ------------------------------------------------------------------
    // 0. Copy device MultiFabs to host-pinned copies so that all CPU-side
    //    particle loops (trilinear_interp, deposit_*) can read/write them
    //    without requiring amrex.the_arena_is_managed=1.
    //    Pinned memory is directly accessible from CPU and DMA-accessible
    //    by the GPU, so the copies back to device are cheap H2D transfers.
    // ------------------------------------------------------------------
    auto make_pinned_copy = [](const amrex::MultiFab& src) {
        amrex::MultiFab h(src.boxArray(), src.DistributionMap(),
                          src.nComp(), 0,
                          amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        amrex::MultiFab::Copy(h, src, 0, 0, src.nComp(), 0);
        return h;
    };

    // Input MFs: device → host-pinned
    amrex::MultiFab macrodata_h = make_pinned_copy(macrodata);
    amrex::MultiFab derived_h   = make_pinned_copy(derived);
    amrex::MultiFab o2_conc_h   = make_pinned_copy(o2_conc_mf);
    amrex::Gpu::streamSynchronize();  // ensure D2H copies are done before CPU reads

    // Output MFs: zeroed host-pinned buffers (will be copied back to device)
    amrex::MultiFab fluid_force_h(fluid_force_mf.boxArray(),
                                  fluid_force_mf.DistributionMap(),
                                  fluid_force_mf.nComp(), 0,
                                  amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    fluid_force_h.setVal(0.0);
    amrex::MultiFab o2_src_h(o2_src_mf.boxArray(),
                             o2_src_mf.DistributionMap(),
                             o2_src_mf.nComp(), 0,
                             amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    o2_src_h.setVal(0.0);

    // Optional phi_mf: copy if provided, otherwise leave null
    std::unique_ptr<amrex::MultiFab> phi_h;
    const amrex::MultiFab* phi_h_ptr = nullptr;
    if (phi_mf != nullptr) {
        phi_h = std::make_unique<amrex::MultiFab>(
            phi_mf->boxArray(), phi_mf->DistributionMap(),
            phi_mf->nComp(), 0,
            amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        amrex::MultiFab::Copy(*phi_h, *phi_mf, 0, 0, phi_mf->nComp(), 0);
        amrex::Gpu::streamSynchronize();
        phi_h_ptr = phi_h.get();
    }

    // ------------------------------------------------------------------
    // 1. Boyle's law diameter correction (reads only particle data — pinned)
    // ------------------------------------------------------------------
    apply_boyles_law(geom);

    // ------------------------------------------------------------------
    // 2. Compute forces at current positions (fills m_forces, stores 
    //    acceleration in particle rdata for velocity Verlet)
    // ------------------------------------------------------------------
    compute_forces(macrodata_h, derived_h, geom);

    // ------------------------------------------------------------------
    // 3. Deposit two-way coupling body forces on fluid.
    //    MUST happen BEFORE any Redistribute() call so that m_forces[]
    //    remains indexed in the same particle-iteration order as step 2.
    //    Redistribute (called below in Verlet, remove_exited_bubbles, and
    //    do_breakup) may reorder particles within tiles, breaking the fi
    //    index alignment if deposits were done after those calls.
    // ------------------------------------------------------------------
    deposit_fluid_forces(fluid_force_h, geom);

    // ------------------------------------------------------------------
    // 4. Mass transfer O2 source + bubble shrinkage.
    //    Also done before Redistribute for the same index-safety reason.
    // ------------------------------------------------------------------
    deposit_o2_sources(o2_src_h, o2_conc_h, derived_h, geom);

    // ------------------------------------------------------------------
    // 5. Velocity Verlet integration
    //    x(t+dt) = x(t) + v(t)*dt_phys + 0.5*a(t)*dt_phys²  [LB cells]
    //    v(t+dt) = v(t) + a(t)*dt_phys                       [m/s, SI]
    //    (Full 2nd-order Verlet; a from compute_forces above)
    //    Redistribute after this call may reorder particle storage.
    // ------------------------------------------------------------------
    {
        for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
            for (auto& kv : m_container.GetParticles(lev)) {
                auto& pbox = kv.second;
                auto& aos  = pbox.GetArrayOfStructs();
                for (auto& p : aos()) {
                    if (!p.id().is_valid()) { continue; }

                    const amrex::Real ax = p.rdata(BubbleIdx::AX);
                    const amrex::Real ay = p.rdata(BubbleIdx::AY);
                    const amrex::Real az = p.rdata(BubbleIdx::AZ);

                    const amrex::Real vx = p.rdata(BubbleIdx::VX);
                    const amrex::Real vy = p.rdata(BubbleIdx::VY);
                    const amrex::Real vz = p.rdata(BubbleIdx::VZ);

                    // Position update (to LB cells): Δx_LB = v_SI * dt_phys / dx_phys
                    p.pos(0) += (vx * dt_phys + 0.5 * ax * dt_phys * dt_phys) / dx_phys;
                    p.pos(1) += (vy * dt_phys + 0.5 * ay * dt_phys * dt_phys) / dx_phys;
                    p.pos(2) += (vz * dt_phys + 0.5 * az * dt_phys * dt_phys) / dx_phys;

                    // Velocity update (SI, m/s)
                    p.rdata(BubbleIdx::VX) = vx + ax * dt_phys;
                    p.rdata(BubbleIdx::VY) = vy + ay * dt_phys;
                    p.rdata(BubbleIdx::VZ) = vz + az * dt_phys;
                }
            }
        }
    }
    m_container.Redistribute();

    // ------------------------------------------------------------------
    // 6. Remove bubbles that have crossed the free surface.
    //    When phi_mf is provided (Chiu & Lin phase-field active), use
    //    Φ < 0.5 at the bubble's cell as the exit criterion — consistent
    //    with the dynamic interface position.  Otherwise fall back to
    //    the fixed free_surface_z z-coordinate threshold.
    // ------------------------------------------------------------------
    remove_exited_bubbles(geom, phi_h_ptr);

    // ------------------------------------------------------------------
    // 7. Breakup check at new positions
    // ------------------------------------------------------------------
    do_breakup(derived_h, geom);

    // ------------------------------------------------------------------
    // 8. Copy host-pinned output buffers back to device MultiFabs.
    //    fluid_force_h and o2_src_h were filled by deposit_fluid_forces
    //    and deposit_o2_sources above; add them into the device MFs so
    //    the LBM step can apply them as body forces.
    // ------------------------------------------------------------------
    amrex::MultiFab::Add(fluid_force_mf, fluid_force_h, 0, 0,
                         fluid_force_mf.nComp(), 0);
    amrex::MultiFab::Add(o2_src_mf, o2_src_h, 0, 0,
                         o2_src_mf.nComp(), 0);
}

// ============================================================================
// write_stats
// ============================================================================
void BubbleManager::write_stats(int step, amrex::Real phys_time)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }

    // Gather statistics across all particles
    int    n_bub   = 0;
    amrex::Real d_sum  = 0.0;
    amrex::Real d_min  = 1.0e30;
    amrex::Real d_max  = 0.0;
    amrex::Real n_O2_total = 0.0;
    amrex::Real dn_total   = 0.0;

    int fi = 0;
    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& aos = kv.second.GetArrayOfStructs();
            for (auto& p : aos()) {
                if (!p.id().is_valid()) { ++fi; continue; }
                ++n_bub;
                const amrex::Real d = p.rdata(BubbleIdx::DIAMETER);
                d_sum     += d;
                d_min      = std::min(d_min, d);
                d_max      = std::max(d_max, d);
                n_O2_total += p.rdata(BubbleIdx::N_O2);
                if (fi < static_cast<int>(m_forces.size())) {
                    dn_total += m_forces[fi].dn_i;
                }
                ++fi;
            }
        }
    }

    const amrex::Real d_mean = (n_bub > 0) ? (d_sum / n_bub) : 0.0;

    if (!m_stats_stream.is_open()) { return; }
    m_stats_stream << std::fixed << std::setprecision(6)
                   << step << ","
                   << phys_time << ","
                   << n_bub << ","
                   << d_mean * 1000.0 << ","   // mm
                   << d_min  * 1000.0 << ","
                   << d_max  * 1000.0 << ","
                   << n_O2_total << ","
                   << dn_total << "\n";
    m_stats_stream.flush();
}

// ============================================================================
// Stats file helpers
// ============================================================================
void BubbleManager::open_stats_file()
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    m_stats_stream.open(m_params.stats_file,
                        std::ios::out | std::ios::trunc);
    m_stats_stream << "step,phys_time_s,n_bubbles,d_mean_mm,"
                      "d_min_mm,d_max_mm,n_O2_total_mol,dn_O2_step_mol_per_s\n";
}

void BubbleManager::close_stats_file()
{
    if (m_stats_stream.is_open()) { m_stats_stream.close(); }
}

} // namespace lbm
