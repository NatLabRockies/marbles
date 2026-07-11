/**
 * @file Bubble.cpp
 * @brief Implementation of Lagrangian bubble tracking for kLa mass transfer.
 *
 * Reference:
 *   Thomas et al. (2021), CES 237, 116538. DOI: 10.1016/j.ces.2021.116538
 *   Rettinger & Rüde (2018), Comput Fluids 172 — point-particle LBM–DEM
 * coupling (cited at Eq. 9 of Thomas et al.) Tenneti, Garg & Subramaniam
 * (2011), Int J Multiph Flow 37(9), 1072–1092 — drag law (Eq. 18) Kawase et al.
 * (1992), Biotechnol. Bioeng. 39(11) — k_L formula Boshenyatov (2012) —
 * coalescence criterion Re_a > 40 Hinze (1955) — breakup equilibrium diameter
 */
#include "Bubble.H"
#include "Constants.H"
#include <cstdint>
#include <cstring>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <sstream>
#include <iomanip>

// -----------------------------------------------------------------------------
// SYCL is stricter than CUDA/HIP about which types can be captured by value
// into a device kernel: it statically checks sycl::is_device_copyable.
// BubbleParams contains amrex::Vector<Real> (sparger_x/y/z) and std::string
// (stats_file) which are not trivially copyable, so the check fails even
// though the device code path only touches the scalar POD members of the
// struct.  Opt in explicitly so the SYCL build accepts capture; CUDA / HIP
// are unaffected because they never see this specialization.
// -----------------------------------------------------------------------------
#ifdef AMREX_USE_SYCL
template <>
struct sycl::is_device_copyable<lbm::BubbleParams> : std::true_type
{};
#endif

namespace lbm {

// -----------------------------------------------------------------------------
// POD subset of BubbleParams captured by the per-bubble ParallelFor lambda in
// BubbleManager::advance().  Capturing the full BubbleParams struct (which
// contains amrex::Vector<Real> sparger_x/y/z and std::string stats_file) by
// value would make the lambda's own type non-trivially-copyable, which SYCL
// rejects at kernel-submit time (sycl::is_device_copyable) even though the
// kernel body only reads these scalar members.  Structs-of-scalars are POD, so
// the check passes cleanly on all backends (CUDA, HIP, SYCL, CPU).
//
// This must live at namespace scope: CUDA rejects function-local types in the
// captured types of an extended __device__ / __host__ __device__ lambda.
// -----------------------------------------------------------------------------
struct BubbleParamsScalar
{
    amrex::Real dx_phys;
    amrex::Real dt_phys;
    amrex::Real nu_fluid;
    amrex::Real rho_fluid;
    amrex::Real g_grav;
    amrex::Real force_cap_factor;
    amrex::Real C_ref;
    amrex::Real D_O2;
    amrex::Real kL_coeff;
    amrex::Real O2_solubility;
    amrex::Real free_surface_z;
    amrex::Real P_atm;
    amrex::Real O2_molar_volume;
    int require_liquid_host;
    int boyle_law_enable;
};

// ============================================================================
// Boyle's law / hydrostatic-pressure helper (June 2026).
//
// Returns the depth-corrected molar volume V_m for a bubble at vertical
// coordinate z_LB (in lattice-Boltzmann cell units).  When the Boyle flag
// is off, this collapses to V_m,ref so all conversions n_O2 ↔ V_b ↔ d
// reduce to the legacy ideal-gas-at-STP model.
//
// When the flag is on:
//   h_phys = max(0, free_surface_z − z_LB) · dx_phys     [m]
//   P(h)   = P_atm + ρ_l · g · h_phys                    [Pa]
//   V_m    = V_m,ref · P_atm / P(h)                      [m³/mol]
//
// h is clamped at 0 so a bubble above the surface doesn't get a negative
// hydrostatic pressure (would have produced V_m > V_m,ref → ρ-expansion).
//
// v_m_at_depth is host-and-device callable (used by inject_bubbles on the
// host and by the per-bubble update kernel on the GPU).
//
// Templated on the parameter struct so the same function can be called with
// either the full BubbleParams (from host code) or a POD scalar-subset copy
// (from the GPU lambda in advance()).  This avoids capturing the whole
// BubbleParams — which contains non-trivially-copyable amrex::Vector /
// std::string members — into the SYCL device kernel.  The lambda parameter
// P must expose these scalar fields: boyle_law_enable, O2_molar_volume,
// free_surface_z, dx_phys, P_atm, rho_fluid, g_grav.
// ============================================================================
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
v_m_at_depth(amrex::Real z_LB, const P& p)
{
    if (p.boyle_law_enable == 0) {
        return p.O2_molar_volume;
    }
    const amrex::Real h_phys =
        amrex::max(amrex::Real(0.0), (p.free_surface_z - z_LB) * p.dx_phys);
    const amrex::Real pressure = p.P_atm + p.rho_fluid * p.g_grav * h_phys;
    return p.O2_molar_volume * (p.P_atm / pressure);
}

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
//              = 0.5 * (24/Re) * (1/eps_f^3) * rho_f * A_b * |u_rel|^2 * eps_f
//              * eps_f^3 * F/F_SN
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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::Real
bubble_drag_cd(amrex::Real Re, amrex::Real eps_p)
{
    if (Re < 1.0e-10) {
        return 0.0;
    }

    // Guard: clamp eps_p to physically meaningful range
    eps_p = std::max(amrex::Real(0.0), std::min(eps_p, amrex::Real(0.99)));
    const amrex::Real eps_f = 1.0 - eps_p;
    const amrex::Real eps_f3 = eps_f * eps_f * eps_f;

    // F_0: Schiller-Naumann corrected by 1/eps_f^3
    const amrex::Real F0 = (1.0 + 0.15 * std::pow(Re, 0.687)) / eps_f3;

    // F_1: dense-packing term
    const amrex::Real F1 = 10.0 * eps_p / (eps_f * eps_f);

    // F_2: cross term (Re × eps_p), Tenneti et al. Eq. 18
    amrex::Real F2 = 0.0;
    if (eps_p > 1.0e-8) {
        const amrex::Real numer =
            eps_f / eps_f3 + 3.0 * eps_p * eps_f + 8.4 * std::pow(Re, -0.343);
        const amrex::Real denom =
            1.0 + std::pow(10.0, 3.0 * eps_p) *
                      std::pow(Re, -(1.0 + 4.0 * eps_p) / 2.0);
        F2 = (eps_p / eps_f3) * (0.413 * Re / 24.0) * (numer / denom);
    }

    // F = F0 + F1 + F2  (Tenneti Eq. 18)
    const amrex::Real F = F0 + F1 + F2;

    // Convert to C_D convention: F_d = 0.5 * C_D * rho_f * A_b * |u_rel|^2
    // F_d(Tenneti) = 3π μ d |u_rel| eps_f F = 0.5 * C_D * rho_f * π*(d/2)^2 *
    // |u_rel|^2 → C_D = 24*eps_f*F / Re
    const amrex::Real Cd = 24.0 * eps_f * F / Re;

    // High-Re cap (Tenneti correlation is validated for Re < ~300)
    return (Re < 1000.0) ? Cd : std::min(Cd, amrex::Real(0.44));
}

// ============================================================================
// Utility: trilinear interpolation in a MultiFab at LB-cell position (px,py,pz)

// ============================================================================
// Interpolation functions
// ============================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real trilinear_interp_device(
    const amrex::Array4<const amrex::Real>& arr,
    int comp,
    const amrex::Real* plo,
    const amrex::Real* dx_arr,
    amrex::Real px,
    amrex::Real py,
    amrex::Real pz)
{
    amrex::Real fi = (px - plo[0]) / dx_arr[0] - 0.5;
    amrex::Real fj = (py - plo[1]) / dx_arr[1] - 0.5;
    amrex::Real fk = (pz - plo[2]) / dx_arr[2] - 0.5;

    int i0 = static_cast<int>(amrex::Math::floor(fi));
    int j0 = static_cast<int>(amrex::Math::floor(fj));
    int k0 = static_cast<int>(amrex::Math::floor(fk));
    int i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1;

    amrex::Real ti = fi - i0;
    amrex::Real tj = fj - j0;
    amrex::Real tk = fk - k0;

    auto lb = amrex::lbound(arr);
    auto ub = amrex::ubound(arr);

    auto clamp = [](int v, int l, int h) {
        return amrex::max(l, amrex::min(h, v));
    };
    i0 = clamp(i0, lb.x, ub.x);
    i1 = clamp(i1, lb.x, ub.x);
    j0 = clamp(j0, lb.y, ub.y);
    j1 = clamp(j1, lb.y, ub.y);
    k0 = clamp(k0, lb.z, ub.z);
    k1 = clamp(k1, lb.z, ub.z);

    auto get = [&](int ii, int jj, int kk) { return arr(ii, jj, kk, comp); };

    return (1 - ti) * (1 - tj) * (1 - tk) * get(i0, j0, k0) +
           ti * (1 - tj) * (1 - tk) * get(i1, j0, k0) +
           (1 - ti) * tj * (1 - tk) * get(i0, j1, k0) +
           ti * tj * (1 - tk) * get(i1, j1, k0) +
           (1 - ti) * (1 - tj) * tk * get(i0, j0, k1) +
           ti * (1 - tj) * tk * get(i1, j0, k1) +
           (1 - ti) * tj * tk * get(i0, j1, k1) +
           ti * tj * tk * get(i1, j1, k1);
}

amrex::Real trilinear_interp(
    const amrex::MultiFab& mf,
    int comp,
    const amrex::Geometry& geom,
    amrex::Real px,
    amrex::Real py,
    amrex::Real pz)
{
    const amrex::Real* plo = geom.ProbLo();
    const amrex::Real* dx_arr = geom.CellSize();

    amrex::Real fi = (px - plo[0]) / dx_arr[0] - 0.5;
    amrex::Real fj = (py - plo[1]) / dx_arr[1] - 0.5;
    amrex::Real fk = (pz - plo[2]) / dx_arr[2] - 0.5;

    int i0 = static_cast<int>(std::floor(fi));
    int j0 = static_cast<int>(std::floor(fj));
    int k0 = static_cast<int>(std::floor(fk));
    int i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1;

    amrex::Real ti = fi - i0;
    amrex::Real tj = fj - j0;
    amrex::Real tk = fk - k0;

    const amrex::Box& domain = geom.Domain();
    auto clamp = [](int v, int lo, int hi) {
        return std::max(lo, std::min(hi, v));
    };
    i0 = clamp(i0, domain.smallEnd(0), domain.bigEnd(0));
    i1 = clamp(i1, domain.smallEnd(0), domain.bigEnd(0));
    j0 = clamp(j0, domain.smallEnd(1), domain.bigEnd(1));
    j1 = clamp(j1, domain.smallEnd(1), domain.bigEnd(1));
    k0 = clamp(k0, domain.smallEnd(2), domain.bigEnd(2));
    k1 = clamp(k1, domain.smallEnd(2), domain.bigEnd(2));

    amrex::Real val = 0.0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        if (!bx.contains(amrex::IntVect(AMREX_D_DECL(i0, j0, k0)))) {
            continue;
        }
        const auto& arr = mf.const_array(mfi);
        auto get = [&](int ii, int jj, int kk) {
            return arr(
                clamp(ii, bx.smallEnd(0), bx.bigEnd(0)),
                clamp(jj, bx.smallEnd(1), bx.bigEnd(1)),
                clamp(kk, bx.smallEnd(2), bx.bigEnd(2)), comp);
        };
        val = (1 - ti) * (1 - tj) * (1 - tk) * get(i0, j0, k0) +
              ti * (1 - tj) * (1 - tk) * get(i1, j0, k0) +
              (1 - ti) * tj * (1 - tk) * get(i0, j1, k0) +
              ti * tj * (1 - tk) * get(i1, j1, k0) +
              (1 - ti) * (1 - tj) * tk * get(i0, j0, k1) +
              ti * (1 - tj) * tk * get(i1, j0, k1) +
              (1 - ti) * tj * tk * get(i0, j1, k1) +
              ti * tj * tk * get(i1, j1, k1);
        break;
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
        pp.query("n_points", p.n_sparger_points);
        pp.query("flow_rate", p.flow_rate);
        pp.query("bubble_diameter", p.d0);

        // Sparger positions — expected as flat list x0 y0 z0 x1 y1 z1 ...
        amrex::Vector<amrex::Real> pos_flat;
        if (pp.queryarr("positions", pos_flat) != 0) {
            p.sparger_x.resize(p.n_sparger_points);
            p.sparger_y.resize(p.n_sparger_points);
            p.sparger_z.resize(p.n_sparger_points);
            for (int i = 0; i < p.n_sparger_points; ++i) {
                p.sparger_x[i] = pos_flat[3 * i + 0];
                p.sparger_y[i] = pos_flat[3 * i + 1];
                p.sparger_z[i] = pos_flat[3 * i + 2];
            }
        } else {
            // Fallback: sparger positions file not provided — will be set to
            // default
            p.sparger_x.assign(p.n_sparger_points, 0.0);
            p.sparger_y.assign(p.n_sparger_points, 0.0);
            p.sparger_z.assign(p.n_sparger_points, 0.0);
        }
    }

    {
        amrex::ParmParse pp("bubble");
        pp.query("O2_molar_volume", p.O2_molar_volume);
        pp.query("O2_solubility", p.O2_solubility);
        pp.query("O2_initial_conc", p.O2_init_conc);
        pp.query("kL_coefficient", p.kL_coeff);
        pp.query("D_O2", p.D_O2);
        pp.query("rho_fluid", p.rho_fluid);
        pp.query("nu_fluid", p.nu_fluid);
        pp.query("surface_tension", p.surf_tension);
        pp.query("P_atm", p.P_atm);
        pp.query("g_grav", p.g_grav);
        pp.query("free_surface_z", p.free_surface_z);

        // Coalescence
        int ienb = p.enable_coalescence ? 1 : 0;
        pp.query("enable_coalescence", ienb);
        p.enable_coalescence = (ienb != 0);
        pp.query("coalescence_Re_crit", p.Re_a_crit);
        pp.query("coalescence_start_time", p.coal_start_time);
        pp.query("coalescence_interval", p.coal_interval);

        // Breakup
        int ienbk = p.enable_breakup ? 1 : 0;
        pp.query("enable_breakup", ienbk);
        p.enable_breakup = (ienbk != 0);
        std::string bmodel = "triangle";
        pp.query("breakup_model", bmodel);
        if (bmodel == "equal") {
            p.breakup_model = 1;
        } else if (bmodel == "random") {
            p.breakup_model = 2;
        } else {
            p.breakup_model = 0;
        } // triangle
        pp.query("min_diameter", p.min_diameter);
        pp.query("eps_max", p.eps_max);

        pp.query("stats_int", p.stats_int);
        pp.query("stats_file", p.stats_file);

        // Force-deposition safety mechanisms (June 2026).
        pp.query("force_cap_factor", p.force_cap_factor);
        pp.query("require_liquid_host", p.require_liquid_host);

        // Boyle's law / hydrostatic-pressure correction (June 2026).
        pp.query("boyle_law_enable", p.boyle_law_enable);

        amrex::Print() << "[bubble] force_cap_factor    = "
                       << p.force_cap_factor
                       << (p.force_cap_factor > 0.0
                               ? "  (|F_couple| <= cap·m_eff·|g|)"
                               : "  (DISABLED)")
                       << "\n";
        amrex::Print()
            << "[bubble] require_liquid_host = " << p.require_liquid_host
            << ((p.require_liquid_host != 0)
                    ? "  (skip force deposit when host cell != CELL_LIQUID)"
                    : "  (DISABLED — deposit at every host cell)")
            << "\n";
        amrex::Print()
            << "[bubble] boyle_law_enable    = " << p.boyle_law_enable
            << ((p.boyle_law_enable != 0)
                    ? "  (V_m = V_m,ref·P_atm/(P_atm + ρ_l·g·h); V_b, d, C_g "
                      "all vary with depth)"
                    : "  (DISABLED — V_m fixed; C_g = 1/V_m identically)")
            << "\n";
    }
}

// ============================================================================
// initialize
// ============================================================================
void BubbleManager::initialize(
    const amrex::Geometry& geom,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    const BubbleParams& params,
    bool append_stats)
{
    m_params = params;
    m_container.Define(geom, dm, ba);
    m_container
        .resizeData(); // resizes m_particles to finestLevel()+1 levels
                       // (Define() only sets GDB; constructors call
                       //  resizeData() automatically but Define() does not)
    // Per-hole residuals, pre-seeded to 1.0 so the first call always injects.
    m_injection_residuals.assign(m_params.n_sparger_points, 1.0);
    m_initialized = true;
    open_stats_file(append_stats);
    amrex::Print() << "[BubbleManager] Initialized. n_sparger_points = "
                   << m_params.n_sparger_points << "\n";
}

// ============================================================================
// inject_bubbles
//   Inject bubbles at sparger locations.
//   dt_phys: physical seconds in this step.
//   All injection is staged on rank 0; Redistribute() scatters to correct
//   ranks.
// ============================================================================
void BubbleManager::inject_bubbles(amrex::Real dt_phys)
{
    if (!m_initialized) {
        return;
    }
    if (m_params.sparger_x.empty()) {
        return;
    }

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
        if (static_cast<int>(m_injection_residuals[ih]) > 0) {
            any_inject = true;
        }
    }

    if (!any_inject) {
        // Still call Redistribute to ensure tile map is populated even on steps
        // where no bubbles are injected (e.g. first step when residual < 1).
        if (!m_particles_ever_injected) {
            return;
        } // m_particles not yet sized
        m_container.Redistribute();
        return;
    }

    // All injection on rank 0; Redistribute() scatters to correct ranks
    if (amrex::ParallelDescriptor::IOProcessor()) {
        auto& ptile = m_container.DefineAndReturnParticleTile(0, 0, 0);
        static std::mt19937 rng(42);
        std::uniform_real_distribution<amrex::Real> jitter(-0.3, 0.3);
        const amrex::Geometry& geom = m_container.Geom(0);
        const amrex::Real* plo = geom.ProbLo();
        const amrex::Real* phi = geom.ProbHi();

        for (int ih = 0; ih < m_params.n_sparger_points; ++ih) {
            int n_here = static_cast<int>(m_injection_residuals[ih]);
            m_injection_residuals[ih] -= static_cast<amrex::Real>(n_here);
            for (int ib = 0; ib < n_here; ++ib) {
                BubbleParticle p;
                p.id() = BubbleParticle::NextID();
                p.cpu() = amrex::ParallelDescriptor::MyProc();
                p.pos(0) = m_params.sparger_x[ih] + jitter(rng);
                p.pos(1) = m_params.sparger_y[ih] + jitter(rng);
                p.pos(2) = m_params.sparger_z[ih] + jitter(rng);
                // Clamp to valid geometry
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    p.pos(d) = std::max(
                        plo[d] + 0.5, std::min(phi[d] - 0.5, p.pos(d)));
                }
                p.rdata(bubble_idx::VX) = 0.0;
                p.rdata(bubble_idx::VY) = 0.0;
                p.rdata(bubble_idx::VZ) = 0.0;
                p.rdata(bubble_idx::DIAMETER) = m_params.d0;
                // Initial moles (Boyle-aware): legacy n_O2 = O2_init_conc · V0
                // implicitly assumes V_m = V_m,ref at the injection depth.
                // When Boyle is on, the bubble at sparger depth is compressed,
                // so the same geometric volume V0 holds more moles by a factor
                // of V_m,ref / V_m(h) = P(h) / P_atm.  Multiplying here makes
                // C_g_inject = n_O2/V0 = O2_init_conc · P(h)/P_atm, i.e. the
                // bubble starts at the local hydrostatic concentration.
                {
                    const amrex::Real Vm_inject =
                        v_m_at_depth(p.pos(2), m_params);
                    const amrex::Real depth_factor =
                        m_params.O2_molar_volume / Vm_inject; // = P(h)/P_atm
                    p.rdata(bubble_idx::N_O2) =
                        m_params.O2_init_conc * V0 * depth_factor;
                }
                p.rdata(bubble_idx::AX) = 0.0;
                p.rdata(bubble_idx::AY) = 0.0;
                p.rdata(bubble_idx::AZ) = 0.0;
                p.rdata(bubble_idx::BREAKUP_COOLDOWN) = 0.0;
                ptile.push_back(p);
            }
        }
    }
    m_container.Redistribute();
    m_particles_ever_injected = true;
}

// ============================================================================
// do_breakup
//   Hinze criterion: break if D > D_e = (sigma/rho_f)^0.6 * epsilon^(-0.4)
//   Epsilon is interpolated from the derived MultiFab (EPSILON_IDX),
//   which is computed in compute_derived() using the Smagorinsky SGS model.
//   LB units → SI: eps_SI = eps_LB * dx²/dt³.
// ============================================================================
void BubbleManager::do_breakup(
    [[maybe_unused]] const amrex::MultiFab& derived,
    [[maybe_unused]] const amrex::Geometry& geom)
{
    if (!m_params.enable_breakup) {
        return;
    }

    static std::mt19937 rng_brk(12345);
    std::uniform_real_distribution<amrex::Real> uniform_dist(0.0, 1.0);

    const amrex::Real sigma = m_params.surf_tension;
    const amrex::Real rho_f = m_params.rho_fluid;
    const amrex::Real dt = m_params.dt_phys;
    const amrex::Real nu_f = m_params.nu_fluid;

    amrex::Vector<BubbleParticle> new_bubbles;

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& pbox = kv.second;
            auto& aos = pbox.GetArrayOfStructs();
            for (auto& p : aos()) {
                if (!p.id().is_valid()) {
                    continue;
                }

                // Skip if bubble is still in breakup cooldown (Kolmogorov
                // timescale)
                if (p.rdata(bubble_idx::BREAKUP_COOLDOWN) > 0.0) {
                    continue;
                }

                // Read cached turbulent dissipation rate interpolated during
                // advance()
                const amrex::Real eps_SI = std::min(
                    amrex::Real(p.rdata(bubble_idx::EPS_CACHED)),
                    m_params.eps_max); // cap: prevent boundary-layer
                                       // over-fragmentation

                const amrex::Real d = p.rdata(bubble_idx::DIAMETER);

                // Skip if bubble is already at or below minimum allowed size
                if (d <= m_params.min_diameter) {
                    continue;
                }

                const amrex::Real D_e =
                    std::pow(sigma / rho_f, 0.6) * std::pow(eps_SI, -0.4);

                if (d <= D_e) {
                    continue;
                } // no breakup

                // Compute daughter sizes BEFORE invalidating the parent so we
                // can skip gracefully if either daughter falls below
                // min_diameter.
                amrex::Real fv = 0.5;
                if (m_params.breakup_model == 0) {
                    // Triangle distribution: lower=0, upper=0.5, mode=0.2
                    // Inverse CDF for triangle distribution
                    const amrex::Real u = uniform_dist(rng_brk);
                    const amrex::Real low = 0.0, high = 0.5, mode = 0.2;
                    const amrex::Real Fc = (mode - low) / (high - low);
                    if (u < Fc) {
                        fv = low + std::sqrt(u * (high - low) * (mode - low));
                    } else {
                        fv =
                            high -
                            std::sqrt((1.0 - u) * (high - low) * (high - mode));
                    }
                } else if (m_params.breakup_model == 2) {
                    fv = 0.5 * uniform_dist(rng_brk);
                } else {
                    fv = 0.5; // equal split
                }

                const amrex::Real V_total =
                    (amrex::Math::pi<amrex::Real>() / 6.0) * std::pow(d, 3.0);
                const amrex::Real V1 = fv * V_total;
                const amrex::Real V2 = (1.0 - fv) * V_total;
                const amrex::Real d1 =
                    std::cbrt(6.0 * V1 / amrex::Math::pi<amrex::Real>());
                const amrex::Real d2 =
                    std::cbrt(6.0 * V2 / amrex::Math::pi<amrex::Real>());
                const amrex::Real n1 = p.rdata(bubble_idx::N_O2) * fv;
                const amrex::Real n2 = p.rdata(bubble_idx::N_O2) * (1.0 - fv);

                // Skip breakup if either daughter would fall below
                // min_diameter. Prevents runaway fragmentation from near-wall
                // epsilon spikes.
                if (d1 < m_params.min_diameter || d2 < m_params.min_diameter) {
                    continue;
                }

                // All checks passed — invalidate the parent and create two
                // daughters
                p.id() = -1;

                // Kolmogorov timescale cooldown: τ_η = sqrt(ν/ε) [s]
                // Convert to LB steps: cooldown = τ_η / dt_phys
                const amrex::Real tau_eta = std::sqrt(nu_f / eps_SI);
                const amrex::Real cooldown_steps = std::ceil(tau_eta / dt);

                // Create daughter bubbles (add to new_bubbles list)
                for (int id = 0; id < 2; ++id) {
                    BubbleParticle daughter;
                    daughter.id() = BubbleParticle::NextID();
                    daughter.cpu() = amrex::ParallelDescriptor::MyProc();
                    daughter.pos(0) = p.pos(0);
                    daughter.pos(1) = p.pos(1);
                    daughter.pos(2) = p.pos(2);
                    daughter.rdata(bubble_idx::VX) = p.rdata(bubble_idx::VX);
                    daughter.rdata(bubble_idx::VY) = p.rdata(bubble_idx::VY);
                    daughter.rdata(bubble_idx::VZ) = p.rdata(bubble_idx::VZ);
                    daughter.rdata(bubble_idx::AX) = p.rdata(bubble_idx::AX);
                    daughter.rdata(bubble_idx::AY) = p.rdata(bubble_idx::AY);
                    daughter.rdata(bubble_idx::AZ) = p.rdata(bubble_idx::AZ);
                    daughter.rdata(bubble_idx::DIAMETER) = (id == 0) ? d1 : d2;
                    daughter.rdata(bubble_idx::N_O2) = (id == 0) ? n1 : n2;
                    daughter.rdata(bubble_idx::BREAKUP_COOLDOWN) =
                        cooldown_steps;
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
    if (!m_params.enable_coalescence) {
        return;
    }
    if (phys_time < m_params.coal_start_time) {
        return;
    }

    const amrex::Real nu_f = m_params.nu_fluid;
    const amrex::Real dx = m_params.dx_phys;

    // Collect all valid bubble data CPU-side
    struct BData
    {
        int lev;
        amrex::Long idx_in_box;
        amrex::Real x, y, z, vx, vy, vz, d, n_O2;
        bool valid = true;
    };
    amrex::Vector<BData> bvec;

    // Map from global index to (lev, grid, tile, local_index) for merging
    struct BRef
    {
        int lev;
        int grid;
        int tile;
        int local_idx;
    };
    amrex::Vector<BRef> brefs;

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& pbox = kv.second;
            auto& aos = pbox.GetArrayOfStructs();
            int grid = kv.first.first;
            int tile = kv.first.second;
            for (int li = 0; li < static_cast<int>(aos().size()); ++li) {
                const auto& p = aos()[li];
                if (!p.id().is_valid()) {
                    continue;
                }
                BData bd;
                bd.x = p.pos(0) * dx; // convert to metres
                bd.y = p.pos(1) * dx;
                bd.z = p.pos(2) * dx;
                bd.vx = p.rdata(bubble_idx::VX);
                bd.vy = p.rdata(bubble_idx::VY);
                bd.vz = p.rdata(bubble_idx::VZ);
                bd.d = p.rdata(bubble_idx::DIAMETER);
                bd.n_O2 = p.rdata(bubble_idx::N_O2);
                bd.valid = true;
                bvec.push_back(bd);
                BRef br{lev, grid, tile, li};
                brefs.push_back(br);
            }
        }
    }

    const int N = static_cast<int>(bvec.size());
    for (int i = 0; i < N; ++i) {
        if (!bvec[i].valid) {
            continue;
        }
        for (int j = i + 1; j < N; ++j) {
            if (!bvec[j].valid) {
                continue;
            }

            const amrex::Real dxi = bvec[i].d;
            const amrex::Real dxj = bvec[j].d;
            const amrex::Real d_h =
                2.0 * dxi * dxj / (dxi + dxj); // harmonic diameter

            // Centre-to-centre distance
            const amrex::Real rx = bvec[i].x - bvec[j].x;
            const amrex::Real ry = bvec[i].y - bvec[j].y;
            const amrex::Real rz = bvec[i].z - bvec[j].z;
            const amrex::Real dist = std::sqrt(rx * rx + ry * ry + rz * rz);

            // Only consider overlapping or very close bubbles
            if (dist > 0.5 * (dxi + dxj) * 1.5) {
                continue;
            }

            // Approach velocity (component along line of centres)
            const amrex::Real dvx = bvec[i].vx - bvec[j].vx;
            const amrex::Real dvy = bvec[i].vy - bvec[j].vy;
            const amrex::Real dvz = bvec[i].vz - bvec[j].vz;
            const amrex::Real dv_mag =
                std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

            const amrex::Real Re_a = dv_mag * d_h / nu_f;

            if (Re_a > m_params.Re_a_crit) {
                // Coalescence: conserve volume and momentum
                const amrex::Real V_i =
                    (amrex::Math::pi<amrex::Real>() / 6.0) * dxi * dxi * dxi;
                const amrex::Real V_j =
                    (amrex::Math::pi<amrex::Real>() / 6.0) * dxj * dxj * dxj;
                const amrex::Real V_new = V_i + V_j;
                const amrex::Real d_new =
                    std::cbrt(6.0 * V_new / amrex::Math::pi<amrex::Real>());
                // Momentum conservation (using m ≈ 0, so mass average is
                // unweighted)
                const amrex::Real vx_new = 0.5 * (bvec[i].vx + bvec[j].vx);
                const amrex::Real vy_new = 0.5 * (bvec[i].vy + bvec[j].vy);
                const amrex::Real vz_new = 0.5 * (bvec[i].vz + bvec[j].vz);
                const amrex::Real n_new = bvec[i].n_O2 + bvec[j].n_O2;
                // Position: centre-of-volume
                const amrex::Real x_new = 0.5 * (bvec[i].x + bvec[j].x);
                const amrex::Real y_new = 0.5 * (bvec[i].y + bvec[j].y);
                const amrex::Real z_new = 0.5 * (bvec[i].z + bvec[j].z);

                // Update particle i with coalesced data
                auto& br_i = brefs[i];
                auto& aos_i =
                    m_container.GetParticles(br_i.lev)[{br_i.grid, br_i.tile}]
                        .GetArrayOfStructs();
                auto& pi = aos_i()[br_i.local_idx];
                pi.pos(0) = x_new / dx;
                pi.pos(1) = y_new / dx;
                pi.pos(2) = z_new / dx;
                pi.rdata(bubble_idx::VX) = vx_new;
                pi.rdata(bubble_idx::VY) = vy_new;
                pi.rdata(bubble_idx::VZ) = vz_new;
                pi.rdata(bubble_idx::DIAMETER) = d_new;
                pi.rdata(bubble_idx::N_O2) = n_new;

                // Invalidate particle j
                auto& br_j = brefs[j];
                auto& aos_j =
                    m_container.GetParticles(br_j.lev)[{br_j.grid, br_j.tile}]
                        .GetArrayOfStructs();
                aos_j()[br_j.local_idx].id() = -1;

                bvec[i].d = d_new;
                bvec[i].vx = vx_new;
                bvec[i].vy = vy_new;
                bvec[i].vz = vz_new;
                bvec[i].n_O2 = n_new;
                bvec[i].x = x_new;
                bvec[i].y = y_new;
                bvec[i].z = z_new;
                bvec[j].valid = false;
            }
            // else: elastic bounce — for now do nothing (simple model)
        }
    }
    m_container.Redistribute();
}

// ============================================================================
// advance  — fused GPU per-step driver
// ============================================================================
void BubbleManager::advance(
    [[maybe_unused]] amrex::Real dt,
    const amrex::MultiFab& macrodata,
    const amrex::MultiFab& derived,
    const amrex::MultiFab& o2_conc_mf,
    const amrex::Geometry& geom,
    amrex::MultiFab& fluid_force_mf,
    amrex::MultiFab& o2_src_mf,
    [[maybe_unused]] amrex::Real phys_time,
    const amrex::MultiFab* phi_mf,
    const amrex::iMultiFab* cell_type_mf)
{
    if (!m_initialized) {
        return;
    }
    if (!m_particles_ever_injected) {
        return;
    }

    // We no longer make host-pinned copies because Particles live in
    // Managed/Device memory implicitly via DefaultAllocator.

    // Reset output MFs to zero natively on device
    fluid_force_mf.setVal(0.0);
    o2_src_mf.setVal(0.0);

    const amrex::Real dx_phys = m_params.dx_phys;

    // Set up iteration
    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {

        const auto& geom_lev = m_container.Geom(lev);
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
            geom_lev.ProbLoArray();
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_arr =
            geom_lev.CellSizeArray();
        amrex::Box domain = geom_lev.Domain();

        for (BubbleContainer::ParIterType pti(m_container, lev); pti.isValid();
             ++pti) {
            auto& aos = pti.GetArrayOfStructs();
            auto* pData = aos().data();
            int num_part = pti.numParticles();
            if (num_part == 0) {
                continue;
            }

            auto md_arr = macrodata.const_array(pti);
            auto der_arr = derived.const_array(pti);
            auto o2_arr = o2_conc_mf.const_array(pti);
            auto force_arr = fluid_force_mf.array(pti);
            auto src_arr = o2_src_mf.array(pti);

            bool has_isf = (cell_type_mf != nullptr);
            auto isf_arr = has_isf ? cell_type_mf->const_array(pti)
                                   : amrex::Array4<const int>{};
            bool has_phi = (phi_mf != nullptr);
            auto phi_arr = has_phi ? phi_mf->const_array(pti)
                                   : amrex::Array4<const amrex::Real>{};

            // Populate the POD scalar-subset of m_params captured by the
            // lambda below.  See BubbleParamsScalar (namespace scope, above)
            // for the rationale.
            const BubbleParamsScalar prms{
                m_params.dx_phys,         m_params.dt_phys,
                m_params.nu_fluid,        m_params.rho_fluid,
                m_params.g_grav,          m_params.force_cap_factor,
                m_params.C_ref,           m_params.D_O2,
                m_params.kL_coeff,        m_params.O2_solubility,
                m_params.free_surface_z,  m_params.P_atm,
                m_params.O2_molar_volume, m_params.require_liquid_host,
                m_params.boyle_law_enable};

            amrex::ParallelFor(num_part, [=] AMREX_GPU_DEVICE(int i) {
                auto& p = pData[i];
                if (!p.id().is_valid()) {
                    return;
                }

                amrex::Real px = p.pos(0);
                amrex::Real py = p.pos(1);
                amrex::Real pz = p.pos(2);

                // Bubble geometry from the current stored diameter.  Diameter
                // is updated at the END of this kernel from the post-update
                // moles n_O2_new via V_b = n_O2 · v_m_at_depth(z, prms), so
                // Boyle's-law depth scaling enters the next step's d / V_b /
                // C_g automatically and self-consistently (no separate inline
                // compression block needed).
                const amrex::Real d_new = p.rdata(bubble_idx::DIAMETER);
                const amrex::Real r = 0.5 * d_new;
                const amrex::Real Vb =
                    (4.0 / 3.0) * amrex::Math::pi<amrex::Real>() * r * r * r;
                const amrex::Real Ab = amrex::Math::pi<amrex::Real>() * r * r;

                // 2. Interpolate fluid velocity -> Bubble Drag
                amrex::Real ufx_lb = trilinear_interp_device(
                    md_arr, constants::VELX_IDX, prob_lo.data(), dx_arr.data(),
                    px, py, pz);
                amrex::Real ufy_lb = trilinear_interp_device(
                    md_arr, constants::VELY_IDX, prob_lo.data(), dx_arr.data(),
                    px, py, pz);
                amrex::Real ufz_lb = trilinear_interp_device(
                    md_arr, constants::VELZ_IDX, prob_lo.data(), dx_arr.data(),
                    px, py, pz);

                const amrex::Real ufx = ufx_lb * dx_phys / prms.dt_phys;
                const amrex::Real ufy = ufy_lb * dx_phys / prms.dt_phys;
                const amrex::Real ufz = ufz_lb * dx_phys / prms.dt_phys;

                amrex::Real vbx = p.rdata(bubble_idx::VX);
                amrex::Real vby = p.rdata(bubble_idx::VY);
                amrex::Real vbz = p.rdata(bubble_idx::VZ);

                amrex::Real urx = vbx - ufx;
                amrex::Real ury = vby - ufy;
                amrex::Real urz = vbz - ufz;
                amrex::Real ur_mag =
                    std::sqrt(urx * urx + ury * ury + urz * urz);

                amrex::Real Re_b =
                    (ur_mag > 0.0) ? (ur_mag * d_new / prms.nu_fluid) : 0.0;
                amrex::Real V_cell = dx_phys * dx_phys * dx_phys;
                amrex::Real eps_p = amrex::min(Vb / V_cell, amrex::Real(0.99));
                amrex::Real Cd = bubble_drag_cd(Re_b, eps_p);

                amrex::Real fd_factor =
                    -0.5 * Cd * prms.rho_fluid * Ab * ur_mag;
                amrex::Real FDx = fd_factor * urx;
                amrex::Real FDy = fd_factor * ury;
                amrex::Real FDz = fd_factor * urz;

                amrex::Real FBx = 0.0, FBy = 0.0,
                            FBz = prms.rho_fluid * Vb * prms.g_grav;
                amrex::Real m_eff =
                    constants::C_ADDED_MASS * prms.rho_fluid * Vb;

                amrex::Real ax_new = (FBx + FDx) / m_eff;
                amrex::Real ay_new = (FBy + FDy) / m_eff;
                amrex::Real az_new = (FBz + FDz) / m_eff;

                amrex::Real FAx = m_eff * ax_new;
                amrex::Real FAy = m_eff * ay_new;
                amrex::Real FAz = m_eff * az_new;

                p.rdata(bubble_idx::AX) = ax_new;
                p.rdata(bubble_idx::AY) = ay_new;
                p.rdata(bubble_idx::AZ) = az_new;

                // 3. Deposit drag+added mass coupling force back to fluid
                // LB force = SI Force / (rho * dx^4 / dt^2)
                amrex::Real conv = (prms.dt_phys * prms.dt_phys) /
                                   (prms.rho_fluid * V_cell * dx_phys);

                // -------- Safety mechanism 1: per-bubble force cap --------
                // Bound |F_couple| at cap_factor · m_eff · |g_grav|.  This
                // prevents a single anomalous fluid-velocity sample (e.g.
                // bubble briefly inside a stale interpolation stencil that
                // includes gas-side or post-impeller transient cells) from
                // injecting unbounded forces into the LB grid.  At
                // cap=50, the ceiling is ~50·m_eff·g — well above the
                // ~10·m_eff·g peak observed in physically clean runs.
                amrex::Real F_couple_x = FAx + FDx;
                amrex::Real F_couple_y = FAy + FDy;
                amrex::Real F_couple_z = FAz + FDz;
                if (prms.force_cap_factor > 0.0) {
                    const amrex::Real F_max =
                        prms.force_cap_factor * m_eff * std::abs(prms.g_grav);
                    const amrex::Real Fmag = std::sqrt(
                        F_couple_x * F_couple_x + F_couple_y * F_couple_y +
                        F_couple_z * F_couple_z);
                    if (Fmag > F_max && Fmag > 0.0) {
                        const amrex::Real s = F_max / Fmag;
                        F_couple_x *= s;
                        F_couple_y *= s;
                        F_couple_z *= s;
                    }
                    // Also catch NaN/Inf — replace with zero so the deposit
                    // is harmless rather than poisoning the fluid grid.
                    if (!std::isfinite(F_couple_x)) {
                        F_couple_x = 0.0;
                    }
                    if (!std::isfinite(F_couple_y)) {
                        F_couple_y = 0.0;
                    }
                    if (!std::isfinite(F_couple_z)) {
                        F_couple_z = 0.0;
                    }
                }

                amrex::Real Ffx = -F_couple_x * conv;
                amrex::Real Ffy = -F_couple_y * conv;
                amrex::Real Ffz = -F_couple_z * conv;

                int ci = static_cast<int>(
                    amrex::Math::floor((px - prob_lo[0]) / dx_arr[0]));
                int cj = static_cast<int>(
                    amrex::Math::floor((py - prob_lo[1]) / dx_arr[1]));
                int ck = static_cast<int>(
                    amrex::Math::floor((pz - prob_lo[2]) / dx_arr[2]));
                if (domain.contains(amrex::IntVect(AMREX_D_DECL(ci, cj, ck)))) {
                    // -------- Safety mechanism 2: liquid-host gating --------
                    // Drag/added-mass forces computed from interpolated fluid
                    // velocity are only physically meaningful when the host
                    // cell is bulk liquid.  When the bubble passes through
                    // CELL_INTERFACE (free-surface band) or briefly sits in
                    // CELL_GAS (above-surface) or CELL_SOLID (impeller-swept
                    // volume), the interpolation stencil pulls in stale or
                    // ABB-equilibrium populations and produces unphysical
                    // u_fluid values.  Skip the deposit in those cases; the
                    // bubble's own kinematic state still updates.
                    bool deposit = true;
                    if ((prms.require_liquid_host != 0) && has_isf) {
                        const int ct_host = isf_arr(ci, cj, ck, 0);
                        if (ct_host != lbm::constants::CELL_LIQUID) {
                            deposit = false;
                        }
                    }
                    if (deposit) {
                        amrex::Gpu::Atomic::AddNoRet(
                            &force_arr(ci, cj, ck, 0), Ffx);
                        amrex::Gpu::Atomic::AddNoRet(
                            &force_arr(ci, cj, ck, 1), Ffy);
                        amrex::Gpu::Atomic::AddNoRet(
                            &force_arr(ci, cj, ck, 2), Ffz);
                    }
                }

                // 4. O2 Volumetric Source
                amrex::Real n_O2 = p.rdata(bubble_idx::N_O2);
                amrex::Real C_g_i = (Vb > 0.0) ? (n_O2 / Vb) : 0.0;

                amrex::Real C_f_lb_raw = trilinear_interp_device(
                    o2_arr, 0, prob_lo.data(), dx_arr.data(), px, py, pz);
                amrex::Real C_f_i =
                    amrex::max(C_f_lb_raw, amrex::Real(0.0)) * prms.C_ref;

                amrex::Real eps_lb = trilinear_interp_device(
                    der_arr, constants::EPSILON_IDX, prob_lo.data(),
                    dx_arr.data(), px, py, pz);
                amrex::Real eps_conv =
                    dx_phys * dx_phys /
                    (prms.dt_phys * prms.dt_phys * prms.dt_phys);
                amrex::Real eps_SI =
                    amrex::max(eps_lb * eps_conv, amrex::Real(1.0e-10));
                p.rdata(bubble_idx::EPS_CACHED) = eps_SI;

                amrex::Real Sc = prms.nu_fluid / prms.D_O2;
                amrex::Real k_L = prms.kL_coeff *
                                  std::pow(eps_SI * prms.nu_fluid, 0.25) *
                                  std::pow(Sc, -0.5);
                amrex::Real S_i = prms.O2_solubility;

                amrex::Real dn_i =
                    k_L * (amrex::Math::pi<amrex::Real>() * d_new * d_new) *
                    (S_i * C_g_i - C_f_i);

                p.rdata(bubble_idx::DN_I) = dn_i; // Cache for stats

                amrex::Real n_O2_new =
                    amrex::max(n_O2 - dn_i * prms.dt_phys, amrex::Real(0.0));
                p.rdata(bubble_idx::N_O2) = n_O2_new;

                // Recompute V_b from moles using Boyle-aware molar volume.
                // When prms.boyle_law_enable = 0, this collapses to the legacy
                //   V_b = n_O2 · O2_molar_volume
                // (isobaric ideal gas at STP).  When on, V_m shrinks with
                // depth so a fixed moles → smaller V_b for deeper bubbles.
                const amrex::Real Vm_now = v_m_at_depth(p.pos(2), prms);
                amrex::Real Vb_new = n_O2_new * Vm_now;
                if (Vb_new > 0.0) {
                    amrex::Real d_new_new = std::cbrt(
                        6.0 * Vb_new / amrex::Math::pi<amrex::Real>());
                    if (d_new_new < 1.0e-5) {
                        p.id() = -1;
                    } else {
                        p.rdata(bubble_idx::DIAMETER) = d_new_new;
                    }
                } else {
                    p.id() = -1;
                }

                amrex::Real src_rate = dn_i / V_cell;
                if (std::isfinite(src_rate) && std::abs(src_rate) < 1.0e4) {
                    if (domain.contains(
                            amrex::IntVect(AMREX_D_DECL(ci, cj, ck)))) {
                        amrex::Gpu::Atomic::AddNoRet(
                            &src_arr(ci, cj, ck, 0), src_rate);
                    }
                }

                // 5. Velocity Verlet Position Update + Solid Body Collisions
                if (p.id().is_valid()) {
                    amrex::Real old_x = px;
                    amrex::Real old_y = py;
                    amrex::Real old_z = pz;

                    amrex::Real new_px =
                        px + (vbx * prms.dt_phys +
                              0.5 * ax_new * prms.dt_phys * prms.dt_phys) /
                                 dx_phys;
                    amrex::Real new_py =
                        py + (vby * prms.dt_phys +
                              0.5 * ay_new * prms.dt_phys * prms.dt_phys) /
                                 dx_phys;
                    amrex::Real new_pz =
                        pz + (vbz * prms.dt_phys +
                              0.5 * az_new * prms.dt_phys * prms.dt_phys) /
                                 dx_phys;

                    p.rdata(bubble_idx::VX) = vbx + ax_new * prms.dt_phys;
                    p.rdata(bubble_idx::VY) = vby + ay_new * prms.dt_phys;
                    p.rdata(bubble_idx::VZ) = vbz + az_new * prms.dt_phys;

                    if (has_isf) {
                        int nci = static_cast<int>(amrex::Math::floor(
                            (new_px - prob_lo[0]) / dx_arr[0]));
                        int ncj = static_cast<int>(amrex::Math::floor(
                            (new_py - prob_lo[1]) / dx_arr[1]));
                        [[maybe_unused]] int nck =
                            static_cast<int>(amrex::Math::floor(
                                (new_pz - prob_lo[2]) / dx_arr[2]));
                        amrex::IntVect niv(AMREX_D_DECL(nci, ncj, nck));

                        // Block bubbles entering MOVING-IMPELLER or
                        // STATIONARY-WALL cells only. GAS cells (free surface
                        // headspace) are NOT solid; bubbles rising through the
                        // surface are handled by the phi-exit check below.
                        // Treating GAS as solid here was responsible for
                        // bubbles freezing at the free surface and visually
                        // "passing into" the tank walls.
                        bool new_in_solid =
                            (!domain.contains(niv) ||
                             isf_arr(niv, 0) == lbm::constants::CELL_SOLID);
                        if (new_in_solid) {
                            amrex::IntVect oiv(AMREX_D_DECL(ci, cj, ck));
                            bool old_in_solid =
                                (!domain.contains(oiv) ||
                                 isf_arr(oiv, 0) == lbm::constants::CELL_SOLID);
                            if (!old_in_solid) {
                                new_px = old_x;
                                new_py = old_y;
                                new_pz = old_z;
                                p.rdata(bubble_idx::VX) = 0;
                                p.rdata(bubble_idx::VY) = 0;
                                p.rdata(bubble_idx::VZ) = 0;
                                p.rdata(bubble_idx::AX) = 0;
                                p.rdata(bubble_idx::AY) = 0;
                                p.rdata(bubble_idx::AZ) = 0;
                            } else {
                                bool found = false;
                                for (int d_i = -1; d_i <= 1 && !found; ++d_i) {
                                    for (int d_j = -1; d_j <= 1 && !found;
                                         ++d_j) {
                                        for (int d_k = -1; d_k <= 1 && !found;
                                             ++d_k) {
                                            amrex::IntVect nniv(AMREX_D_DECL(
                                                ci + d_i, cj + d_j, ck + d_k));
                                            if (domain.contains(nniv) &&
                                                isf_arr(nniv, 0) !=
                                                    lbm::constants::
                                                        CELL_SOLID) {
                                                new_px = (nniv[0] + 0.5) *
                                                             dx_arr[0] +
                                                         prob_lo[0];
                                                new_py = (nniv[1] + 0.5) *
                                                             dx_arr[1] +
                                                         prob_lo[1];
                                                new_pz = AMREX_D_PICK(
                                                    amrex::Real(0.0),
                                                    amrex::Real(0.0),
                                                    (nniv[2] + 0.5) *
                                                            dx_arr[2] +
                                                        prob_lo[2]);
                                                p.rdata(bubble_idx::VX) = 0;
                                                p.rdata(bubble_idx::VY) = 0;
                                                p.rdata(bubble_idx::VZ) = 0;
                                                p.rdata(bubble_idx::AX) = 0;
                                                p.rdata(bubble_idx::AY) = 0;
                                                p.rdata(bubble_idx::AZ) = 0;
                                                found = true;
                                            }
                                        }
                                    }
                                }
                                if (!found) {
                                    p.id() = -1;
                                }
                            }
                        }
                    }
                    p.pos(0) = new_px;
                    p.pos(1) = new_py;
                    p.pos(2) = new_pz;

                    // 6. Free Surface Outgassing Exit
                    if (has_phi) {
                        int fsi = static_cast<int>(amrex::Math::floor(
                            (new_px - prob_lo[0]) / dx_arr[0]));
                        int fsj = static_cast<int>(amrex::Math::floor(
                            (new_py - prob_lo[1]) / dx_arr[1]));
                        [[maybe_unused]] int fsk =
                            static_cast<int>(amrex::Math::floor(
                                (new_pz - prob_lo[2]) / dx_arr[2]));
                        amrex::IntVect fsiv(AMREX_D_DECL(fsi, fsj, fsk));
                        if (domain.contains(fsiv) && phi_arr(fsiv, 0) < 0.5) {
                            p.id() = -1;
                        }
                    } else {
                        if (new_pz >= prms.free_surface_z) {
                            p.id() = -1;
                        }
                    }
                }
            });
        }
    }

    // Sort and compact
    m_container.Redistribute();

    // 7. AMReX CPU fallback loops for collision/breakup via Managed Memory
    // These run natively on host via the unified memory wrapper tracking
    amrex::Gpu::synchronize();

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& aos = kv.second.GetArrayOfStructs();
            for (auto& p : aos()) {
                if (!p.id().is_valid()) {
                    continue;
                }
                amrex::ParticleReal& cd = p.rdata(bubble_idx::BREAKUP_COOLDOWN);
                if (cd > 0.0) {
                    cd -= 1.0;
                }
            }
        }
    }
    do_breakup(derived, geom);
}
// ============================================================================
// write_stats
// ============================================================================
void BubbleManager::write_stats(int step, amrex::Real phys_time)
{
    // Backward-compatible entry point: emit NaN for the C_L / V_liq fields
    // when the caller hasn't computed them (e.g. early init / shutdown
    // diagnostics).  Distinguishable from a true zero in the CSV.
    write_stats(
        step, phys_time, std::numeric_limits<amrex::Real>::quiet_NaN(),
        std::numeric_limits<amrex::Real>::quiet_NaN());
}

void BubbleManager::write_stats(
    int step,
    amrex::Real phys_time,
    amrex::Real C_L_mol_m3,
    amrex::Real V_liq_m3)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    // Gather statistics across all particles
    int n_bub = 0;
    amrex::Real d_sum = 0.0;
    amrex::Real d_min = 1.0e30;
    amrex::Real d_max = 0.0;
    amrex::Real n_O2_total = 0.0;
    amrex::Real dn_total = 0.0;

    for (int lev = 0; lev <= m_container.finestLevel(); ++lev) {
        for (auto& kv : m_container.GetParticles(lev)) {
            auto& aos = kv.second.GetArrayOfStructs();
            for (auto& p : aos()) {
                if (!p.id().is_valid()) {
                    continue;
                }
                ++n_bub;
                const amrex::Real d = p.rdata(bubble_idx::DIAMETER);
                d_sum += d;
                d_min = std::min(d_min, d);
                d_max = std::max(d_max, d);
                n_O2_total += p.rdata(bubble_idx::N_O2);
                dn_total += p.rdata(bubble_idx::DN_I);
            }
        }
    }

    const amrex::Real d_mean = (n_bub > 0) ? (d_sum / n_bub) : 0.0;

    if (!m_stats_stream.is_open()) {
        return;
    }
    m_stats_stream << std::fixed << std::setprecision(6) << step << ","
                   << phys_time << "," << n_bub << "," << d_mean * 1000.0
                   << "," // mm
                   << d_min * 1000.0 << "," << d_max * 1000.0 << ","
                   << n_O2_total << "," << dn_total << "," << std::scientific
                   << std::setprecision(6) << C_L_mol_m3
                   << ","               // volume-averaged dissolved O₂ [mol/m³]
                   << V_liq_m3 << "\n"; // liquid volume LIQ + φ·IFC [m³]
    m_stats_stream.flush();
}

// ============================================================================
// Stats file helpers
// ============================================================================
void BubbleManager::open_stats_file(bool append)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }
    if (append) {
        // Detect whether the file already has a header.  When the working
        // directory was wiped (or the previous CSV moved out) the restart
        // path still requests append=true so existing per-bubble data is
        // not clobbered, but the destination file is empty — we must write
        // the header in that case, otherwise the new run leaves an
        // un-parseable, header-less CSV.
        std::ifstream probe(m_params.stats_file, std::ios::binary);
        bool file_is_empty = true;
        if (probe.is_open()) {
            file_is_empty = (probe.peek() == std::ifstream::traits_type::eof());
        }
        m_stats_stream.open(m_params.stats_file, std::ios::out | std::ios::app);
        if (file_is_empty) {
            m_stats_stream
                << "step,phys_time_s,n_bubbles,d_mean_mm,"
                   "d_min_mm,d_max_mm,n_O2_total_mol,dn_O2_step_mol_per_s,"
                   "C_L_mol_m3,V_liq_m3\n";
        }
    } else {
        m_stats_stream.open(
            m_params.stats_file, std::ios::out | std::ios::trunc);
        m_stats_stream
            << "step,phys_time_s,n_bubbles,d_mean_mm,"
               "d_min_mm,d_max_mm,n_O2_total_mol,dn_O2_step_mol_per_s,"
               "C_L_mol_m3,V_liq_m3\n";
    }
}

void BubbleManager::close_stats_file()
{
    if (m_stats_stream.is_open()) {
        m_stats_stream.close();
    }
}

// ============================================================================
// Checkpoint / Restart
// ============================================================================
void BubbleManager::checkpoint(
    const std::string& dir, const std::string& name) const
{
    if (m_initialized && m_particles_ever_injected) {
        m_container.Checkpoint(dir, name);
    }
}

void BubbleManager::restart(const std::string& dir, const std::string& name)
{
    if (m_initialized) {
        m_container.Restart(dir, name);
        // Ensure m_particles_ever_injected correctly turns on if loaded
        if (m_container.TotalNumberOfParticles() > 0) {
            m_particles_ever_injected = true;
        }
    }
}

} // namespace lbm
