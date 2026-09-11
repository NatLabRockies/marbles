#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <AMReX_REAL.H>

namespace lbm::constants {

static constexpr int N_MACRO_STATES = 19;
static constexpr int RHO_IDX = 0;
static constexpr int VELX_IDX = 1;
static constexpr int VELY_IDX = 2;
static constexpr int VELZ_IDX = 3;
static constexpr int VMAG_IDX = 4;

static constexpr int TWO_RHO_E_IDX = 5;
static constexpr int Q_CORR_X_IDX = 6;
static constexpr int Q_CORR_Y_IDX = 7;
static constexpr int Q_CORR_Z_IDX = 8;

static constexpr int PXX_IDX = 9;
static constexpr int PYY_IDX = 10;
static constexpr int PZZ_IDX = 11;
static constexpr int PXY_IDX = 12;
static constexpr int PXZ_IDX = 13;
static constexpr int PYZ_IDX = 14;

static constexpr int QX_IDX = 15;
static constexpr int QY_IDX = 16;
static constexpr int QZ_IDX = 17;

static constexpr int TEMPERATURE_IDX = 18;

#if AMREX_SPACEDIM == 2
static constexpr int N_MICRO_STATES = 9;
#else
static constexpr int N_MICRO_STATES = 27;
#endif

static constexpr int N_DERIVED = 8;
static constexpr int VORTX_IDX = 0;
static constexpr int VORTY_IDX = 1;
static constexpr int VORTZ_IDX = 2;
static constexpr int VORTM_IDX = 3;

static constexpr int D_Q_CORR_X_IDX = 4;
static constexpr int D_Q_CORR_Y_IDX = 5;
static constexpr int D_Q_CORR_Z_IDX = 6;

// Energy dissipation rate ε (LB units: dx²/dt³) stored in derived[lev]
static constexpr int EPSILON_IDX = 7;

static constexpr int N_IS_FLUID = 4;
// component indices for m_is_fluid (iMultiFab)
static constexpr int IS_FLUID_IDX = 0;               // 0 = fluid/solid
static constexpr int EB_BOUNDARY_IDX = 1;            // 1 = EB boundary flag
static constexpr int IS_FLUID_SIDE_IDX = 2;          // 2 = fluid-side boundary
static constexpr int IS_FLUID_SIDE_BOUNDARY_IDX = 3; // 3 = next layer boundary

// ---------------------------------------------------------------
// FSLBM cell types (stored in m_cell_type iMultiFab, 1 component)
// Körner et al. (2005), Schwarzmeier et al. (2023, JCP 473 111753)
// ---------------------------------------------------------------
static constexpr int CELL_SOLID = 0; // EB or moving body — no PDFs
static constexpr int CELL_GAS = 1;   // headspace — no PDFs, no LBM update
static constexpr int CELL_INTERFACE =
    2; // free-surface cell — has PDFs, fill ∈ (0,1)
static constexpr int CELL_LIQUID = 3; // fully liquid — has PDFs, fill = 1

// Conversion thresholds for fill level φ
static constexpr amrex::Real FSLBM_PHI_LO = amrex::Real(1.0e-4); // below → GAS
static constexpr amrex::Real FSLBM_PHI_HI =
    amrex::Real(1.0) - amrex::Real(1.0e-4); // above → LIQUID

static constexpr amrex::Real ROOT3 =
    static_cast<amrex::Real>(1.732050807568877293527446341);

static constexpr amrex::Real LOW_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1e-10);

static constexpr amrex::Real SMALL_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::epsilon() * 1e10);

static constexpr amrex::Real LARGE_NUM =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::max() * 1e-10);

static constexpr int DATWIDTH = 24;
static constexpr int DATPRECISION = 16;

static constexpr amrex::Real AIR_MOLAR_MASS = 28.96;

// ---------------------------------------------------------------
// Physical constants for two-phase kLa mass transfer
// ---------------------------------------------------------------
static constexpr amrex::Real O2_MOLAR_VOLUME = 0.0224; // m³/mol at STP
static constexpr amrex::Real O2_SOLUBILITY = 0.032;    // dimensionless (Bunsen)
static constexpr amrex::Real O2_DIFFUSIVITY = 2.0e-9;  // m²/s in water at 25°C
static constexpr amrex::Real O2_INIT_CONC_GAS = 44.6;  // mol/m³ pure O2 at STP
static constexpr amrex::Real KL_COEFFICIENT = 0.301; // C in Kawase k_L formula
static constexpr amrex::Real SURF_TENSION_AIRWATER =
    0.072;                                       // N/m (air-water at 25°C)
static constexpr amrex::Real RHO_WATER = 1000.0; // kg/m³
static constexpr amrex::Real NU_WATER = 1.0e-6;  // m²/s kinematic viscosity
static constexpr amrex::Real P_ATM = 101325.0;   // Pa
static constexpr amrex::Real G_GRAV = 9.81; // m/s² gravitational acceleration
static constexpr amrex::Real C_ADDED_MASS =
    0.5; // added-mass coefficient (sphere)
static constexpr amrex::Real SMAGORINSKY_CS =
    0.1; // Smagorinsky constant (calibrated)

} // namespace lbm::constants

#endif
