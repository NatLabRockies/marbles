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

static constexpr int N_DERIVED = 7;
static constexpr int VORTX_IDX = 0;
static constexpr int VORTY_IDX = 1;
static constexpr int VORTZ_IDX = 2;
static constexpr int VORTM_IDX = 3;

static constexpr int D_Q_CORR_X_IDX = 4;
static constexpr int D_Q_CORR_Y_IDX = 5;
static constexpr int D_Q_CORR_Z_IDX = 6;

static constexpr int N_IS_FLUID = 3;

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

} // namespace lbm::constants

#endif
