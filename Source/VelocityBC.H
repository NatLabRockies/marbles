#ifndef VELOCITYBC_H
#define VELOCITYBC_H

#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include "Constants.H"

namespace lbm::bc {

struct NoOp
{
    struct DeviceOp
    {

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& /*iv*/,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& /*geom*/,
            amrex::Real& /*rho*/,
            amrex::RealVect& /*vel*/,
            amrex::Real& /*R*/,
            amrex::Real& /*T*/,
            amrex::Real& /*gamma*/) const
        {}
    };

    using DeviceType = DeviceOp;

    explicit NoOp(const std::string& prefix = "");

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct Constant
{
    struct DeviceOp
    {
        int dir{1};
        amrex::Real u0{0.0};

        amrex::Real initial_density = 1.0;
        amrex::Real initial_temperature = 1.0 / 3.0;
        amrex::Real adiabatic_exponent = 5.0 / 3.0;
        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;
        amrex::Real Mach_ref{0.005};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& /*iv*/,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& /*geom*/,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& R,
            amrex::Real& T,
            amrex::Real& gamma) const
        {
            rho = initial_density;
            vel[dir] = u0;
            R = (R_u / m_bar);
            T = initial_temperature;
            gamma = adiabatic_exponent;
        }
    };

    using DeviceType = DeviceOp;

    explicit Constant(const std::string& prefix = "velocity_bc_constant");

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct Channel
{
    struct DeviceOp
    {
        amrex::Real u_ref{1.0};

        amrex::Real initial_density = 1.0;
        amrex::Real initial_temperature = 1.0 / 3.0;
        amrex::Real adiabatic_exponent = 5.0 / 3.0;
        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;
        amrex::Real Mach_ref{0.005};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& R,
            amrex::Real& T,
            amrex::Real& gamma) const
        {
            rho = initial_density;
            const int* domhi = geom.Domain().hiVect();
            const amrex::Real c1 = iv[1] * (domhi[1] - iv[1]);
            const amrex::Real c2 =
                AMREX_D_PICK(, c1, iv[2] * (domhi[2] - iv[2]));
            vel[0] = 16.0 * u_ref * c1 * c2 / std::pow((domhi[1] + 1), 4);
            R = (R_u / m_bar);
            T = initial_temperature;
            gamma = adiabatic_exponent;
        }
    };

    using DeviceType = DeviceOp;

    explicit Channel(const std::string& prefix = "velocity_bc_channel");

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct Parabolic
{
    struct DeviceOp
    {
        int normal_dir{0};
        int tangential_dir{1};
        amrex::Real um{1.0};

        amrex::Real initial_density = 1.0;
        amrex::Real initial_temperature = 1.0 / 3.0;
        amrex::Real adiabatic_exponent = 5.0 / 3.0;

        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;
        amrex::Real Mach_ref{0.05};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            const int /*idir*/,
            const int /*ndir*/,
            const amrex::Real /*time*/,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& R,
            amrex::Real& T,
            amrex::Real& gamma) const
        {
            rho = initial_density;
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* prob_hi = geom.ProbHi();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real height =
                prob_hi[normal_dir] - prob_lo[normal_dir];
            const amrex::Real x =
                prob_lo[normal_dir] + (iv[normal_dir] + 0.5) * dx[normal_dir];

            vel[tangential_dir] =
                4.0 * um * x * (height - x) / (height * height);

            R = (R_u / m_bar);
            T = initial_temperature;
            gamma = adiabatic_exponent;
        }
    };

    using DeviceType = DeviceOp;

    explicit Parabolic(const std::string& prefix = "velocity_bc_parabolic");

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

} // namespace lbm::bc
#endif
