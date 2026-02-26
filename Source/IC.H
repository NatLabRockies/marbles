#ifndef IC_H
#define IC_H
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include "Stencil.H"
#include "Utilities.H"

namespace lbm::ic {

struct Constant
{
    struct DeviceOp
    {
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect mach_components = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initial_temperature = 1.0 / 3.0;
        amrex::Real adiabatic_exponent = 5.0 / 3.0;

        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;

        // Shape parameters
        int region_type{0}; // 0 = all, 1 = sphere, 2 = cylinder, 3 = box
        amrex::RealVect region_center = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::Real region_radius{1.0};
        amrex::Real cylinder_height{1.0};
        int cylinder_dir{2}; // 0 = x, 1 = y, 2 = z
        amrex::RealVect region_lo = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect region_hi = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::Real outside_density{0.0};

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5) * dx[1];
            const amrex::Real z =
                AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5) * dx[2]);

            bool is_inside = true;
            if (region_type == 1) { // sphere
                const amrex::Real r2 =
                    AMREX_D_TERM((x - region_center[0]) * (x - region_center[0]),
                                 + (y - region_center[1]) * (y - region_center[1]),
                                 + (z - region_center[2]) * (z - region_center[2]));
                if (r2 > region_radius * region_radius) {
                    is_inside = false;
                }
            } else if (region_type == 2) { // cylinder
                amrex::Real r2 = 0.0;
                amrex::Real h = 0.0;
                if (cylinder_dir == 0) {
                    r2 = AMREX_D_TERM(0.0,
                                      + (y - region_center[1]) * (y - region_center[1]),
                                      + (z - region_center[2]) * (z - region_center[2]));
                    h = std::abs(x - region_center[0]);
                } else if (cylinder_dir == 1) {
                    r2 = AMREX_D_TERM((x - region_center[0]) * (x - region_center[0]),
                                      + 0.0,
                                      + (z - region_center[2]) * (z - region_center[2]));
                    h = std::abs(y - region_center[1]);
                } else {
                    r2 = (x - region_center[0]) * (x - region_center[0]) +
                         (y - region_center[1]) * (y - region_center[1]);
                    h = AMREX_D_PICK(0.0, 0.0, std::abs(z - region_center[2]));
                }
                if (r2 > region_radius * region_radius ||
                    h > 0.5 * cylinder_height) {
                    is_inside = false;
                }
            } else if (region_type == 3) { // box
                if (AMREX_D_TERM(x < region_lo[0] || x > region_hi[0],
                                 || y < region_lo[1] || y > region_hi[1],
                                 || z < region_lo[2] || z > region_hi[2])) {
                    is_inside = false;
                }
            }

            if (is_inside) {
                rho = density;
            } else {
                rho = outside_density;
            }
            vel = velocity;
        }

        AMREX_GPU_DEVICE void thermal_initialize(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5) * dx[1];
            const amrex::Real z =
                AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5) * dx[2]);

            bool is_inside = true;
            if (region_type == 1) { // sphere
                const amrex::Real r2 =
                    AMREX_D_TERM((x - region_center[0]) * (x - region_center[0]),
                                 + (y - region_center[1]) * (y - region_center[1]),
                                 + (z - region_center[2]) * (z - region_center[2]));
                if (r2 > region_radius * region_radius) {
                    is_inside = false;
                }
            } else if (region_type == 2) { // cylinder
                amrex::Real r2 = 0.0;
                amrex::Real h = 0.0;
                if (cylinder_dir == 0) {
                    r2 = AMREX_D_TERM(0.0,
                                      + (y - region_center[1]) * (y - region_center[1]),
                                      + (z - region_center[2]) * (z - region_center[2]));
                    h = std::abs(x - region_center[0]);
                } else if (cylinder_dir == 1) {
                    r2 = AMREX_D_TERM((x - region_center[0]) * (x - region_center[0]),
                                      + 0.0,
                                      + (z - region_center[2]) * (z - region_center[2]));
                    h = std::abs(y - region_center[1]);
                } else {
                    r2 = (x - region_center[0]) * (x - region_center[0]) +
                         (y - region_center[1]) * (y - region_center[1]);
                    h = AMREX_D_PICK(0.0, 0.0, std::abs(z - region_center[2]));
                }
                if (r2 > region_radius * region_radius ||
                    h > 0.5 * cylinder_height) {
                    is_inside = false;
                }
            } else if (region_type == 3) { // box
                if (AMREX_D_TERM(x < region_lo[0] || x > region_hi[0],
                                 || y < region_lo[1] || y > region_hi[1],
                                 || z < region_lo[2] || z > region_hi[2])) {
                    is_inside = false;
                }
            }

            if (is_inside) {
                rho = density;
            } else {
                rho = outside_density;
            }
            vel = velocity;

            temperature = initial_temperature;
            R = R_u / m_bar;
            gamma = adiabatic_exponent;
        }
    };

    using DeviceType = DeviceOp;

    explicit Constant(const std::string& prefix = "ic_constant");

    static std::string identifier() { return "ic_constant"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct TaylorGreen
{
    struct DeviceOp
    {
        amrex::Real rho0{1.0};
        amrex::Real v0{1.0};
        amrex::RealVect omega = {AMREX_D_DECL(1.0, 1.0, 1.0)};
        amrex::Real L{1.0 / M_PI};

        amrex::Real initial_temperature = 1.0 / 3.0;
        const amrex::Real adiabatic_exponent = 5.0 / 3.0;

        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5) * dx[1];
            const amrex::Real z =
                AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5) * dx[2]);
            const amrex::Real omega0 = omega[0];
            const amrex::Real omega1 = omega[1];
            const amrex::Real omega2 = AMREX_D_PICK(, 0.0, omega[2]);

            rho = rho0 + rho0 * v0 * v0 / 16.0 *
                             (std::cos(2.0 * omega0 * x / L) +
                              std::cos(2.0 * omega1 * y / L)) *
                             (std::cos(2.0 * omega2 * z / L) + 2.0);
            AMREX_D_DECL(
                vel[0] = v0 * std::sin(omega0 * x / L) *
                         std::cos(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[1] = -v0 * std::cos(omega0 * x / L) *
                         std::sin(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[2] = 0.0);
        }

        AMREX_GPU_DEVICE
        void thermal_initialize(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();
            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5) * dx[1];
            const amrex::Real z =
                AMREX_D_PICK(, 0.0, prob_lo[2] + (iv[2] + 0.5) * dx[2]);
            const amrex::Real omega0 = omega[0];
            const amrex::Real omega1 = omega[1];
            const amrex::Real omega2 = AMREX_D_PICK(, 0.0, omega[2]);

            rho = rho0 + rho0 * v0 * v0 / 16.0 *
                             (std::cos(2.0 * omega0 * x / L) +
                              std::cos(2.0 * omega1 * y / L)) *
                             (std::cos(2.0 * omega2 * z / L) + 2.0);
            AMREX_D_DECL(
                vel[0] = v0 * std::sin(omega0 * x / L) *
                         std::cos(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[1] = -v0 * std::cos(omega0 * x / L) *
                         std::sin(omega1 * y / L) * std::cos(omega2 * z / L),
                vel[2] = 0.0);

            temperature = initial_temperature;
            R = R_u / m_bar;
            gamma = adiabatic_exponent;
        }
    };

    using DeviceType = DeviceOp;

    explicit TaylorGreen();

    static std::string identifier() { return "ic_taylorgreen"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct ViscosityTest
{
    struct DeviceOp
    {
        amrex::Real wave_length{1.0};
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect mach_components = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initial_temperature = 1.0 / 3.0;
        amrex::Real adiabatic_exponent = 5.0 / 3.0;

        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];

            rho = density;

            const amrex::Real a0 = 0.010 * speed_of_sound_ref;

            vel[0] = velocity[0] + a0 * std::sin(2.0 * M_PI * y / wave_length);
            vel[1] = velocity[1];
#if AMREX_SPACEDIM == 3
            vel[2] = velocity[2];
#endif
        }

        AMREX_GPU_DEVICE void thermal_initialize(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {

            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];

            rho = density;

            const amrex::Real a0 = 0.010 * speed_of_sound_ref;

            vel[0] = velocity[0] + a0 * std::sin(2.0 * M_PI * y / wave_length);
            vel[1] = velocity[1];
#if AMREX_SPACEDIM == 3
            vel[2] = velocity[2];
#endif

            temperature = initial_temperature;
            R = R_u / m_bar;
            gamma = adiabatic_exponent;
        }
    };

    using DeviceType = DeviceOp;

    explicit ViscosityTest();

    static std::string identifier() { return "ic_viscosity_test"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct ThermalDiffusivityTest
{
    struct DeviceOp
    {
        amrex::Real wave_length{1.0};
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect mach_components = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initial_temperature = 1.0 / 3.0;
        amrex::Real adiabatic_exponent = 5.0 / 3.0;
        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];

            const amrex::Real a0 = 0.0010 * initial_temperature;

            rho = density + a0 * std::sin(2.0 * M_PI * y / wave_length);

            vel[0] = velocity[0];
            vel[1] = velocity[1];
#if AMREX_SPACEDIM == 3
            vel[2] = velocity[2];
#endif
        }

        AMREX_GPU_DEVICE
        void thermal_initialize(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {

            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real y = prob_lo[1] + (iv[1] + 0.5 * 0.0) * dx[1];

            const amrex::Real a0 = 0.0010 * initial_temperature;

            R = R_u / m_bar;
            const amrex::Real pressure = density * R * initial_temperature;
            rho = density + a0 * std::sin(2.0 * M_PI * y / wave_length);

            vel[0] = velocity[0];
            vel[1] = velocity[1];
#if AMREX_SPACEDIM == 3
            vel[2] = velocity[2];
#endif

            gamma = adiabatic_exponent;

            temperature = pressure / (rho * R);
        }
    };

    using DeviceType = DeviceOp;

    explicit ThermalDiffusivityTest();

    static std::string identifier() { return "ic_thermaldiffusivity_test"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

struct SodTest
{
    struct DeviceOp
    {
        amrex::Real density_ratio{1.0};
        amrex::Real temperature_ratio{1.0};
        amrex::Real x_discontinuity{10.0};
        amrex::Real density{1.0};
        amrex::RealVect velocity = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        amrex::RealVect mach_components = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        amrex::Real initial_temperature = 1.0 / 3.0;

        amrex::Real adiabatic_exponent = 5.0 / 3.0;
        const amrex::Real R_u = constants::AIR_MOLAR_MASS;
        amrex::Real m_bar = constants::AIR_MOLAR_MASS;
        amrex::Real speed_of_sound_ref =
            std::sqrt(adiabatic_exponent * (R_u / m_bar) * initial_temperature);

        int model_type = 1;

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel) const
        {
            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];

            vel[0] = velocity[0];
            vel[1] = velocity[1];
#if AMREX_SPACEDIM == 3
            vel[2] = velocity[2];
#endif
            rho = density + 0.5 *
                                (1.0 + std::tanh((x - x_discontinuity) * 3.0)) *
                                (density_ratio * density - density);
        }

        AMREX_GPU_DEVICE
        void thermal_initialize(
            const amrex::IntVect& iv,
            amrex::GeometryData const& geom,
            amrex::Real& rho,
            amrex::RealVect& vel,
            amrex::Real& temperature,
            amrex::Real& R,
            amrex::Real& gamma) const
        {

            const amrex::Real* prob_lo = geom.ProbLo();
            const amrex::Real* dx = geom.CellSize();

            const amrex::Real x = prob_lo[0] + (iv[0] + 0.5 * 0.0) * dx[0];

            R = R_u / m_bar;

            vel[0] = velocity[0];
            vel[1] = velocity[1];
#if AMREX_SPACEDIM == 3
            vel[2] = velocity[2];
#endif

            gamma = adiabatic_exponent;

            rho = density + 0.5 *
                                (1.0 + std::tanh((x - x_discontinuity) * 3.0)) *
                                (density_ratio * density - density);
            temperature = initial_temperature +
                          0.5 * (1.0 + std::tanh((x - x_discontinuity) * 3.0)) *
                              (temperature_ratio * initial_temperature -
                               initial_temperature);
        }
    };

    using DeviceType = DeviceOp;

    explicit SodTest();

    static std::string identifier() { return "ic_sod"; }

    DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

class InitializerBase
{
public:
    InitializerBase() = default;

    virtual ~InitializerBase() = default;

    virtual void initialize(const int lev, amrex::GeometryData const& geom) = 0;
    virtual void initialize_lattice(const int lev, amrex::GeometryData const& geom, amrex::MultiFab& mf) = 0;
};

template <typename ICOp>
class Initializer : public InitializerBase
{
public:
    using ICDeviceOp = typename ICOp::DeviceType;

    Initializer(
        const amrex::Real& mesh_speed,
        const ICOp ic_op,
        amrex::Vector<amrex::MultiFab>& f)
        : m_mesh_speed(mesh_speed), m_op(ic_op), m_f(f), m_g(nullptr)
    {}

    Initializer(
        const amrex::Real& mesh_speed,
        const ICOp& ic_op,
        amrex::Vector<amrex::MultiFab>& f,
        amrex::Vector<amrex::MultiFab>& g)
        : m_mesh_speed(mesh_speed), m_op(ic_op), m_f(f), m_g(&g)
    {}

    ~Initializer() override = default;

    void initialize(const int lev, amrex::GeometryData const& geom) override
    {
        AMREX_ASSERT(m_g != nullptr);
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        const auto& weight = stencil.weights;
        auto const& f_arrs = m_f[lev].arrays();
        auto const& g_arrs = (*m_g)[lev].arrays();
        const amrex::Real l_mesh_speed = m_mesh_speed;
        ICDeviceOp ic(m_op.device_instance());

        amrex::RealVect zero_vec = {AMREX_D_DECL(0.0, 0.0, 0.0)};
        const amrex::Real theta0 = stencil::Stencil::THETA0;

        if (ic.model_type == 1) {
            amrex::ParallelFor(
                m_f[lev], m_f[lev].nGrowVect(), constants::N_MICRO_STATES,
                [=] AMREX_GPU_DEVICE(
                    int nbx, int i, int j, int k, int q) noexcept {
                    const amrex::Real wt = weight[q];
                    const auto& ev = evs[q];

                    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                    amrex::Real rho = 1.0;
                    amrex::RealVect vel(AMREX_D_DECL(0.0, 0.0, 0.0));
                    amrex::Real temperature = 1.0 / 3.0;
                    amrex::Real specific_gas_constant =
                        ic.R_u / constants::AIR_MOLAR_MASS;
                    amrex::Real gamma = 1.667;

                    ic.thermal_initialize(
                        iv, geom, rho, vel, temperature, specific_gas_constant,
                        gamma);

                    f_arrs[nbx](i, j, k, q) = set_equilibrium_value(
                        rho, vel, specific_gas_constant * temperature,
                        l_mesh_speed, wt, ev);

                    amrex::Real cv = specific_gas_constant / (gamma - 1.0);

                    amrex::Real two_rho_e =
                        get_energy(temperature, rho, vel, cv);

                    amrex::Real /*qxEq, qyEq, qzEq,*/ rxx_eq(0.0), ryy_eq(0.0),
                        rzz_eq(0.0), rxy_eq(0.0), rxz_eq(0.0), ryz_eq(0.0);

                    amrex::RealVect heat_flux = {AMREX_D_DECL(0.0, 0.0, 0.0)};

                    get_equilibrium_moments(
                        rho, vel, two_rho_e, cv, specific_gas_constant,
                        heat_flux, rxx_eq, ryy_eq, rzz_eq, rxy_eq, rxz_eq,
                        ryz_eq);

                    amrex::GpuArray<amrex::Real, 6> flux_of_heat_flux = {
                        rxx_eq, ryy_eq, rzz_eq, rxy_eq, rxz_eq, ryz_eq};

                    g_arrs[nbx](i, j, k, q) =
                        set_extended_grad_expansion_generic(
                            two_rho_e, heat_flux, flux_of_heat_flux,
                            l_mesh_speed, wt, ev, theta0, zero_vec, 1.0);
                });
        } else {
            amrex::Abort("Invalid model_type");
        }

        amrex::Gpu::synchronize();
    }

    void initialize_lattice(const int lev, amrex::GeometryData const& geom, amrex::MultiFab& mf) override
    {
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        const auto& weight = stencil.weights;
        auto const& f_arrs = mf.arrays();
        const amrex::Real l_mesh_speed = m_mesh_speed;
        ICDeviceOp ic(m_op.device_instance());

        if (ic.model_type == 1) {
            amrex::ParallelFor(
                mf, mf.nGrowVect(), constants::N_MICRO_STATES,
                [=] AMREX_GPU_DEVICE(
                    int nbx, int i, int j, int k, int q) noexcept {
                    const amrex::Real wt = weight[q];
                    const auto& ev = evs[q];

                    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                    amrex::Real rho = 1.0;
                    amrex::RealVect vel(AMREX_D_DECL(0.0, 0.0, 0.0));
                    amrex::Real temperature = 1.0 / 3.0;
                    amrex::Real specific_gas_constant =
                        ic.R_u / constants::AIR_MOLAR_MASS;
                    amrex::Real gamma = 1.667;

                    ic.thermal_initialize(
                        iv, geom, rho, vel, temperature, specific_gas_constant,
                        gamma);

                    f_arrs[nbx](i, j, k, q) = set_equilibrium_value(
                        rho, vel, specific_gas_constant * temperature,
                        l_mesh_speed, wt, ev);
                });
        } else {
            amrex::Abort("Invalid model_type");
        }

        amrex::Gpu::synchronize();
    }

    const amrex::Real m_mesh_speed;
    const ICOp m_op;
    amrex::Vector<amrex::MultiFab>& m_f;
    amrex::Vector<amrex::MultiFab>* m_g;
};

} // namespace lbm::ic
#endif
