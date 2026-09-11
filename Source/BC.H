#ifndef BC_H
#define BC_H
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include "Constants.H"
#include "Utilities.H"
#include "Stencil.H"
#include "VelocityBC.H"

namespace lbm::bc {

static constexpr int PERIODIC = 0;
static constexpr int NOSLIPWALL = 1;
static constexpr int VELOCITY = 2;
static constexpr int PRESSURE = 3;
static constexpr int OUTFLOW_ZEROTH_ORDER = 5;
static constexpr int SLIPWALLXNORMAL = 6;
static constexpr int SLIPWALLYNORMAL = 7;
static constexpr int SLIPWALLZNORMAL = 8;

template <typename VelBCOp>
struct BCFill
{

    using VelBCOpType = VelBCOp;

    const amrex::GpuArray<amrex::IntVect, constants::N_MICRO_STATES> m_evs =
        stencil::Stencil().evs;
    const amrex::GpuArray<int, constants::N_MICRO_STATES> m_bounce_dirs =
        stencil::Stencil().bounce_dirs;
    const amrex::GpuArray<int, constants::N_MICRO_STATES> m_bounce_Xdir =
        stencil::Stencil().bounce_Xdir;
    const amrex::GpuArray<int, constants::N_MICRO_STATES> m_bounce_Ydir =
        stencil::Stencil().bounce_Ydir;
    const amrex::GpuArray<int, constants::N_MICRO_STATES> m_bounce_Zdir =
        stencil::Stencil().bounce_Zdir;
    const amrex::GpuArray<amrex::Real, constants::N_MICRO_STATES> m_weights =
        stencil::Stencil().weights;

    const amrex::Real m_theta0 = stencil::Stencil::THETA0;
    amrex::RealVect m_zero_vec = {AMREX_D_DECL(0.0, 0.0, 0.0)};

    const amrex::Real m_mesh_speed;
    const amrex::GpuArray<int, AMREX_SPACEDIM * 2> m_bc_type = {0};
    const amrex::IntVect m_ng = {AMREX_D_DECL(0, 0, 0)};
    const VelBCOpType m_vel_bc_op;
    bool m_is_energy_lattice;

    AMREX_GPU_HOST
    constexpr explicit BCFill(
        const amrex::Real& mesh_speed,

        const amrex::GpuArray<int, AMREX_SPACEDIM * 2>& bc_type,
        const amrex::IntVect& ng,
        const VelBCOpType& vel_bc_op,
        const bool& m_is_energy_lattice_ = false)
        : m_mesh_speed(mesh_speed)
        , m_bc_type(bc_type)
        , m_ng(ng)
        , m_vel_bc_op(vel_bc_op)
        , m_is_energy_lattice(m_is_energy_lattice_)
    {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bounce_back_bc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to bounce back. 2D e.g., for
        // the ghost cell at at (-1, -1), store the
        // value f_(0, 0) associated with e=[-1, -1] and
        // orient it such that it points e=[1, 1]. This
        // reproduces the bounce back effect directly
        // when the streaming step is performed.

        data(iv, q) = data(iv + ev, m_bounce_dirs[q]);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bounce_back_xbc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to bounce back. 2D e.g., for
        // the ghost cell at at (-1, -1), store the
        // value f_(0, 0) associated with e=[-1, -1] and
        // orient it such that it points e=[1, 1]. This
        // reproduces the bounce back effect directly
        // when the streaming step is performed.

        data(iv, q) = data(iv + ev, m_bounce_Xdir[q]);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bounce_back_ybc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to bounce back. 2D e.g., for
        // the ghost cell at at (-1, -1), store the
        // value f_(0, 0) associated with e=[-1, -1] and
        // orient it such that it points e=[1, 1]. This
        // reproduces the bounce back effect directly
        // when the streaming step is performed.

        data(iv, q) = data(iv + ev, m_bounce_Ydir[q]);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bounce_back_zbc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to bounce back. 2D e.g., for
        // the ghost cell at at (-1, -1), store the
        // value f_(0, 0) associated with e=[-1, -1] and
        // orient it such that it points e=[1, 1]. This
        // reproduces the bounce back effect directly
        // when the streaming step is performed.

        data(iv, q) = data(iv + ev, m_bounce_Zdir[q]);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void velocity_bc_g(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        const int /*idir*/,
        const int /*ndir*/,
        const amrex::Real& rho_bc,
        const amrex::RealVect& vel_bc,
        const amrex::Real& R,
        const amrex::Real& temperature,
        const amrex::Real& adiabatic_exponent,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to a prescribed velocity. This
        // is done by assigning all the values of f_ on
        // this ghost layer to be the equilibrium
        // distributions associated with the desired
        // velocity. They will naturally exert their
        // effect on the interior fluid points during
        // the streaming step.

        // In this case, a value for u, v, and w must be
        // specified, either hard coded or set by a
        // user-specified function

        const amrex::Real cv = R / (adiabatic_exponent - 1.0);

        const amrex::Real two_rho_e =
            get_energy(temperature, rho_bc, vel_bc, cv);

        amrex::Real rxx_eq(0.0), ryy_eq(0.0), rzz_eq(0.0), rxy_eq(0.0),
            rxz_eq(0.0), ryz_eq(0.0);

        amrex::RealVect heat_flux = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        get_equilibrium_moments(
            rho_bc, vel_bc, two_rho_e, cv, R, heat_flux, rxx_eq, ryy_eq, rzz_eq,
            rxy_eq, rxz_eq, ryz_eq);

        amrex::GpuArray<amrex::Real, 6> flux_of_heat_flux = {
            rxx_eq, ryy_eq, rzz_eq, rxy_eq, rxz_eq, ryz_eq};

        data(iv, q) = set_extended_grad_expansion_generic(
            two_rho_e, heat_flux, flux_of_heat_flux, m_mesh_speed, m_weights[q],
            ev, m_theta0, m_zero_vec, 1.0);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void velocity_bc_f(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        const int /*idir*/,
        const int /*ndir*/,
        const amrex::Real& rho_bc,
        const amrex::RealVect& vel_bc,
        const amrex::Real& RT,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to a prescribed velocity. This
        // is done by assigning all the values of f_ on
        // this ghost layer to be the equilibrium
        // distributions associated with the desired
        // velocity. They will naturally exert their
        // effect on the interior fluid points during
        // the streaming step.

        // In this case, a value for u, v, and w must be
        // specified, either hard coded or set by a
        // user-specified function
        data(iv, q) = set_equilibrium_value(
            rho_bc, vel_bc, RT, m_mesh_speed, m_weights[q], ev);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void pressure_bc_f(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        const int idir,
        const int ndir,
        const amrex::Real rho_bc,
        const amrex::Real& RT,
        amrex::Array4<amrex::Real> const& data,
        amrex::Real& rho_bc_out) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to a prescribed
        // pressure/density. This is done by assigning
        // all the values of f_ on this ghost layer to
        // be the equilibrium distributions associated
        // with the desired density. They will naturally
        // exert their effect on the interior fluid
        // points during the streaming step.

        amrex::Real rho_out = 0.0;
        amrex::Real rho_tan = 0.0;

        const amrex::IntVect ivn(iv + ev);

        // Build the 3-part rho calculation at the
        // (iv+ev) fluid cell
        for (int qq = 0; qq < constants::N_MICRO_STATES; qq++) {

            const auto& evn_out = m_evs[qq];
            const int bounce_qq = m_bounce_dirs[qq];
            const auto& evn_in = m_evs[bounce_qq];

            if (evn_in[idir] == -ndir) {
                rho_out += 2.0 * data(ivn + evn_out, bounce_qq);
            } else if (evn_in[idir] == 0) {
                rho_tan += 1.0 * data(ivn + evn_out, bounce_qq);
            }
        }

        amrex::RealVect vel_bc(AMREX_D_DECL(0.0, 0.0, 0.0));
        vel_bc[idir] = ndir * (1.0 - (rho_out + rho_tan) / rho_bc);

        data(iv, q) = set_equilibrium_value(
            rho_bc, vel_bc, RT, m_mesh_speed, m_weights[q], ev);

        rho_bc_out = rho_out + rho_tan;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void pressure_bc_g(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& ev,
        const int idir,
        const int ndir,
        const amrex::Real rho_bc_out,
        const amrex::Real rho_bc,
        const amrex::Real& R,
        const amrex::Real& temperature,
        const amrex::Real& adiabatic_exponent,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to a prescribed
        // pressure/density. This is done by assigning
        // all the values of f_ on this ghost layer to
        // be the equilibrium distributions associated
        // with the desired density. They will naturally
        // exert their effect on the interior fluid
        // points during the streaming step.

        // Build the 3-part rho calculation at the
        // (iv+ev) fluid cell

        amrex::RealVect vel_bc(AMREX_D_DECL(0.0, 0.0, 0.0));
        vel_bc[idir] = ndir * (1.0 - (rho_bc_out) / rho_bc);

        const amrex::Real cv = R / (adiabatic_exponent - 1.0);

        const amrex::Real two_rho_e =
            get_energy(temperature, rho_bc, vel_bc, cv);

        amrex::Real rxx_eq(0.0), ryy_eq(0.0), rzz_eq(0.0), rxy_eq(0.0),
            rxz_eq(0.0), ryz_eq(0.0);

        amrex::RealVect heat_flux = {AMREX_D_DECL(0.0, 0.0, 0.0)};

        get_equilibrium_moments(
            rho_bc, vel_bc, two_rho_e, cv, R, heat_flux, rxx_eq, ryy_eq, rzz_eq,
            rxy_eq, rxz_eq, ryz_eq);

        amrex::GpuArray<amrex::Real, 6> flux_of_heat_flux = {
            rxx_eq, ryy_eq, rzz_eq, rxy_eq, rxz_eq, ryz_eq};

        data(iv, q) = set_extended_grad_expansion_generic(
            two_rho_e, heat_flux, flux_of_heat_flux, m_mesh_speed, m_weights[q],
            ev, m_theta0, m_zero_vec, 1.0);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void simple_outflow_bc(
        const amrex::IntVect& iv,
        const int q,
        const amrex::IntVect& /*ev*/,
        const int idir,
        const int ndir,
        amrex::Array4<amrex::Real> const& data) const
    {
        // Prepare the ghost cell layer of f_ for
        // streaming due to an outflow condition This
        // copies the rate of change from the last two
        // interior cell onto the ghost layer such that
        // when streaming you have the effect of "no
        // change" in the direction across the boundary.
        amrex::IntVect iv_shift(amrex::IntVect::TheDimensionVector(idir));

        data(iv, q) = data(iv + 1 * ndir * iv_shift, q);
    }

    AMREX_GPU_DEVICE
    void operator()(
        const amrex::IntVect& iv,
        amrex::Array4<amrex::Real> const& data,
        const int /*dcomp*/,
        const int /*numcomp*/,
        amrex::GeometryData const& geom,
        const amrex::Real time,
        const amrex::BCRec* /*bcr*/,
        const int /*bcomp*/,
        const int /*orig_comp*/) const
    {
        // do something for external Dirichlet (BCType::ext_dir)
        const int* domlo = geom.Domain().loVect();
        const int* domhi = geom.Domain().hiVect();

        // Define an inward-pointing normal vector
        // that is, positive on the x-, y-, z- wall and
        // negative on the x+, y+, z+ wall
        const amrex::GpuArray<int, 2> normal_dirs = {1, -1};

        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
            for (int lohi = 0; lohi < 2; lohi++) {

                const int normal_dir = normal_dirs[lohi];
                const int bc_idx = idir + lohi * AMREX_SPACEDIM;

                if ((((lohi == 0) && (iv[idir] < domlo[idir])) ||
                     ((lohi == 1) && (iv[idir] > domhi[idir])))) {

                    const auto& lb = amrex::lbound(data);
                    const auto& ub = amrex::ubound(data);
                    const amrex::IntVect liv(
                        amrex::IntVect(AMREX_D_DECL(lb.x, lb.y, lb.z)) + m_ng);
                    const amrex::IntVect uiv(
                        amrex::IntVect(AMREX_D_DECL(ub.x, ub.y, ub.z)) - m_ng);
                    const amrex::Box inside(liv, uiv);

                    for (int q = 0; q < constants::N_MICRO_STATES; q++) {

                        const auto& ev = m_evs[q];
                        const amrex::IntVect ivn(iv + ev);

                        // If this distribution points into (enters) the
                        // domain, and we are accessing a cell within the data
                        // bounds
                        if (inside.contains(ivn)) {

                            if (m_bc_type[bc_idx] == NOSLIPWALL) {
                                bounce_back_bc(iv, q, ev, data);
                            } else if (m_bc_type[bc_idx] == VELOCITY) {
                                amrex::RealVect vel_bc(
                                    AMREX_D_DECL(0.0, 0.0, 0.0));
                                amrex::Real specific_gas_constant(1.0),
                                    temperature(1.0 / 3.0);
                                amrex::Real gamma(5.0 / 3.0);

                                amrex::Real rho_bc(0.0);

                                m_vel_bc_op(
                                    iv, idir, normal_dir, time, geom, rho_bc,
                                    vel_bc, specific_gas_constant, temperature,
                                    gamma);

                                if (m_is_energy_lattice) {
                                    velocity_bc_g(
                                        iv, q, ev, idir, normal_dir, rho_bc,
                                        vel_bc, specific_gas_constant,
                                        temperature, gamma, data);
                                } else {
                                    velocity_bc_f(
                                        iv, q, ev, idir, normal_dir, rho_bc,
                                        vel_bc,
                                        specific_gas_constant * temperature,
                                        data);
                                }

                            } else if (m_bc_type[bc_idx] == PRESSURE) {

                                amrex::RealVect vel_bc(
                                    AMREX_D_DECL(0.0, 0.0, 0.0));
                                amrex::Real specific_gas_constant(1.0),
                                    temperature(1.0 / 3.0);
                                amrex::Real gamma(5.0 / 3.0);

                                amrex::Real rho_bcout(1.0);
                                amrex::Real rho_bc(1.0);

                                m_vel_bc_op(
                                    iv, idir, normal_dir, time, geom, rho_bc,
                                    vel_bc, specific_gas_constant, temperature,
                                    gamma);

                                if (m_is_energy_lattice) {

                                    rho_bcout = 1.0;
                                    pressure_bc_g(
                                        iv, q, ev, idir, normal_dir, rho_bcout,
                                        rho_bc, specific_gas_constant,
                                        temperature, gamma, data);

                                } else {
                                    pressure_bc_f(
                                        iv, q, ev, idir, normal_dir, rho_bc,
                                        specific_gas_constant * temperature,
                                        data, rho_bcout);
                                }

                            } else if (
                                m_bc_type[bc_idx] == OUTFLOW_ZEROTH_ORDER) {
                                simple_outflow_bc(
                                    iv, q, ev, idir, normal_dir, data);
                            } else if (m_bc_type[bc_idx] == SLIPWALLXNORMAL) {
                                bounce_back_xbc(iv, q, ev, data);
                            } else if (m_bc_type[bc_idx] == SLIPWALLYNORMAL) {
                                bounce_back_ybc(iv, q, ev, data);
                            } else if (m_bc_type[bc_idx] == SLIPWALLZNORMAL) {
                                bounce_back_zbc(iv, q, ev, data);
                            }
                        } else {
                            // to avoid spurious NaNs fill with a dummy value
                            // Use 0.0 instead of -1.0 to avoid negative
                            // densities in component lattices (passive scalars)
                            data(iv, q) = 0.0;
                        }
                    }
                }
            }
        }
    }
};

template <typename VelBCOp>
struct BCOpCreator
{
    using VelBCOpType = typename VelBCOp::DeviceType;
    using FunctorType = BCFill<VelBCOpType>;

    explicit BCOpCreator(
        const amrex::Real& mesh_speed,
        const amrex::GpuArray<int, AMREX_SPACEDIM * 2>& bc_type,
        const amrex::IntVect& ng,
        bool m_is_energy_lattice_ = false)
        : m_mesh_speed(mesh_speed)
        , m_bc_type(bc_type)
        , m_ng(ng)
        , m_vel_bc_op()
        , m_is_energy_lattice(m_is_energy_lattice_)
    {}

    explicit BCOpCreator(
        const amrex::Real& mesh_speed,
        const amrex::GpuArray<int, AMREX_SPACEDIM * 2>& bc_type,
        const amrex::IntVect& ng,
        const std::string& prefix,
        bool m_is_energy_lattice_ = false)
        : m_mesh_speed(mesh_speed)
        , m_bc_type(bc_type)
        , m_ng(ng)
        , m_vel_bc_op(prefix)
        , m_is_energy_lattice(m_is_energy_lattice_)
    {}

    FunctorType operator()() const
    {
        return FunctorType{
            m_mesh_speed, m_bc_type, m_ng, m_vel_bc_op.device_instance(),
            m_is_energy_lattice};
    }

    const amrex::Real m_mesh_speed;
    const amrex::GpuArray<int, AMREX_SPACEDIM * 2> m_bc_type = {0};
    const amrex::IntVect m_ng = {AMREX_D_DECL(0, 0, 0)};
    VelBCOp m_vel_bc_op;
    bool m_is_energy_lattice;
};

} // namespace lbm::bc
#endif
