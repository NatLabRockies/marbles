#include "VelocityBC.H"

namespace lbm::bc {
NoOp::NoOp(const std::string& /*prefix*/) {}

Constant::Constant(const std::string& prefix)
{
    amrex::ParmParse pp(prefix);
    pp.query("dir", m_op.dir);

    m_op.model_type = 1;
    pp.query("Mach_ref", m_op.Mach_ref);
    pp.query("initial_density", m_op.initial_density);
    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);
    m_op.u0 = m_op.Mach_ref * m_op.speed_of_sound_ref;
}

Channel::Channel(const std::string& prefix)
{
    amrex::ParmParse pp(prefix);

    m_op.model_type = 1;
    pp.query("Mach_ref", m_op.Mach_ref);
    pp.query("initial_density", m_op.initial_density);
    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);
    m_op.u_ref = m_op.Mach_ref * m_op.speed_of_sound_ref;
}

Parabolic::Parabolic(const std::string& prefix)
{
    amrex::ParmParse pp(prefix);
    pp.query("normal_dir", m_op.normal_dir);
    pp.query("tangential_dir", m_op.tangential_dir);

    m_op.model_type = 1;
    pp.query("Mach_ref", m_op.Mach_ref);
    pp.query("initial_density", m_op.initial_density);
    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);
    m_op.um = m_op.Mach_ref * m_op.speed_of_sound_ref;
}

} // namespace lbm::bc
