#include "IC.H"
namespace lbm::ic {

Constant::Constant(const std::string& prefix)
{
    amrex::ParmParse pp(prefix);
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    amrex::ParmParse ppl("lbm");
    ppl.query("initial_temperature", m_op.initial_temperature);
    ppl.query("adiabatic_exponent", m_op.adiabatic_exponent);
    ppl.query("mean_molecular_mass", m_op.m_bar);

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
    }

    std::string region_type_str = "all";
    pp.query("region_type", region_type_str);
    if (region_type_str == "all") {
        m_op.region_type = 0;
    } else if (region_type_str == "sphere") {
        m_op.region_type = 1;
    } else if (region_type_str == "cylinder") {
        m_op.region_type = 2;
    } else if (region_type_str == "box") {
        m_op.region_type = 3;
    } else {
        amrex::Abort(
            "Constant::Constant(): Unknown region_type: " + region_type_str);
    }

    amrex::Vector<amrex::Real> region_center{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("region_center", region_center, 0, AMREX_SPACEDIM);
    for (int n = 0; n < region_center.size(); n++) {
        m_op.region_center[n] = region_center[n];
    }

    pp.query("region_radius", m_op.region_radius);

    amrex::Vector<amrex::Real> region_lo_v{AMREX_D_DECL(0.0, 0.0, 0.0)};
    amrex::Vector<amrex::Real> region_hi_v{AMREX_D_DECL(0.0, 0.0, 0.0)};
    pp.queryarr("region_lo", region_lo_v, 0, AMREX_SPACEDIM);
    pp.queryarr("region_hi", region_hi_v, 0, AMREX_SPACEDIM);
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        m_op.region_lo[n] = region_lo_v[n];
        m_op.region_hi[n] = region_hi_v[n];
    }

    pp.query("cylinder_height", m_op.cylinder_height);
    std::string cylinder_dir_str = "z";
    pp.query("cylinder_direction", cylinder_dir_str);
    if (cylinder_dir_str == "x") {
        m_op.cylinder_dir = 0;
    } else if (cylinder_dir_str == "y") {
        m_op.cylinder_dir = 1;
    } else if (cylinder_dir_str == "z") {
        m_op.cylinder_dir = 2;
    } else {
        amrex::Abort(
            "Constant::Constant(): Unknown cylinder_direction: " +
            cylinder_dir_str);
    }

    pp.query("outside_density", m_op.outside_density);
}

TaylorGreen::TaylorGreen()
{
    amrex::ParmParse pp(identifier());
    pp.query("rho0", m_op.rho0);
    pp.query("v0", m_op.v0);

    amrex::Vector<amrex::Real> omega{AMREX_D_DECL(1.0, 1.0, 1.0)};
    pp.queryarr("omega", omega, 0, AMREX_SPACEDIM);
    for (int n = 0; n < omega.size(); n++) {
        m_op.omega[n] = omega[n];
    }
}

ViscosityTest::ViscosityTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("wave_length", m_op.wave_length);
    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    amrex::ParmParse ppl("lbm");
    ppl.query("initial_temperature", m_op.initial_temperature);
    ppl.query("adiabatic_exponent", m_op.adiabatic_exponent);
    ppl.query("mean_molecular_mass", m_op.m_bar);

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);

    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
    }
}

ThermalDiffusivityTest::ThermalDiffusivityTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("wave_length", m_op.wave_length);
    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    amrex::ParmParse ppl("lbm");
    ppl.query("initial_temperature", m_op.initial_temperature);
    ppl.query("adiabatic_exponent", m_op.adiabatic_exponent);
    ppl.query("mean_molecular_mass", m_op.m_bar);

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
    }
}

SodTest::SodTest()
{
    amrex::ParmParse pp(identifier());
    pp.query("density", m_op.density);

    amrex::Vector<amrex::Real> velocity{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
    for (int n = 0; n < velocity.size(); n++) {
        m_op.velocity[n] = velocity[n];
    }

    pp.query("density_ratio", m_op.density_ratio);
    pp.query("temperature_ratio", m_op.temperature_ratio);
    pp.query("x_discontinuity", m_op.x_discontinuity);

    m_op.model_type = 1;
    amrex::Vector<amrex::Real> mach_components{AMREX_D_DECL(0, 0, 0)};
    pp.queryarr("mach_components", mach_components, 0, AMREX_SPACEDIM);
    for (int n = 0; n < mach_components.size(); n++) {
        m_op.mach_components[n] = mach_components[n];
    }

    amrex::ParmParse ppl("lbm");
    ppl.query("initial_temperature", m_op.initial_temperature);
    ppl.query("adiabatic_exponent", m_op.adiabatic_exponent);
    ppl.query("mean_molecular_mass", m_op.m_bar);

    pp.query("initial_temperature", m_op.initial_temperature);
    pp.query("adiabatic_exponent", m_op.adiabatic_exponent);
    pp.query("mean_molecular_mass", m_op.m_bar);
    m_op.speed_of_sound_ref = std::sqrt(
        m_op.adiabatic_exponent * (m_op.R_u / m_op.m_bar) *
        m_op.initial_temperature);

    for (int n = 0; n < mach_components.size(); n++) {
        m_op.velocity[n] = m_op.mach_components[n] * m_op.speed_of_sound_ref;
    }
}

} // namespace lbm::ic
