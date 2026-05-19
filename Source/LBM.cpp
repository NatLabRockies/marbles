#include <fstream>
#include <iomanip>
#include <memory>
#include <AMReX_Parser.H>

#include "LBM.H"

namespace lbm {
LBM::LBM()
{
    BL_PROFILE("LBM::LBM()");
    read_parameters();

    int nlevs_max = max_level + 1;
    initialize_eb(Geom(maxLevel()), maxLevel());

    m_macrodata_varnames.push_back("rho");
    m_macrodata_varnames.push_back("vel_x");
    m_macrodata_varnames.push_back("vel_y");
    m_macrodata_varnames.push_back("vel_z");
    m_macrodata_varnames.push_back("vel_mag");

    m_macrodata_varnames.push_back("two_rho_e");
    m_macrodata_varnames.push_back("QCorrX");
    m_macrodata_varnames.push_back("QCorrY");
    m_macrodata_varnames.push_back("QCorrZ");
    m_macrodata_varnames.push_back("pxx");
    m_macrodata_varnames.push_back("pyy");
    m_macrodata_varnames.push_back("pzz");
    m_macrodata_varnames.push_back("pxy");
    m_macrodata_varnames.push_back("pxz");
    m_macrodata_varnames.push_back("pyz");
    m_macrodata_varnames.push_back("qx");
    m_macrodata_varnames.push_back("qy");
    m_macrodata_varnames.push_back("qz");
    m_macrodata_varnames.push_back("temperature");

    const size_t n_zero = 2;
    for (int q = 0; q < constants::N_MICRO_STATES; q++) {
        const auto num_str = std::to_string(q);
        const auto zero_padded_str =
            std::string(n_zero - std::min(n_zero, num_str.length()), '0') +
            num_str;
        m_microdata_varnames.push_back("f_" + zero_padded_str);
    }

    for (int q = 0; q < constants::N_MICRO_STATES; q++) {
        const auto num_str = std::to_string(q);
        const auto zero_padded_str =
            std::string(n_zero - std::min(n_zero, num_str.length()), '0') +
            num_str;
        m_microdata_g_varnames.push_back("g_" + zero_padded_str);
    }

    m_deriveddata_varnames.push_back("vort_x");
    m_deriveddata_varnames.push_back("vort_y");
    m_deriveddata_varnames.push_back("vort_z");
    m_deriveddata_varnames.push_back("vort_mag");

    m_deriveddata_varnames.push_back("dQCorrX");
    m_deriveddata_varnames.push_back("dQCorrY");
    m_deriveddata_varnames.push_back("dQCorrZ");

    // Energy dissipation rate (added for kLa two-phase model)
    m_deriveddata_varnames.push_back("epsilon");

    m_idata_varnames.push_back("is_fluid");
    m_idata_varnames.push_back("eb_boundary");
    m_idata_varnames.push_back("eb_fluid_boundary");
    // placeholder name for the new 4th component of m_is_fluid
    m_idata_varnames.push_back("eb_fluid_boundary_2");
    // fractional field is stored separately as a Real MultiFab
    m_fracdata_varnames.push_back("is_fluid_fraction");
    for (const auto& vname : m_macrodata_varnames) {
        m_lbm_varnames.push_back(vname);
    }
    if (m_save_streaming) {
        for (const auto& vname : m_microdata_varnames) {
            m_lbm_varnames.push_back(vname);
        }

        for (const auto& vname : m_microdata_g_varnames) {
            m_lbm_varnames.push_back(vname);
        }
    }
    if (m_save_derived) {
        for (const auto& vname : m_deriveddata_varnames) {
            m_lbm_varnames.push_back(vname);
        }
    }
    for (const auto& vname : m_idata_varnames) {
        m_lbm_varnames.push_back(vname);
    }
    for (const auto& vname : m_fracdata_varnames) {
        m_lbm_varnames.push_back(vname);
    }

    for (int i = 0; i < m_n_components; ++i) {
        m_lbm_varnames.push_back("Y_" + std::to_string(i));
    }

    read_tagging_parameters();

    m_isteps.resize(nlevs_max, 0);
    m_nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        m_nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    m_ts_new.resize(nlevs_max, 0.0);
    m_ts_old.resize(nlevs_max, constants::LOW_NUM);
    m_dts.resize(nlevs_max, constants::LARGE_NUM);

    m_macrodata.resize(nlevs_max);
    m_f.resize(nlevs_max);
    m_component_lattices.resize(m_n_components);
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i].resize(nlevs_max);
    }
    m_g.resize(nlevs_max);
    m_eq.resize(nlevs_max);
    m_eq_g.resize(nlevs_max);
    m_derived.resize(nlevs_max);
    m_is_fluid.resize(nlevs_max);
    m_is_fluid_fraction.resize(nlevs_max);
    m_plt_mf.resize(nlevs_max);
    m_mask.resize(nlevs_max);
    m_stationary_mask.resize(nlevs_max);
    m_cell_type.resize(nlevs_max);
    m_phi_fslbm.resize(nlevs_max);

    m_factory.resize(nlevs_max);
    // BCs
    m_bcs.resize(constants::N_MICRO_STATES);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // lo-side BCs
        if (m_bc_lo[idim] == bc::PERIODIC) {
            for (auto& bc : m_bcs) {
                bc.setLo(idim, amrex::BCType::int_dir);
            }
        } else if (
            (m_bc_lo[idim] == bc::NOSLIPWALL) ||
            (m_bc_lo[idim] == bc::SLIPWALLXNORMAL) ||
            (m_bc_lo[idim] == bc::SLIPWALLYNORMAL) ||
            (m_bc_lo[idim] == bc::SLIPWALLZNORMAL) ||
            (m_bc_lo[idim] == bc::VELOCITY) ||
            (m_bc_lo[idim] == bc::PRESSURE) ||
            (m_bc_lo[idim] == bc::OUTFLOW_ZEROTH_ORDER)) {
            for (auto& bc : m_bcs) {
                bc.setLo(idim, amrex::BCType::ext_dir);
            }
        } else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCs
        if (m_bc_hi[idim] == bc::PERIODIC) {
            for (auto& bc : m_bcs) {
                bc.setHi(idim, amrex::BCType::int_dir);
            }
        } else if (
            (m_bc_hi[idim] == bc::NOSLIPWALL) ||
            (m_bc_hi[idim] == bc::SLIPWALLXNORMAL) ||
            (m_bc_hi[idim] == bc::SLIPWALLYNORMAL) ||
            (m_bc_hi[idim] == bc::SLIPWALLZNORMAL) ||
            (m_bc_hi[idim] == bc::VELOCITY) ||
            (m_bc_hi[idim] == bc::PRESSURE) ||
            (m_bc_hi[idim] == bc::OUTFLOW_ZEROTH_ORDER)) {
            for (auto& bc : m_bcs) {
                bc.setHi(idim, amrex::BCType::ext_dir);
            }
        } else {
            amrex::Abort("Invalid bc_hi");
        }
    }
}

LBM::~LBM() = default;

void LBM::init_data()
{
    BL_PROFILE("LBM::init_data()");

    stencil::check_stencil();

    if (m_restart_chkfile.empty()) {
        // start simulation from the beginning
        const amrex::Real time = 0.0;
        set_ics();
        InitFromScratch(time);

        average_down(amrex::IntVect(0));

        compute_dt();

        // Initialize Lagrangian bubble container (after grids are finalized)
        if (m_enable_bubbles) {
            m_bubbles.initialize(Geom(0), grids[0], dmap[0], m_bubble_params);
        }

        if (m_chk_int > 0) {
            write_checkpoint_file();
        }

        open_forces_file(true);
        open_species_stats_file(true);
        compute_eb_forces();
    } else {
        // restart from a checkpoint
        read_checkpoint_file();

        open_forces_file(false);
        open_species_stats_file(false);
    }

    if (m_plot_int > 0) {
        write_plot_file();
    }

    set_bcs();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Grid summary: " << std::endl;
        printGridSummary(amrex::OutStream(), 0, finest_level);
    }
}

void LBM::read_parameters()
{
    BL_PROFILE("LBM::read_parameters()");

    {
        amrex::ParmParse pp;
        pp.query("max_step", m_max_step);
        pp.query("stop_time", m_stop_time);
    }

    {
        amrex::ParmParse pp("amr");

        pp.query("regrid_int", m_regrid_int);
        pp.query("plot_file", m_plot_file);
        pp.query("plot_int", m_plot_int);
        m_print_int = (m_plot_int > 0) ? amrex::max(1, m_plot_int / 10) : m_print_int;  // default to plot_int/10
        pp.query("print_int", m_print_int);
        pp.query("chk_file", m_chk_file);
        pp.query("chk_int", m_chk_int);
        pp.query("restart", m_restart_chkfile);
        pp.query("file_name_digits", m_file_name_digits);
    }

    {
        amrex::ParmParse pp("lbm");
        pp.query("n_components", m_n_components);
        pp.queryarr("bc_lo", m_bc_lo, 0, AMREX_SPACEDIM);
        pp.queryarr("bc_hi", m_bc_hi, 0, AMREX_SPACEDIM);
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            m_bc_type[i] = m_bc_lo[i];
            m_bc_type[i + AMREX_SPACEDIM] = m_bc_hi[i];
        }

        // Check bcs against possible periodic geometry
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            // if it's periodic, it must have internal BC marked.
            if (amrex::DefaultGeometry().isPeriodic(dir)) {
                if (m_bc_lo[dir] != bc::PERIODIC) {
                    amrex::Abort(
                        "BC is periodic in direction " + std::to_string(dir) +
                        " but low BC is not 0");
                }
                if (m_bc_hi[dir] != bc::PERIODIC) {
                    amrex::Abort(
                        "BC is periodic in direction " + std::to_string(dir) +
                        " but high BC is not 0");
                }
            } else {
                // If not periodic, should not be interior.
                if (m_bc_lo[dir] == bc::PERIODIC) {
                    amrex::Abort(
                        "BC is interior in direction " + std::to_string(dir) +
                        " but not periodic");
                }
                if (m_bc_hi[dir] == bc::PERIODIC) {
                    amrex::Abort(
                        "BC is interior in direction " + std::to_string(dir) +
                        " but not periodic");
                }
            }
        }

        const std::string vel_bc_key = "velocity_bc_type";
        bool has_vel_bc = false;
        // if it is velocity BC, make sure you have a velocity BC type
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            if ((m_bc_lo[dir] == bc::VELOCITY) ||
                (m_bc_hi[dir] == bc::VELOCITY)) {
                has_vel_bc = true;
            }
        }
        if (!(pp.contains(vel_bc_key.c_str())) && has_vel_bc) {
            amrex::Abort(
                "LBM::read_paramaters: velocity BC is used without specifying "
                "the type to be used");
        }
        pp.query(vel_bc_key.c_str(), m_velocity_bc_type);

        pp.get("ic_type", m_ic_type);

        pp.query("dx_outer", m_dx_outer);
        pp.query("dt_outer", m_dt_outer);

        pp.query("nu", m_nu);
        m_alpha = m_nu;
        pp.query("alpha", m_alpha);

        m_component_diffusivities.resize(m_n_components);
        for (int i = 0; i < m_n_components; ++i) {
            std::string diff_key = "diffusivity_component_" + std::to_string(i);
            m_component_diffusivities[i] = m_nu;
            pp.query(diff_key.c_str(), m_component_diffusivities[i]);
        }

        pp.query("save_streaming", m_save_streaming);
        pp.query("save_derived", m_save_derived);

        pp.query("compute_forces", m_compute_forces);
        pp.query("forces_file", m_forces_file);
        pp.query("clamp_component_densities", m_clamp_component_densities);
        pp.query("use_entropic_components", m_use_entropic_components);
        pp.query("use_entropic_f", m_use_entropic_f);

        pp.query("initial_temperature", m_initialTemperature);

        pp.query("body_is_isothermal", m_bodyIsIsothermal);
        pp.query("body_temperature", m_bodyTemperature);

        pp.query("is_fluid_fraction_threshold", m_is_fluid_fraction_threshold);

        // Reaction parameters
        pp.query("enable_reactions",  m_enable_reactions);
        pp.query("rxn_k_forward",     m_rxn_k_forward);
        pp.query("rxn_k_reverse",     m_rxn_k_reverse);
        pp.query("rxn_k_product",     m_rxn_k_product);

        // Timed catalyst injection
        pp.query("cat_inject_step",    m_cat_inject_step);
        pp.query("cat_inject_density", m_cat_inject_density);
        {
            amrex::Vector<amrex::Real> lo_tmp(AMREX_SPACEDIM, 0.0);
            amrex::Vector<amrex::Real> hi_tmp(AMREX_SPACEDIM, 0.0);
            pp.queryarr("cat_inject_box_lo", lo_tmp, 0, AMREX_SPACEDIM);
            pp.queryarr("cat_inject_box_hi", hi_tmp, 0, AMREX_SPACEDIM);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                m_cat_inject_box_lo[d] = lo_tmp[d];
                m_cat_inject_box_hi[d] = hi_tmp[d];
            }
        }

        pp.query("species_stats_file", m_species_stats_file);
    }

    // Moving body parameters
    {
        amrex::ParmParse pp("body");
        pp.query("is_moving", m_body_is_moving);
        
        // Use temporary vectors for queryarr, then copy to GpuArray
        amrex::Vector<amrex::Real> vel_tmp(3, 0.0);
        amrex::Vector<amrex::Real> omega_tmp(3, 0.0);
        amrex::Vector<amrex::Real> center_tmp(3, 0.0);
        
        pp.queryarr("velocity", vel_tmp, 0, 3);
        pp.queryarr("angular_velocity", omega_tmp, 0, 3);
        pp.queryarr("center", center_tmp, 0, 3);
        
        for (int i = 0; i < 3; ++i) {
            m_body_velocity[i] = vel_tmp[i];
            m_body_angular_velocity[i] = omega_tmp[i];
            m_body_center[i] = center_tmp[i];
        }
        
        if (m_body_is_moving) {
            amrex::Print() << "\n=== Moving Body Configuration ===" << std::endl;
            amrex::Print() << "Body velocity: (" << m_body_velocity[0] << ", "
                          << m_body_velocity[1] << ", " << m_body_velocity[2] << ")" << std::endl;
            amrex::Print() << "Angular velocity: (" << m_body_angular_velocity[0] << ", "
                          << m_body_angular_velocity[1] << ", " << m_body_angular_velocity[2] << ")" << std::endl;
            amrex::Print() << "Rotation center: (" << m_body_center[0] << ", "
                          << m_body_center[1] << ", " << m_body_center[2] << ")" << std::endl;
        }
    }

    // Free-surface parameters (Chiu & Lin 2011 conservative phase-field)
    {
        amrex::ParmParse pp("lbm");
        pp.query("free_surface",          m_free_surface);
        pp.query("free_surface_z",        m_free_surface_z);
        pp.query("free_surface_gamma",    m_phi_gamma_coeff);
        if (m_free_surface) {
            // Read reference density from the initial-condition block so that
            // FSLBM seeding, ABB, and repair thresholds scale with the actual
            // bulk density rather than assuming ρ = 1.
            amrex::ParmParse ppic("ic_constant");
            ppic.query("density", m_fslbm_rho_ref);
            pp.query("fslbm_sigma",          m_fslbm_sigma);
            pp.query("fslbm_contact_angle",  m_fslbm_contact_angle_deg);

            amrex::Print() << "\n=== Free Surface Configuration (FSLBM) ===" << std::endl;
            amrex::Print() << "  Interface z (LB cells)   : " << m_free_surface_z << std::endl;
            amrex::Print() << "  Reference density ρ_ref  : " << m_fslbm_rho_ref << std::endl;
            amrex::Print() << "  Surface tension σ (LB)   : " << m_fslbm_sigma
                           << (m_fslbm_sigma == 0.0 ? "  (flat interface)" : "") << std::endl;
            amrex::Print() << "  Contact angle θ (deg)    : " << m_fslbm_contact_angle_deg
                           << (std::abs(m_fslbm_contact_angle_deg - 90.0) < 0.01 ? "  (neutral wetting)" : "")
                           << std::endl;
        }
    }

    // Get geometry type for SDF reconstruction
    {
        amrex::ParmParse pp("eb2");
        pp.query("geom_type", m_body_geom_type);
        amrex::Print() << "Read body geom_type: '" << m_body_geom_type << "'" << std::endl;
    }

    {
        amrex::ParmParse pp("lbm");
    // threshold for converting fractional mask to integer is_fluid
    pp.query("is_fluid_fraction_threshold", m_is_fluid_fraction_threshold);

        pp.query("adiabatic_exponent", m_adiabaticExponent);
        pp.query("mean_molecular_mass", m_m_bar);

        m_speedOfSound_Ref = std::sqrt(
            m_adiabaticExponent * (m_R_u / m_m_bar) * m_initialTemperature);

        m_mesh_speed = m_dx_outer / m_dt_outer;
        m_cs = m_mesh_speed / constants::ROOT3;

        m_cs_2 = m_cs * m_cs;
    }

    // ---------------------------------------------------------------
    // Lagrangian bubble (kLa) parameters
    // ---------------------------------------------------------------
    {
        amrex::ParmParse pp("lbm");
        pp.query("enable_bubbles", m_enable_bubbles);
    }

    if (m_enable_bubbles) {
        BubbleManager::read_params(m_bubble_params);

        // Physical unit conversions — must be set explicitly in input file.
        // lbm.dx_outer and lbm.dt_outer are dimensionless LB units (= 1.0);
        // dx_phys and dt_phys must be the actual SI values.
        m_bubble_params.dx_phys = m_dx_outer;  // overridden below if lbm.dx_phys provided
        m_bubble_params.dt_phys = m_dt_outer;  // overridden below if lbm.dt_phys provided
        m_bubble_params.nu_lb   = m_nu;
        {
            amrex::ParmParse pp_lbm("lbm");
            // Read actual SI values (meters, seconds) — REQUIRED for physical accuracy.
            // dx_phys = physical cell size [m];  dt_phys = physical time step [s].
            pp_lbm.query("dx_phys", m_bubble_params.dx_phys);
            pp_lbm.query("dt_phys", m_bubble_params.dt_phys);
        }

        // Concentration reference scale: 1 LB_rho ≡ m_bubble_o2_C_ref mol/m³
        {
            amrex::ParmParse pp("bubble");
            pp.query("O2_concentration_reference", m_bubble_o2_C_ref);
        }
        // Propagate C_ref into BubbleParams so deposit_o2_sources can use it
        m_bubble_params.C_ref = m_bubble_o2_C_ref;

        amrex::Print() << "[BubbleManager] Bubble physics enabled.\n"
                       << "  dx_phys = " << m_bubble_params.dx_phys << " m\n"
                       << "  dt_phys = " << m_bubble_params.dt_phys << " s\n"
                       << "  O2 C_ref = " << m_bubble_o2_C_ref << " mol/m3 per LB_rho\n";
    }
}

void LBM::read_tagging_parameters()
{
    BL_PROFILE("LBM::read_ragging_parameters()");

    const std::string tag_prefix = "tagging";
    amrex::ParmParse pp(tag_prefix);
    amrex::Vector<std::string> refinement_indicators;
    pp.queryarr(
        "refinement_indicators", refinement_indicators, 0,
        pp.countval("refinement_indicators"));
    for (const auto& refinement_indicator : refinement_indicators) {
        const std::string ref_prefix = tag_prefix + "." + refinement_indicator;
        amrex::ParmParse ppr(ref_prefix);

        // Tag a given box
        amrex::RealBox realbox;
        if (ppr.countval("in_box_lo") > 0) {
            amrex::Vector<amrex::Real> box_lo(AMREX_SPACEDIM);
            amrex::Vector<amrex::Real> box_hi(AMREX_SPACEDIM);
            ppr.getarr("in_box_lo", box_lo, 0, static_cast<int>(box_lo.size()));
            ppr.getarr("in_box_hi", box_hi, 0, static_cast<int>(box_hi.size()));
            realbox = amrex::RealBox(box_lo.data(), box_hi.data());
        }

        amrex::AMRErrorTagInfo info;

        if (realbox.ok()) {
            info.SetRealBox(realbox);
        }

        if (ppr.countval("start_time") > 0) {
            amrex::Real min_time;
            ppr.get("start_time", min_time);
            info.SetMinTime(min_time);
        }

        if (ppr.countval("end_time") > 0) {
            amrex::Real max_time;
            ppr.get("end_time", max_time);
            info.SetMaxTime(max_time);
        }

        if (ppr.countval("max_level") > 0) {
            int tag_max_level;
            ppr.get("max_level", tag_max_level);
            info.SetMaxLevel(tag_max_level);
        }

        bool itexists = false;
        if (ppr.countval("value_greater") > 0) {
            amrex::Real value;
            ppr.get("value_greater", value);
            std::string field;
            ppr.get("field_name", field);
            m_err_tags.push_back(
                amrex::AMRErrorTag(
                    value, amrex::AMRErrorTag::GREATER, field, info));
            itexists = check_field_existence(field);
        } else if (ppr.countval("value_less") > 0) {
            amrex::Real value;
            ppr.get("value_less", value);
            std::string field;
            ppr.get("field_name", field);
            m_err_tags.push_back(
                amrex::AMRErrorTag(
                    value, amrex::AMRErrorTag::LESS, field, info));
            itexists = check_field_existence(field);
        } else if (ppr.countval("adjacent_difference_greater") > 0) {
            amrex::Real value;
            ppr.get("adjacent_difference_greater", value);
            std::string field;
            ppr.get("field_name", field);
            m_err_tags.push_back(
                amrex::AMRErrorTag(
                    value, amrex::AMRErrorTag::GRAD, field, info));
            itexists = check_field_existence(field);
        } else if (realbox.ok()) {
            m_err_tags.push_back(amrex::AMRErrorTag(info));
            itexists = true;
        } else {
            amrex::Abort(
                "LBM::read_tagging_parameters(): unrecognized refinement "
                "indicator for " +
                refinement_indicator);
        }

        if (!itexists) {
            amrex::Error(
                "LBM::read_tagging_parameters(): unknown variable field for "
                "tagging "
                "criteria " +
                refinement_indicator);
        }
    }
}

void LBM::evolve()
{
    BL_PROFILE("LBM::evolve()");

    amrex::Real cur_time = m_ts_new[0];
    int last_plot_file_step = 0;

    for (int step = m_isteps[0]; step < m_max_step && cur_time < m_stop_time;
         ++step) {
        compute_dt();

        if (m_print_int > 0 && step % m_print_int == 0) {
            amrex::Print() << "\n==============================================="
                              "==============================="
                           << std::endl;
            amrex::Print() << "Step: " << step << " dt : " << m_dts[0]
                           << " time: " << cur_time << " to " << cur_time + m_dts[0]
                           << std::endl;
        }

        m_fillpatch_op->fillpatch(0, cur_time, m_f[0]);
        for (int i = 0; i < m_n_components; ++i) {
            m_component_fillpatch_ops[i]->fillpatch(0, cur_time, m_component_lattices[i][0]);
        }

        m_fillpatch_g_op->fillpatch(0, cur_time, m_g[0]);

        time_step(0, cur_time, 1);

        post_time_step();

        cur_time += m_dts[0];

        // sync up time
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_ts_new[lev] = cur_time;
        }

        if (m_plot_int > 0 && (step + 1) % m_plot_int == 0) {
            last_plot_file_step = step + 1;
            write_plot_file();
        }

        if (m_chk_int > 0 && (step + 1) % m_chk_int == 0) {
            write_checkpoint_file();
        }

        if (cur_time >= m_stop_time - 1.e-6 * m_dts[0]) {
            break;
        }
    }
    if (m_plot_int > 0 && m_isteps[0] > last_plot_file_step) {
        write_plot_file();
    }
    close_species_stats_file();
    close_forces_file();
}

// advance a level by dt
// includes a recursive call for finer levels
void LBM::time_step(const int lev, const amrex::Real time, const int iteration)
{
    BL_PROFILE("LBM::time_step()");
    if (m_regrid_int > 0) // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static amrex::Vector<int> last_regrid_step(max_level + 1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && m_isteps[lev] > last_regrid_step[lev]) {
            if (m_isteps[lev] % m_regrid_int == 0) {
                // regrid could add newly refine levels (if finest_level <
                // max_level) so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = m_isteps[k];
                }

                // if there are newly created levels, set the time step
                // dt gets halved here
                for (int k = old_finest + 1; k <= finest_level; ++k) {
                    m_dts[k] = m_dts[k - 1] / MaxRefRatio(k - 1);
                }
                if (amrex::ParallelDescriptor::IOProcessor()) {
                    amrex::Print()
                        << "Grid summary after regrid: " << std::endl;
                    printGridSummary(amrex::OutStream(), 0, finest_level);
                }
            }
        }
    }

    if (Verbose() != 0) {
        amrex::Print() << "[Level " << lev << " step " << m_isteps[lev] + 1
                       << "] ";
        amrex::Print() << "Advance with time = " << m_ts_new[lev]
                       << " dt = " << m_dts[lev] << std::endl;
    }

    if (lev < finest_level) {
        m_fillpatch_op->fillpatch(lev + 1, m_ts_new[lev + 1], m_f[lev + 1]);
        for (int i = 0; i < m_n_components; ++i) {
            m_component_fillpatch_ops[i]->fillpatch(lev + 1, m_ts_new[lev + 1], m_component_lattices[i][lev + 1]);
        }

        m_fillpatch_g_op->fillpatch(lev + 1, m_ts_new[lev + 1], m_g[lev + 1]);

        for (int i = 1; i <= m_nsubsteps[lev + 1]; ++i) {
            m_fillpatch_op->physbc(lev + 1, m_ts_new[lev + 1], m_f[lev + 1]);
            for (int c = 0; c < m_n_components; ++c) {
                m_component_fillpatch_ops[c]->physbc(lev + 1, m_ts_new[lev + 1], m_component_lattices[c][lev + 1]);
            }

            m_fillpatch_g_op->physbc(lev + 1, m_ts_new[lev + 1], m_g[lev + 1]);

            time_step(lev + 1, time + (i - 1) * m_dts[lev + 1], i);
        }
    }

    advance(lev, time, m_dts[lev], iteration, m_nsubsteps[lev]);

    ++m_isteps[lev];

    if (Verbose() != 0) {
        amrex::Print() << "[Level " << lev << " step " << m_isteps[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells"
                       << std::endl;
    }
}

void LBM::advance(
    const int lev,
    const amrex::Real /*time*/,
    const amrex::Real dt_lev,
    const int /*iteration*/,
    const int /*ncycle*/)
{
    BL_PROFILE("LBM::advance()");

    m_ts_old[lev] = m_ts_new[lev]; // old time is now current time (time)
    m_ts_new[lev] += dt_lev;       // new time is ahead

    // --- O2 mass tracking diagnostic (per-step) --- DISABLED for performance
    // To re-enable: uncomment this block and the corresponding measurement blocks below.
#if 0
    auto sum_comp0_mass = [&]() -> amrex::Real {
        if (m_n_components < 1) return 0.0;
        amrex::Real total = 0.0;
        for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
            total += m_component_lattices[0][lev].sum(q);
        }
        return total;
    };
    const bool o2_diag = (m_n_components > 0 && m_isteps[lev] >= 100);
    amrex::Real o2_mass_A0 = 0.0;
    if (o2_diag) {
        o2_mass_A0 = sum_comp0_mass();
    }
#endif

    // Update moving body position and reconstruct fluid/solid boundaries
    if (m_body_is_moving) {
        reconstruct_body_sdf(lev, m_ts_new[lev]);
        
        // Fill ghost cells BEFORE refill_and_spill so that spill algorithm
        // doesn't read stale/uninitialized values from ghost cells
        m_f[lev].FillBoundary(Geom(lev).periodicity());
        m_g[lev].FillBoundary(Geom(lev).periodicity());
        for (int i = 0; i < m_n_components; ++i) {
            m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
        }
        
        refill_and_spill(lev);
    }

#if 0  // O2 mass diagnostic — disabled for performance
    amrex::Real o2_mass_A = 0.0;
    if (o2_diag) {
        o2_mass_A = sum_comp0_mass();
        amrex::Real loss_spill = (o2_mass_A0 > 1e-20) ? (o2_mass_A0 - o2_mass_A) / o2_mass_A0 : 0.0;
        if (loss_spill > 0.001 || m_isteps[lev] % 200 == 0) {
            amrex::Print() << "[O2_mass step=" << m_isteps[lev]
                           << "] A0(start)=" << o2_mass_A0
                           << " A(after_spill)=" << o2_mass_A
                           << " loss_spill=" << loss_spill*100 << "%\n";
        }
    }
#endif

    // Free-surface advance: FSLBM (Körner 2005) replaces both advance_phi and
    // stream(lev, m_f).  When m_free_surface is false the standard stream() runs.
    if (m_free_surface) {
        fslbm_advance_surface(lev);  // streams m_f + updates φ + converts cells
    } else {
        stream(lev, m_f);
    }

#if 0  // O2 mass diagnostic — disabled for performance
    amrex::Real o2_mass_B = 0.0;
    if (o2_diag) {
        o2_mass_B = sum_comp0_mass();
    }
#endif

    for (int i = 0; i < m_n_components; ++i) {
        stream(lev, m_component_lattices[i]);
    }

#if 0  // O2 mass diagnostic — disabled for performance
    amrex::Real o2_mass_C = 0.0;
    if (o2_diag) {
        o2_mass_C = sum_comp0_mass();
        // Print detailed if significant loss detected at any stage
        amrex::Real loss_B = (o2_mass_A > 1e-20) ? (o2_mass_A - o2_mass_B) / o2_mass_A : 0.0;
        amrex::Real loss_C = (o2_mass_B > 1e-20) ? (o2_mass_B - o2_mass_C) / o2_mass_B : 0.0;
        if (loss_B > 0.01 || loss_C > 0.01 || m_isteps[lev] % 200 == 0) {
            amrex::Print() << "[O2_mass step=" << m_isteps[lev]
                           << "] B(after_fslbm)=" << o2_mass_B
                           << " C(after_stream)=" << o2_mass_C
                           << " loss_fslbm=" << loss_B*100 << "%"
                           << " loss_stream=" << loss_C*100 << "%\n";
        }
    }
#endif

    stream(lev, m_g);

    if (lev < finest_level) {
        average_down_to(lev, amrex::IntVect(1));
    }

    // Clamp negative component densities BEFORE collision so the collision
    // step always acts on clean (non-negative) populations.
    // Activated via lbm.clamp_component_densities = 1 in the input file.
    if (m_clamp_component_densities) {
        clamp_negative_component_densities(lev);
    }

    collide(lev);

#if 0  // O2 mass diagnostic — disabled for performance
    if (o2_diag) {
        amrex::Real o2_mass_D = sum_comp0_mass();
        amrex::Real loss_D = (o2_mass_C > 1e-20) ? (o2_mass_C - o2_mass_D) / o2_mass_C : 0.0;
        if (loss_D > 0.01 || m_isteps[lev] % 200 == 0) {
            amrex::Print() << "[O2_mass step=" << m_isteps[lev]
                           << "] D(after_collide)=" << o2_mass_D
                           << " loss_collide=" << loss_D*100 << "%\n";
        }
    }
#endif

    // Catalyst injection: executed exactly once on level 0 when the
    // configured step is reached.  After injection the populations are
    // filled and m_cat_inject_done prevents any repeated application.
    if (lev == 0) {
        apply_timed_catalyst_injection(lev);
    }

    // Operator-split chemistry: add/remove mass from the four scalar
    // fields according to the two-step catalytic reaction kinetics.
    if (m_enable_reactions && m_n_components >= 4) {
        apply_reaction_source_terms(lev);
    }

    // ------------------------------------------------------------------
    // Lagrangian bubble physics — two-phase kLa mass transfer
    // Execute only on the base level to keep a single particle container.
    // ------------------------------------------------------------------
    if (m_enable_bubbles && lev == 0) {
        // Physical time in seconds (m_ts_new is in LB steps, dt_phys is s/step)
        const amrex::Real phys_time = m_ts_new[lev] * m_bubble_params.dt_phys;

        // Temporary MultiFabs for bubble↔fluid coupling (zeroed each step)
        amrex::MultiFab bubble_force(
            m_f[lev].boxArray(), m_f[lev].DistributionMap(), 3, 0);
        amrex::MultiFab o2_src(
            m_f[lev].boxArray(), m_f[lev].DistributionMap(), 1, 0);
        bubble_force.setVal(0.0);
        o2_src.setVal(0.0);

        // Sparger injection (every step)
        // Must pass physical seconds per step, not the dimensionless LB m_dt_outer.
        m_bubbles.inject_bubbles(m_bubble_params.dt_phys);

        // Determine O2 concentration MultiFab (component 0 if available)
        // A valid kLa run requires at least 1 component for dissolved O2.
        if (m_n_components < 1) {
            amrex::Abort("lbm.enable_bubbles = 1 requires lbm.n_components >= 1 "
                         "(component 0 = dissolved O2).");
        }

        // Precompute macroscopic O2 density (sum over all N_MICRO_STATES populations)
        // so that BubbleManager::deposit_o2_sources can interpolate the correct C_f.
        // BUG FIX: previously passed m_component_lattices[0][lev] directly and
        // trilinear_interp used comp=0 (q=0 rest population only, ≈ w_0 × rho_O2 ≈ rho_O2/3),
        // underestimating C_f by ~3× and overestimating the driving force.
        // Use MultiFab::Add in a loop to avoid __device__ lambdas in a private method.
        amrex::MultiFab rho_o2(
            m_f[lev].boxArray(), m_f[lev].DistributionMap(), 1, 0);
        rho_o2.setVal(0.0);
        for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
            amrex::MultiFab::Add(rho_o2, m_component_lattices[0][lev], q, 0, 1, 0);
        }

        // Advance bubbles: forces, Verlet integration, mass transfer
        m_bubbles.advance(
            dt_lev,
            m_macrodata[lev],
            m_derived[lev],
            rho_o2,       // 1-component macroscopic O2 density [LB_rho]
            Geom(lev),
            bubble_force,
            o2_src,
            phys_time,
            // BUG FIX: was m_is_fluid_fraction[lev] (EB SDF, < 0.5 inside impeller EB cells)
            // which falsely removed live bubbles passing through the impeller swept volume.
            // Correct field: m_phi_fslbm[lev] (gas-liquid phase field, < 0.5 only in gas headspace).
            m_free_surface ? &m_phi_fslbm[lev] : nullptr);

        // Coalescence check (every coal_interval steps)
        ++m_bubble_step_counter;
        if (m_bubble_params.enable_coalescence &&
            m_bubble_step_counter % m_bubble_params.coal_interval == 0) {
            m_bubbles.do_coalescence(phys_time);
        }

        // Apply bubble body force to fluid distributions (He-Luo scheme)
        // Diagnostic: print max force magnitude to catch anomalous deposits.
        if (m_print_int > 0 && m_isteps[lev] % m_print_int == 0) {
            const amrex::Real Fx_max = bubble_force.norm0(0);
            const amrex::Real Fy_max = bubble_force.norm0(1);
            const amrex::Real Fz_max = bubble_force.norm0(2);
            amrex::Print() << "[bubble_force step=" << m_isteps[lev]
                           << "] max|Fx|=" << Fx_max
                           << "  max|Fy|=" << Fy_max
                           << "  max|Fz|=" << Fz_max << "\n";
        }
        apply_bubble_body_force(lev, bubble_force);

        // Apply O2 source to component-0 lattice
        if (m_n_components > 0) {
            // --- O2 diagnostic: print every print_int steps ---
            if (m_print_int > 0 && m_isteps[lev] % m_print_int == 0) {
                const amrex::Real src_max = o2_src.norm0();
                const amrex::Real rho_o2_before = m_component_lattices[0][lev].norm0();
                amrex::Print() << "[O2_debug step=" << m_isteps[lev]
                               << "] o2_src.norm0=" << src_max
                               << "  rho_O2_before=" << rho_o2_before << "\n";
            }
            apply_bubble_o2_source(lev, o2_src);
            if (m_print_int > 0 && m_isteps[lev] % m_print_int == 0) {
                const amrex::Real rho_o2_after = m_component_lattices[0][lev].norm0();
                amrex::Print() << "[O2_debug step=" << m_isteps[lev]
                               << "] rho_O2_after=" << rho_o2_after << "\n";
            }
        }

        // Statistics output
        if (m_bubble_params.stats_int > 0 &&
            m_isteps[lev] % m_bubble_params.stats_int == 0) {
            m_bubbles.write_stats(m_isteps[lev], phys_time);
        }
    }
}

void LBM::post_time_step()
{
    BL_PROFILE("LBM::post_time_step()");

    for (int lev = 0; lev <= finest_level; ++lev) {
        compute_derived(lev);
    }

    compute_eb_forces();

    // Write per-step mean species concentrations when reactions are enabled.
    if (m_enable_reactions && m_n_components >= 4) {
        write_species_stats();
    }
}

// Stream the information to the neighbor particles
void LBM::stream(const int lev, amrex::Vector<amrex::MultiFab>& fs)
{
    BL_PROFILE("LBM::stream()");

    amrex::MultiFab f_star(
        boxArray(lev), DistributionMap(lev), constants::N_MICRO_STATES,
        fs[lev].nGrow(), amrex::MFInfo(), *(m_factory[lev]));
    f_star.setVal(0.0);

    auto const& fs_arrs = f_star.arrays();
    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& f_arrs = fs[lev].const_arrays();

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& bounce_dirs = stencil.bounce_dirs;
    amrex::ParallelFor(
        fs[lev], fs[lev].nGrowVect(), constants::N_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k),
            int q) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const auto& ev = evs[q];
            const amrex::IntVect ivn(iv + ev);
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) == 1) {
                const auto f_arr = f_arrs[nbx];
                const auto fs_arr = fs_arrs[nbx];
                const auto& lb = amrex::lbound(f_arr);
                const auto& ub = amrex::ubound(f_arr);
                const amrex::Box fbox(
                    amrex::IntVect(AMREX_D_DECL(lb.x, lb.y, lb.z)),
                    amrex::IntVect(AMREX_D_DECL(ub.x, ub.y, ub.z)));
                if (fbox.contains(ivn)) {
                    if (is_fluid_arrs[nbx](ivn, lbm::constants::IS_FLUID_IDX) != 0) {
                        fs_arr(ivn, q) = f_arr(iv, q);
                    } else {
                        fs_arr(iv, bounce_dirs[q]) = f_arr(iv, q);
                    }
                }
            }
        });
    amrex::Gpu::synchronize();

    amrex::MultiFab::Copy(
        fs[lev], f_star, 0, 0, constants::N_MICRO_STATES, fs[lev].nGrowVect());
    fs[lev].FillBoundary(Geom(lev).periodicity());
}

// Clamp negative component densities to zero
void LBM::clamp_negative_component_densities(const int lev)
{
    BL_PROFILE("LBM::clamp_negative_component_densities()");
    
    // Component lattices can develop negative populations through:
    // 1. Streaming - when cells with zero density stream from boundary ghost cells,
    //    floating point errors can create small negatives
    // 2. Collision - relaxation of near-zero density toward equilibrium can produce negatives
    // 
    // Unlike m_f and m_g which are always initialized with positive equilibrium distributions,
    // component lattices start with zeros in most of the domain and only gain mass through
    // initial conditions in specific regions. This makes them susceptible to negative values
    // during streaming and collision operations.
    //
    // This clamp runs on ALL cells (interior, ghost, solid) to prevent negatives from
    // propagating between timesteps via ghost cells.
    
    for (int c = 0; c < m_n_components; ++c) {
        auto const& f_comp_arrs = m_component_lattices[c][lev].arrays();
        
        amrex::ParallelFor(
            m_component_lattices[c][lev],
            m_component_lattices[c][lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Compute total density
                amrex::Real rho_comp = 0.0;
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    rho_comp += f_comp_arrs[nbx](i, j, k, q);
                }
                
                // Zero out all populations if total density is negative
                if (rho_comp < 0.0) {
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        f_comp_arrs[nbx](i, j, k, q) = 0.0;
                    }
                }
            });
        amrex::Gpu::synchronize();
    }
}

// Collide the particles
void LBM::collide(const int lev)
{
    BL_PROFILE("LBM::collide()");

    f_to_macrodata(lev);

    compute_q_corrections(lev);

    macrodata_to_equilibrium(lev);

    relax_f_to_equilibrium(lev);
}

// convert macrodata to equilibrium.
void LBM::macrodata_to_equilibrium(const int lev)
{
    BL_PROFILE("LBM::macrodata_to_equilibrium()");
    AMREX_ASSERT(m_macrodata[lev].nGrow() >= m_eq[lev].nGrow());
    auto const& md_arrs = m_macrodata[lev].const_arrays();
    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& eq_arrs = m_eq[lev].arrays();
    auto const& eq_arrs_g = m_eq_g[lev].arrays();
    const amrex::Real l_mesh_speed = m_mesh_speed;

    AMREX_ASSERT(m_macrodata[lev].nGrow() > m_derived[lev].nGrow());
    auto const& d_arrs = m_derived[lev].const_arrays();

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& weight = stencil.weights;

    const amrex::RealVect zero_vec = {AMREX_D_DECL(0.0, 0.0, 0.0)};
    const amrex::Real specific_gas_constant = m_R_u / m_m_bar;
    const amrex::Real cv = specific_gas_constant / (m_adiabaticExponent - 1.0);
    const amrex::Real nu = m_nu;
    const amrex::Real dt = m_dts[lev];
    const amrex::Real alpha = m_alpha;
    const amrex::Real theta0 = stencil::Stencil::THETA0;

    amrex::ParallelFor(
        m_eq[lev], m_eq[lev].nGrowVect(), constants::N_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k),
            int q) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) == 1) {

                const auto md_arr = md_arrs[nbx];
                const auto eq_arr = eq_arrs[nbx];
                const auto eq_arr_g = eq_arrs_g[nbx];
                const auto d_arr = d_arrs[nbx];

                const amrex::Real rho = md_arr(iv, constants::RHO_IDX);
                const amrex::RealVect vel = {AMREX_D_DECL(
                    md_arr(iv, constants::VELX_IDX),
                    md_arr(iv, constants::VELY_IDX),
                    md_arr(iv, constants::VELZ_IDX))};

                const amrex::Real two_rho_e =
                    md_arr(iv, constants::TWO_RHO_E_IDX);

                const amrex::Real wt = weight[q];

                const auto& ev = evs[q];

                const amrex::Real temperature =
                    md_arr(iv, constants::TEMPERATURE_IDX);

                const amrex::Real omega =
                    1.0 /
                    (nu / (specific_gas_constant * temperature * dt) + 0.5);
                const amrex::Real omega_one =
                    1.0 /
                    (alpha / (specific_gas_constant * temperature * dt) + 0.5);
                const amrex::Real omega_one_by_omega = omega_one / omega;
                const amrex::Real omega_corr =
                    (2.0 - omega) / (2.0 * omega * rho);
                
                const amrex::Real pxx_ext =
                    vel[0] * vel[0] + specific_gas_constant * temperature +
                    dt * (omega_corr)*d_arr(iv, constants::D_Q_CORR_X_IDX);
                const amrex::Real pyy_ext =
                    vel[1] * vel[1] + specific_gas_constant * temperature +
                    dt * (omega_corr)*d_arr(iv, constants::D_Q_CORR_Y_IDX);
                const amrex::Real pzz_ext = AMREX_D_PICK(
                    0.0, 0.0,
                    vel[2] * vel[2] + specific_gas_constant * temperature +
                        dt * (omega_corr)*d_arr(iv, constants::D_Q_CORR_Z_IDX));

                eq_arr(iv, q) = set_extended_equilibrium_value(
                    rho, vel, pxx_ext, pyy_ext, pzz_ext, l_mesh_speed, wt, ev);

                amrex::Real AMREX_D_DECL(qx_eq = 0.0, qy_eq = 0.0, qz_eq = 0.0);
                amrex::Real rxx_eq(0.0), ryy_eq(0.0), rzz_eq(0.0), rxy_eq(0.0),
                    rxz_eq(0.0), ryz_eq(0.0);

                amrex::RealVect heat_flux = {AMREX_D_DECL(0.0, 0.0, 0.0)};
                get_equilibrium_moments(
                    rho, vel, two_rho_e, cv, specific_gas_constant, heat_flux,
                    rxx_eq, ryy_eq, rzz_eq, rxy_eq, rxz_eq, ryz_eq);

                qx_eq = heat_flux[0];
                qy_eq = heat_flux[1];
                AMREX_3D_ONLY(qz_eq = heat_flux[2]);

                const amrex::Real pxx = md_arr(iv, constants::PXX_IDX);
                const amrex::Real pyy = md_arr(iv, constants::PYY_IDX);
                AMREX_3D_ONLY(
                    const amrex::Real pzz = md_arr(iv, constants::PZZ_IDX));
                const amrex::Real pxy = md_arr(iv, constants::PXY_IDX);
                AMREX_3D_ONLY(
                    const amrex::Real pxz = md_arr(iv, constants::PXZ_IDX));
                AMREX_3D_ONLY(
                    const amrex::Real pyz = md_arr(iv, constants::PYZ_IDX));

                const amrex::Real qx = md_arr(iv, constants::QX_IDX);
                const amrex::Real qy = md_arr(iv, constants::QY_IDX);
                AMREX_3D_ONLY(
                    const amrex::Real qz = md_arr(iv, constants::QZ_IDX));

                qx_eq *= omega_one_by_omega;
                qy_eq *= omega_one_by_omega;
                AMREX_3D_ONLY(qz_eq *= omega_one_by_omega);

                qx_eq += (1.0 - omega_one_by_omega) *
                         (qx AMREX_D_TERM(
                              -2.0 * vel[0] * pxx, -2.0 * vel[1] * pxy,
                              -2.0 * vel[2] * pxz) -
                          vel[0] * dt * d_arr(iv, constants::D_Q_CORR_X_IDX));

                qy_eq += (1.0 - omega_one_by_omega) *
                         (qy AMREX_D_TERM(
                              -2.0 * vel[0] * pxy, -2.0 * vel[1] * pyy,
                              -2.0 * vel[2] * pyz) -
                          vel[1] * dt * d_arr(iv, constants::D_Q_CORR_Y_IDX));

                AMREX_3D_ONLY(
                    qz_eq +=
                    (1.0 - omega_one_by_omega) *
                    (qz - 2.0 * vel[0] * pxz - 2.0 * vel[1] * pyz -
                     2.0 * vel[2] * pzz -
                     vel[2] * dt * d_arr(iv, constants::D_Q_CORR_Z_IDX)));

                amrex::RealVect heat_flux_mrt = {
                    AMREX_D_DECL(qx_eq, qy_eq, qz_eq)};

                amrex::GpuArray<amrex::Real, 6> flux_of_heat_flux = {
                    rxx_eq, ryy_eq, rzz_eq, rxy_eq, rxz_eq, ryz_eq};

                eq_arr_g(iv, q) = set_extended_grad_expansion_generic(
                    two_rho_e, heat_flux_mrt, flux_of_heat_flux, l_mesh_speed,
                    wt, ev, theta0, zero_vec, 1.0);
            }
        });
    amrex::Gpu::synchronize();
}

// Relax the particles toward the equilibrium state.
void LBM::relax_f_to_equilibrium(const int lev)
{
    BL_PROFILE("LBM::relax_f_to_equilibrium()");
    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& eq_arrs = m_eq[lev].const_arrays();
    auto const& eq_arrs_g = m_eq_g[lev].const_arrays();
    auto const& f_arrs = m_f[lev].arrays();
    auto const& g_arrs = m_g[lev].arrays();
    auto const& md_arrs = m_macrodata[lev].arrays();

    amrex::Real specific_gas_constant = (m_R_u / m_m_bar);
    amrex::Real nu = m_nu;
    amrex::Real dt = m_dts[lev];

    const bool body_is_isothermal = m_bodyIsIsothermal;
    const bool use_entropic_f     = m_use_entropic_f;

    const amrex::Real l_mesh_speed = m_mesh_speed;
    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& weight = stencil.weights;

    // Per-direction BGK pass: always updates g; updates f only when entropic is OFF.
    amrex::ParallelFor(
        m_f[lev], m_eq[lev].nGrowVect(), constants::N_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k),
            int q) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) == 1) {
                const auto f_arr = f_arrs[nbx];
                const auto eq_arr = eq_arrs[nbx];
                const auto md_arr = md_arrs[nbx];

                const auto g_arr = g_arrs[nbx];
                const auto eq_arr_g = eq_arrs_g[nbx];

                amrex::Real temperature =
                    md_arr(iv, constants::TEMPERATURE_IDX);
                amrex::Real omega =
                    1.0 /
                    (nu / (specific_gas_constant * temperature * dt) + 0.5);

                // f is updated here only for plain BGK; the entropic path
                // handles f in a separate cell-loop below.
                if (!use_entropic_f) {
                    f_arr(iv, q) += omega * (eq_arr(iv, q) - f_arr(iv, q));
                }

                g_arr(iv, q) += omega * (eq_arr_g(iv, q) - g_arr(iv, q));

                if (body_is_isothermal) {
                    if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_SIDE_IDX) == 1) {
                        g_arr(iv, q) = eq_arr_g(iv, q);
                    }
                }
            }
        });
    amrex::Gpu::synchronize();  // catch any CUDA error before bubble CPU code runs

    // --- Entropic alpha solve for m_f (Ansumali & Karlin, Phys. Rev. E 2002) ---
    // Finds alpha* in (0, 2] s.t. H(f + alpha*(f_eq - f)) = H(f), where
    //   H(f) = sum_q f_q * ln(f_q / f_ref_q)
    //   f_ref_q = f^eq(rho, 0, T)  (zero-velocity Maxwellian reference).
    //
    // alpha_use hierarchy:
    //   Newton succeeds  :  min(omega, alpha*)   [tightest H-theorem bound]
    //   Negative pops    :  min(omega, 1.0)      [safe fallback, monotone]
    //   Disabled         :  (this block not executed)
    //
    // Controlled by lbm.use_entropic_f (default = 0).
    if (use_entropic_f) {
        constexpr int NQ = constants::N_MICRO_STATES;
        amrex::ParallelFor(
            m_f[lev], m_eq[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) != 1) {
                    return;
                }

                const auto f_arr  = f_arrs[nbx];
                const auto eq_arr = eq_arrs[nbx];
                const auto md_arr = md_arrs[nbx];

                const amrex::Real temperature =
                    md_arr(iv, constants::TEMPERATURE_IDX);
                const amrex::Real omega =
                    1.0 / (nu / (specific_gas_constant * temperature * dt) + 0.5);
                const amrex::Real p_by_rho = specific_gas_constant * temperature;

                // BGK target (q-corrected, already stored in m_eq)
                amrex::GpuArray<amrex::Real, NQ> eq_all;
                for (int q = 0; q < NQ; ++q) {
                    eq_all[q] = eq_arr(iv, q);
                }

                // Zero-velocity unit-density reference for the H-function:
                //   f_ref[q] = f^eq(1, 0, T)  — identical to eq_unit_arr(iv, q+NQ)
                //   used for the component entropic solve.  Unit density makes the
                //   logarithm independent of rho (matching the component convention).
                const amrex::RealVect zero_vel = {AMREX_D_DECL(0.0, 0.0, 0.0)};
                amrex::GpuArray<amrex::Real, NQ> f_ref;
                for (int q = 0; q < NQ; ++q) {
                    f_ref[q] = set_extended_equilibrium_value(
                        1.0, zero_vel, p_by_rho, p_by_rho, p_by_rho,
                        l_mesh_speed, weight[q], evs[q]);
                }

                // Safe fallback: min(omega, 1.0) guarantees monotone approach
                amrex::Real alpha_use = amrex::min(omega, 1.0);

                // Attempt Newton only when all pre-collision pops are positive
                amrex::Real H0 = 0.0;
                bool all_positive = true;
                for (int q = 0; q < NQ; ++q) {
                    const amrex::Real fq = f_arr(iv, q);
                    if (fq <= 0.0) { all_positive = false; break; }
                    H0 += fq * log(fq / f_ref[q]);
                }

                if (all_positive) {
                    // Newton: g(alpha) = H(f + alpha*s) - H0 = 0
                    amrex::Real alpha = 2.0;
                    for (int iter = 0; iter < 10; ++iter) {
                        amrex::Real g = -H0, dg = 0.0;
                        bool fhat_positive = true;
                        for (int q = 0; q < NQ; ++q) {
                            const amrex::Real sq   = eq_all[q] - f_arr(iv, q);
                            const amrex::Real fhat = f_arr(iv, q) + alpha * sq;
                            if (fhat <= 0.0) { fhat_positive = false; break; }
                            const amrex::Real ln_ratio = log(fhat / f_ref[q]);
                            g  += fhat * ln_ratio;
                            dg += sq * (ln_ratio + 1.0);
                        }
                        if (!fhat_positive || fabs(dg) < 1.0e-14) { break; }
                        alpha -= g / dg;
                        alpha = amrex::min(2.0, amrex::max(0.0, alpha));
                        if (fabs(g) < 1.0e-12 * (fabs(H0) + 1.0e-30)) { break; }
                    }
                    alpha_use = amrex::min(omega, alpha);
                }

                // Apply entropic collision to f
                for (int q = 0; q < NQ; ++q) {
                    f_arr(iv, q) += alpha_use * (eq_all[q] - f_arr(iv, q));
                }
            });
        amrex::Gpu::synchronize();
    }

    const bool use_entropic_components = m_use_entropic_components;

    // --------------------------------------------------------------------------
    // Pre-compute two unit-density equilibrium shapes once, shared across all
    // components.  set_extended_equilibrium_value is linear in rho.
    //
    // Layout (2 * N_MICRO_STATES components):
    //   [0 .. N)   : eq_flow(q)  = f^eq(1, u_local, T_local)
    //                Used as the BGK target: f^eq_comp = rho_comp * eq_flow(q)
    //   [N .. 2N)  : eq_ref(q)   = f^eq(1, 0, T_local)
    //                Used as the H-function reference in the entropic solve.
    //                Generalises the lattice weights w_q (which equal eq_ref
    //                only at the isothermal point cs^2 = 1/3).
    // --------------------------------------------------------------------------
    const int NQ = constants::N_MICRO_STATES;
    amrex::MultiFab eq_unit(
        m_macrodata[lev].boxArray(),
        m_macrodata[lev].DistributionMap(),
        2 * NQ,
        m_eq[lev].nGrowVect());
    {
        auto const& eq_unit_arrs = eq_unit.arrays();
        amrex::ParallelFor(
            eq_unit, eq_unit.nGrowVect(),
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) == 1) {
                    const auto md_arr = md_arrs[nbx];
                    const auto eq_unit_arr = eq_unit_arrs[nbx];

                    const amrex::Real temperature =
                        md_arr(iv, constants::TEMPERATURE_IDX);
                    const amrex::Real p_by_rho =
                        specific_gas_constant * temperature;

                    // Flow equilibrium: rho=1, local velocity
                    const amrex::RealVect vel = {AMREX_D_DECL(
                        md_arr(iv, constants::VELX_IDX),
                        md_arr(iv, constants::VELY_IDX),
                        md_arr(iv, constants::VELZ_IDX))};
                    const amrex::Real pxx_eq =
                        vel[0] * vel[0] + p_by_rho;
                    const amrex::Real pyy_eq =
                        vel[1] * vel[1] + p_by_rho;
                    const amrex::Real pzz_eq = AMREX_D_PICK(
                        0.0, 0.0, vel[2] * vel[2] + p_by_rho);

                    // Reference equilibrium: rho=1, zero velocity, local T
                    // p_ii^ref = 0^2 + p_by_rho = p_by_rho  (all diagonal components)
                    const amrex::RealVect zero_vel = {AMREX_D_DECL(0.0, 0.0, 0.0)};

                    for (int q = 0; q < NQ; ++q) {
                        eq_unit_arr(iv, q) = set_extended_equilibrium_value(
                            1.0, vel, pxx_eq, pyy_eq, pzz_eq,
                            l_mesh_speed, weight[q], evs[q]);
                        eq_unit_arr(iv, q + NQ) = set_extended_equilibrium_value(
                            1.0, zero_vel, p_by_rho, p_by_rho, p_by_rho,
                            l_mesh_speed, weight[q], evs[q]);
                    }
                }
            });
        amrex::Gpu::synchronize();
    }
    auto const& eq_unit_arrs = eq_unit.const_arrays();

    for (int c = 0; c < m_n_components; ++c) {
        auto const& f_comp_arrs = m_component_lattices[c][lev].arrays();
        amrex::Real diff = m_component_diffusivities[c];

        amrex::ParallelFor(
            m_component_lattices[c][lev], m_eq[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) == 1) {
                    const auto f_comp_arr = f_comp_arrs[nbx];
                    const auto md_arr = md_arrs[nbx];
                    const auto eq_unit_arr = eq_unit_arrs[nbx];

                    amrex::Real rho_comp = 0.0;
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        rho_comp += f_comp_arr(iv, q);
                    }

                    const amrex::Real temperature =
                        md_arr(iv, constants::TEMPERATURE_IDX);
                    const amrex::Real omega_comp =
                        1.0 /
                        (diff / (specific_gas_constant * temperature * dt) +
                         0.5);

                    // Scale cached unit-density shape by rho_comp
                    amrex::GpuArray<amrex::Real, constants::N_MICRO_STATES> eq_all;
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        eq_all[q] = rho_comp * eq_unit_arr(iv, q);
                    }

                    // --- Entropic alpha solve (Ansumali & Karlin, Phys. Rev. E 2002) ---
                    // Find alpha* in (0, 2] s.t. H(f + alpha*(f_eq - f)) = H(f),
                    // where H(f) = sum_q f_q * ln(f_q / eq_ref_q)  and
                    //   eq_ref_q = f^eq(1, 0, T_local)  (zero-velocity reference).
                    // This generalises the standard H = sum f ln(f/w_q), which holds
                    // only at the isothermal point cs^2 = 1/3.
                    //
                    // alpha_use hierarchy (when entropic is enabled):
                    //   Newton succeeds  :  min(omega_comp, alpha*)   [tightest bound]
                    //   Newton skipped   :  min(omega_comp, 1.0)      [safe fallback —
                    //     (negative pops,   monotone approach to eq., never overshoots]
                    //     near-zero rho)
                    //   Entropic disabled:  omega_comp                [pure BGK]
                    //
                    // The safe fallback is critical at low diffusivity where
                    // omega_comp ≈ 2.0: sharp IC edges and refill ghost cells can
                    // develop negative populations that escape the Newton solve and
                    // would otherwise receive the worst-case BGK relaxation rate.
                    // Controlled by lbm.use_entropic_components (default = 0).
                    amrex::Real alpha_use = omega_comp; // entropic disabled: pure BGK
                    if (use_entropic_components) {
                        // Safe fallback for any cell where Newton cannot run.
                        // min(omega_comp, 1.0) guarantees monotone approach
                        // toward equilibrium regardless of population state.
                        alpha_use = amrex::min(omega_comp, 1.0);

                        // Attempt the full entropic solve only when all pre-collision
                        // populations are strictly positive (H0 is well-defined).
                        amrex::Real H0 = 0.0;
                        bool all_positive = true;
                        for (int q = 0; q < NQ; ++q) {
                            amrex::Real fq = f_comp_arr(iv, q);
                            if (fq <= 0.0) { all_positive = false; break; }
                            H0 += fq * log(fq / eq_unit_arr(iv, q + NQ));
                        }

                        if (all_positive) {
                            // Newton iteration: g(alpha) = H(f + alpha*s) - H0 = 0
                            // g'(alpha) = sum_q s_q * (ln(fhat_q / eq_ref_q) + 1)
                            amrex::Real alpha = 2.0; // start at BGK mirror point
                            for (int iter = 0; iter < 10; ++iter) {
                                amrex::Real g = -H0, dg = 0.0;
                                bool fhat_positive = true;
                                for (int q = 0; q < NQ; ++q) {
                                    amrex::Real sq =
                                        eq_all[q] - f_comp_arr(iv, q);
                                    amrex::Real fhat =
                                        f_comp_arr(iv, q) + alpha * sq;
                                    if (fhat <= 0.0) {
                                        fhat_positive = false;
                                        break;
                                    }
                                    amrex::Real ln_fhat_w =
                                        log(fhat / eq_unit_arr(iv, q + NQ));
                                    g  += fhat * ln_fhat_w;
                                    dg += sq * (ln_fhat_w + 1.0);
                                }
                                if (!fhat_positive || fabs(dg) < 1.0e-14) break;
                                alpha -= g / dg;
                                // clamp to [0, 2] for safety
                                alpha = amrex::min(2.0, amrex::max(0.0, alpha));
                                if (fabs(g) < 1.0e-12 * (fabs(H0) + 1.0e-30)) break;
                            }
                            // Tightest bound: Newton result further caps omega_comp
                            alpha_use = amrex::min(omega_comp, alpha);
                        }
                    }

                    // --- Apply collision with entropic alpha ---
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        f_comp_arr(iv, q) +=
                            alpha_use * (eq_all[q] - f_comp_arr(iv, q));
                    }
                }
            });
    }

    amrex::Gpu::synchronize();
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
    }
    m_g[lev].FillBoundary(Geom(lev).periodicity());
}

// calculate the macro fluid properties from the distributions
void LBM::f_to_macrodata(const int lev)
{
    BL_PROFILE("LBM::f_to_macrodata()");
    auto const& md_arrs = m_macrodata[lev].arrays();
    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& f_arrs = m_f[lev].const_arrays();
    auto const& g_arrs = m_g[lev].const_arrays();
    const amrex::Real l_mesh_speed = m_mesh_speed;
    amrex::Real specific_gas_constant = m_R_u / m_m_bar;
    amrex::Real cv = specific_gas_constant / (m_adiabaticExponent - 1.0);

    const bool body_is_isothermal = m_bodyIsIsothermal;
    const amrex::Real body_temperature = m_bodyTemperature;

    const bool body_is_moving = m_body_is_moving;
    const auto body_velocity = m_body_velocity;
    const auto body_angular_velocity = m_body_angular_velocity;
    const auto body_center = m_body_center;
    const amrex::Real current_time = m_ts_new[lev];
    const auto prob_lo = Geom(lev).ProbLoArray();
    const auto dx = Geom(lev).CellSizeArray();

    const bool has_stationary_body = m_has_stationary_body;
    auto const& stat_mask_arrs = m_stationary_mask[lev].const_arrays();

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    amrex::ParallelFor(
        m_macrodata[lev], m_macrodata[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) == 1) {

                const auto f_arr = f_arrs[nbx];
                const auto g_arr = g_arrs[nbx];
                const auto md_arr = md_arrs[nbx];

                amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0;

                amrex::Real pxx(0.0), pyy(0.0), pzz(0.0), pxy(0.0), pxz(0.0),
                    pyz(0.0);
                amrex::Real two_rho_e = 0.0,
                            AMREX_D_DECL(qx = 0.0, qy = 0.0, qz = 0.0);

                for (int q = 0; q < constants::N_MICRO_STATES; q++) {
                    rho += f_arr(iv, q);
                    const auto& ev = evs[q];
                    AMREX_D_DECL(
                        u += ev[0] * f_arr(iv, q), v += ev[1] * f_arr(iv, q),
                        w += ev[2] * f_arr(iv, q));

                    pxx += ev[0] * ev[0] * f_arr(iv, q);
                    pyy += ev[1] * ev[1] * f_arr(iv, q);
                    pxy += ev[0] * ev[1] * f_arr(iv, q);
#if AMREX_SPACEDIM == 3
                    pzz += ev[2] * ev[2] * f_arr(iv, q);
                    pxz += ev[0] * ev[2] * f_arr(iv, q);
                    pyz += ev[1] * ev[2] * f_arr(iv, q);
#endif

                    two_rho_e += g_arr(iv, q);

                    AMREX_D_DECL(
                        qx += ev[0] * g_arr(iv, q), qy += ev[1] * g_arr(iv, q),
                        qz += ev[2] * g_arr(iv, q));
                }
                AMREX_D_DECL(
                    u *= l_mesh_speed / rho, v *= l_mesh_speed / rho,
                    w *= l_mesh_speed / rho);

                // Guard: if rho collapsed to zero (isolated newly-fluid cell with no donor)
                // zero the velocity so we don't produce NaN.  f_to_macrodata will be called
                // again after the next refill so the cell recovers in the next step.
                if (rho <= 0.0) {
                    AMREX_D_DECL(u = 0.0, v = 0.0, w = 0.0);
                }

                if (body_is_moving) {
                    if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_SIDE_IDX) == 1) {
                        bool apply_velocity = true;
                        if (has_stationary_body) {
                            apply_velocity = false;
                            // Check if any neighbor is a moving solid
                            // Moving solid = Solid in is_fluid AND Fluid in stationary_mask
                            for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                                const auto& ev = evs[q];
                                amrex::IntVect iv_nb = iv + ev;
                                if (is_fluid_arrs[nbx](iv_nb, lbm::constants::IS_FLUID_IDX) == 0) {
                                    // It is solid. Is it stationary?
                                    // stationary_mask: 1=Fluid, 0=Solid
                                    if (stat_mask_arrs[nbx](iv_nb) == 1) {
                                        // It is NOT stationary solid, so it must be moving solid
                                        apply_velocity = true;
                                        break;
                                    }
                                }
                            }
                        }

                        if (apply_velocity) {
                            // Calculate body center at current time
                            amrex::Real cx = body_center[0] + body_velocity[0] * current_time;
                            amrex::Real cy = body_center[1] + body_velocity[1] * current_time;
                            amrex::Real cz = body_center[2] + body_velocity[2] * current_time;

                            // Calculate cell center coordinates
                            amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                            amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                            amrex::Real z = 0.0;
#if AMREX_SPACEDIM == 3
                            z = prob_lo[2] + (k + 0.5) * dx[2];
#endif

                            // Calculate velocity due to translation and rotation
                            // v = v_trans + omega x r
                            // r = (x,y,z) - (cx,cy,cz)
                            amrex::Real rx = x - cx;
                            amrex::Real ry = y - cy;
                            amrex::Real rz = z - cz;

                            u = body_velocity[0] + (body_angular_velocity[1] * rz - body_angular_velocity[2] * ry);
                            v = body_velocity[1] + (body_angular_velocity[2] * rx - body_angular_velocity[0] * rz);
#if AMREX_SPACEDIM == 3
                            w = body_velocity[2] + (body_angular_velocity[0] * ry - body_angular_velocity[1] * rx);
#endif
                        }
                    }
                }

                md_arr(iv, constants::RHO_IDX) = rho;
                AMREX_D_DECL(
                    md_arr(iv, constants::VELX_IDX) = u,
                    md_arr(iv, constants::VELY_IDX) = v,
                    md_arr(iv, constants::VELZ_IDX) = w);
                md_arr(iv, constants::VMAG_IDX) =
                    std::sqrt(AMREX_D_TERM(u * u, +v * v, +w * w));

                md_arr(iv, constants::PXX_IDX) = pxx;
                md_arr(iv, constants::PYY_IDX) = pyy;
                md_arr(iv, constants::PZZ_IDX) = pzz;
                md_arr(iv, constants::PXY_IDX) = pxy;
                md_arr(iv, constants::PXZ_IDX) = pxz;
                md_arr(iv, constants::PYZ_IDX) = pyz;

                md_arr(iv, constants::TWO_RHO_E_IDX) = two_rho_e;
                AMREX_D_DECL(
                    md_arr(iv, constants::QX_IDX) = qx,
                    md_arr(iv, constants::QY_IDX) = qy,
                    md_arr(iv, constants::QZ_IDX) = qz);

                amrex::Real temperature =
                    get_temperature(two_rho_e, rho, u, v, w, cv);

                if (body_is_isothermal) {
                    if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_SIDE_IDX) == 1) {
                        temperature = body_temperature;
                    }
                }

                md_arr(iv, constants::TEMPERATURE_IDX) = temperature;

                md_arr(iv, constants::Q_CORR_X_IDX) =
                    rho * u *
                    ((1.0 - 3.0 * specific_gas_constant * temperature) - u * u);
                md_arr(iv, constants::Q_CORR_Y_IDX) =
                    rho * v *
                    ((1.0 - 3.0 * specific_gas_constant * temperature) - v * v);
                md_arr(iv, constants::Q_CORR_Z_IDX) =
                    rho * w *
                    ((1.0 - 3.0 * specific_gas_constant * temperature) - w * w);
            }
        });
    amrex::Gpu::synchronize();
    m_macrodata[lev].FillBoundary(Geom(lev).periodicity());
}

// Compute derived quantities
void LBM::compute_derived(const int lev)
{
    BL_PROFILE("LBM::compute_derived()");
    AMREX_ASSERT(m_macrodata[lev].nGrow() > m_derived[lev].nGrow());
    const auto& idx = geom[lev].InvCellSizeArray();

    // Smagorinsky constant and LB kinematic viscosity captured for epsilon
    const amrex::Real Cs     = constants::SMAGORINSKY_CS;
    const amrex::Real nu_lb  = m_nu;   // LB kinematic viscosity (dx²/step)

    auto const& md_arrs = m_macrodata[lev].const_arrays();
    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& d_arrs = m_derived[lev].arrays();
    const amrex::Box& dbox = geom[lev].Domain();
    amrex::ParallelFor(
        m_derived[lev], m_derived[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
            const auto md_arr = md_arrs[nbx];
            const auto if_arr = is_fluid_arrs[nbx];
            const auto d_arr = d_arrs[nbx];
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));

            if (if_arr(iv, lbm::constants::IS_FLUID_IDX) == 1) {
                // Off-diagonal velocity gradients (for vorticity)
                const amrex::Real vx = gradient(
                    0, constants::VELY_IDX, iv, idx, dbox, if_arr, md_arr);
                const amrex::Real wx = gradient(
                    0, constants::VELZ_IDX, iv, idx, dbox, if_arr, md_arr);
                const amrex::Real uy = gradient(
                    1, constants::VELX_IDX, iv, idx, dbox, if_arr, md_arr);
                const amrex::Real wy = gradient(
                    1, constants::VELZ_IDX, iv, idx, dbox, if_arr, md_arr);
                const amrex::Real uz = AMREX_D_PICK(
                    0, 0,
                    gradient(
                        2, constants::VELX_IDX, iv, idx, dbox, if_arr, md_arr));
                const amrex::Real vz = AMREX_D_PICK(
                    0, 0,
                    gradient(
                        2, constants::VELY_IDX, iv, idx, dbox, if_arr, md_arr));

                d_arr(iv, constants::VORTX_IDX) = wy - vz;
                d_arr(iv, constants::VORTY_IDX) = uz - wx;
                d_arr(iv, constants::VORTZ_IDX) = vx - uy;
                d_arr(iv, constants::VORTM_IDX) = std::sqrt(
                    (wy - vz) * (wy - vz) + (uz - wx) * (uz - wx) +
                    (vx - uy) * (vx - uy));

                // Diagonal velocity gradients for strain rate
                const amrex::Real ux = gradient(
                    0, constants::VELX_IDX, iv, idx, dbox, if_arr, md_arr);
                const amrex::Real vy = gradient(
                    1, constants::VELY_IDX, iv, idx, dbox, if_arr, md_arr);
                const amrex::Real wz = AMREX_D_PICK(
                    0, 0,
                    gradient(
                        2, constants::VELZ_IDX, iv, idx, dbox, if_arr, md_arr));

                // Strain-rate tensor S_ij; S_mag² = 2 * S_ij * S_ij
                const amrex::Real Sxy = 0.5 * (uy + vx);
                const amrex::Real Sxz = 0.5 * (uz + wx);
                const amrex::Real Syz = 0.5 * (vz + wy);
                const amrex::Real S_mag2 = 2.0 * (ux*ux + vy*vy + wz*wz +
                                                   2.0*(Sxy*Sxy + Sxz*Sxz + Syz*Syz));
                const amrex::Real S_mag  = std::sqrt(S_mag2);

                // Smagorinsky SGS viscosity (LB units, dx_LB = 1)
                const amrex::Real nu_sgs = Cs * Cs * S_mag;
                const amrex::Real nu_T   = nu_lb + nu_sgs;

                // Energy dissipation rate ε = ν_T * S_mag² (LB units: dx²/step³)
                d_arr(iv, constants::EPSILON_IDX) = nu_T * S_mag2;
            } else {
                d_arr(iv, constants::EPSILON_IDX) = 0.0;
            }
        });
    amrex::Gpu::synchronize();
}

// Compute derived quantities

void LBM::compute_q_corrections(const int lev)
{
    BL_PROFILE("LBM::compute_derived()");
    AMREX_ASSERT(m_macrodata[lev].nGrow() > m_derived[lev].nGrow());
    const auto& idx = geom[lev].InvCellSizeArray();

    auto const& md_arrs = m_macrodata[lev].const_arrays();
    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& d_arrs = m_derived[lev].arrays();
    const amrex::Box& dbox = geom[lev].Domain();
    amrex::ParallelFor(
        m_derived[lev], m_derived[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
            const auto md_arr = md_arrs[nbx];
            const auto if_arr = is_fluid_arrs[nbx];
            const auto d_arr = d_arrs[nbx];
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));

            if (if_arr(iv, lbm::constants::IS_FLUID_IDX) == 1) {
                d_arr(iv, constants::D_Q_CORR_X_IDX) = gradient(
                    0, constants::Q_CORR_X_IDX, iv, idx, dbox, if_arr, md_arr);
                d_arr(iv, constants::D_Q_CORR_Y_IDX) = gradient(
                    1, constants::Q_CORR_Y_IDX, iv, idx, dbox, if_arr, md_arr);

#if AMREX_SPACEDIM == 3
                d_arr(iv, constants::D_Q_CORR_Z_IDX) = gradient(
                    2, constants::Q_CORR_Z_IDX, iv, idx, dbox, if_arr, md_arr);
#endif
            }
        });
    amrex::Gpu::synchronize();
}

// Compute forces on EB
void LBM::compute_eb_forces()
{
    BL_PROFILE("LBM::compute_eb_forces()");

    amrex::Vector<amrex::Real> forces(AMREX_SPACEDIM, 0);

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    const auto& bounce_dirs = stencil.bounce_dirs;
    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const& f_arrs = m_f[lev].const_arrays();
        auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
        auto const& mask_arrs = m_mask[lev].const_arrays();

        const auto cf = amrex::ParReduce(
            amrex::TypeList<AMREX_D_DECL(
                amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum)>{},
            amrex::TypeList<AMREX_D_DECL(
                amrex::Real, amrex::Real, amrex::Real)>{},
            m_f[lev], amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept
                -> amrex::GpuTuple<AMREX_D_DECL(
                    amrex::Real, amrex::Real, amrex::Real)> {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> fs = {0.0};
                if ((is_fluid_arrs[nbx](iv, lbm::constants::EB_BOUNDARY_IDX) == 1) &&
                    (mask_arrs[nbx](iv) == 0)) {
                    for (int q = 0; q < constants::N_MICRO_STATES; q++) {
                        const auto& ev = evs[q];
                        const amrex::IntVect ivr(iv + evs[bounce_dirs[q]]);

                        for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
                            fs[idir] += 2.0 * ev[idir] * f_arrs[nbx](ivr, q) *
                                        is_fluid_arrs[nbx](ivr, lbm::constants::IS_FLUID_IDX);
                        }
                    }
                }
                return {AMREX_D_DECL(fs[0], fs[1], fs[2])};
            });

        AMREX_D_DECL(
            forces[0] += amrex::get<0>(cf), forces[1] += amrex::get<1>(cf),
            forces[2] += amrex::get<2>(cf));
    }

    amrex::ParallelDescriptor::ReduceRealSum(
        forces.data(), static_cast<int>(forces.size()));

    output_forces_file(forces);
}

// a wrapper for EstTimeStep
void LBM::compute_dt()
{
    BL_PROFILE("LBM::compute_dt()");
    amrex::Vector<amrex::Real> dt_tmp(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = est_time_step(lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(
        dt_tmp.data(), static_cast<int>(dt_tmp.size()));

    constexpr amrex::Real change_max = 1.1;
    amrex::Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max * m_dts[lev]);
        n_factor *= m_nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor * dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const amrex::Real eps = 1.e-3 * dt_0;
    if (m_ts_new[0] + dt_0 > m_stop_time - eps) {
        dt_0 = m_stop_time - m_ts_new[0];
    }

    m_dts[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        m_dts[lev] = m_dts[lev - 1] / m_nsubsteps[lev];
    }
}

// compute dt
amrex::Real LBM::est_time_step(const int /*lev*/)
{
    BL_PROFILE("LBM::est_time_step()");
    return 1.0;
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
void LBM::MakeNewLevelFromCoarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("LBM::MakeNewLevelFromCoarse()");

    m_factory[lev] = amrex::makeEBFabFactory(
        Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);

    m_macrodata[lev].define(
        ba, dm, m_macrodata[lev - 1].nComp(), m_macrodata[lev - 1].nGrow(),
        amrex::MFInfo(), *(m_factory[lev]));
    m_f[lev].define(
        ba, dm, m_f[lev - 1].nComp(), m_f[lev - 1].nGrow(), amrex::MFInfo(),
        *(m_factory[lev]));
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].define(
            ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
            *(m_factory[lev]));
        m_component_lattices[i][lev].setVal(0.0);
    }
    m_g[lev].define(
        ba, dm, m_g[lev - 1].nComp(), m_g[lev - 1].nGrow(), amrex::MFInfo(),
        *(m_factory[lev]));
    m_is_fluid[lev].define(
        ba, dm, m_is_fluid[lev - 1].nComp(), m_is_fluid[lev - 1].nGrow());
    m_is_fluid_fraction[lev].define(ba, dm, 1, m_is_fluid_fraction[lev - 1].nGrow());
    m_eq[lev].define(
        ba, dm, m_eq[lev - 1].nComp(), m_eq[lev - 1].nGrow(), amrex::MFInfo(),
        *(m_factory[lev]));
    m_eq_g[lev].define(
        ba, dm, m_eq_g[lev - 1].nComp(), m_eq_g[lev - 1].nGrow(),
        amrex::MFInfo(), *(m_factory[lev]));
    m_derived[lev].define(
        ba, dm, m_derived[lev - 1].nComp(), m_derived[lev - 1].nGrow(),
        amrex::MFInfo(), *(m_factory[lev]));
    m_mask[lev].define(
        ba, dm, m_mask[lev - 1].nComp(), m_mask[lev - 1].nGrow());
    m_cell_type[lev].define(ba, dm, 1, m_f_nghost);
    m_cell_type[lev].setVal(constants::CELL_LIQUID);
    m_phi_fslbm[lev].define(ba, dm, 1, m_f_nghost);
    m_phi_fslbm[lev].setVal(amrex::Real(1.0));

    m_ts_new[lev] = time;
    m_ts_old[lev] = constants::LOW_NUM;

    initialize_is_fluid(lev);
    // initialize fractional field from integer mask (component 0)
    {
        auto const& if_arrs = m_is_fluid[lev].const_arrays();
        auto const& frac_arrs = m_is_fluid_fraction[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                frac_arrs[nbx](i, j, k, 0) = static_cast<amrex::Real>(
                    if_arrs[nbx](i, j, k, 0));
            });
        amrex::Gpu::synchronize();
    }
    initialize_mask(lev);
    m_fillpatch_op->fillpatch_from_coarse(lev, time, m_f[lev]);
    for (int i = 0; i < m_n_components; ++i) {
        m_component_fillpatch_ops[i]->fillpatch_from_coarse(lev, time, m_component_lattices[i][lev]);
    }

    m_fillpatch_g_op->fillpatch_from_coarse(lev, time, m_g[lev]);

    m_macrodata[lev].setVal(0.0);
    m_eq[lev].setVal(0.0);
    m_eq_g[lev].setVal(0.0);
    m_derived[lev].setVal(0.0);

    f_to_macrodata(lev);

    compute_q_corrections(lev);

    macrodata_to_equilibrium(lev);

    compute_derived(lev);
}

// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
void LBM::MakeNewLevelFromScratch(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("LBM::MakeNewLevelFromScratch()");

    m_factory[lev] = amrex::makeEBFabFactory(
        Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);

    m_macrodata[lev].define(
        ba, dm, constants::N_MACRO_STATES, m_macrodata_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_f[lev].define(
        ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].define(
            ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
            *(m_factory[lev]));
        m_component_lattices[i][lev].setVal(0.0);
    }
    m_g[lev].define(
        ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_is_fluid[lev].define(ba, dm, constants::N_IS_FLUID, m_f[lev].nGrow());
    m_is_fluid_fraction[lev].define(ba, dm, 1, m_is_fluid[lev].nGrow());
    m_eq[lev].define(
        ba, dm, constants::N_MICRO_STATES, m_eq_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_eq_g[lev].define(
        ba, dm, constants::N_MICRO_STATES, m_eq_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_derived[lev].define(
        ba, dm, constants::N_DERIVED, m_derived_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_mask[lev].define(ba, dm, 1, 0);
    m_stationary_mask[lev].define(ba, dm, 1, m_is_fluid[lev].nGrow());
    m_cell_type[lev].define(ba, dm, 1, m_f_nghost);
    m_cell_type[lev].setVal(constants::CELL_LIQUID);
    m_phi_fslbm[lev].define(ba, dm, 1, m_f_nghost);
    m_phi_fslbm[lev].setVal(amrex::Real(1.0));

    m_ts_new[lev] = time;
    m_ts_old[lev] = constants::LOW_NUM;

    // Initialize the data
    init_stationary_body(lev);
    initialize_is_fluid(lev);
    // initialize fractional field from integer mask (component 0)
    {
        auto const& if_arrs = m_is_fluid[lev].const_arrays();
        auto const& frac_arrs = m_is_fluid_fraction[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                frac_arrs[nbx](i, j, k, 0) = static_cast<amrex::Real>(
                    if_arrs[nbx](i, j, k, 0));
            });
        amrex::Gpu::synchronize();
    }

    // FSLBM (Körner 2005): sharp-interface cell-type + fill-level initialization.
    if (m_free_surface) {
        fslbm_init_cell_type(lev);
    }

    initialize_mask(lev);
    initialize_f(lev);
    m_macrodata[lev].setVal(0.0);
    m_eq[lev].setVal(0.0);
    m_eq_g[lev].setVal(0.0);
    m_derived[lev].setVal(0.0);

    f_to_macrodata(lev);

    compute_q_corrections(lev);

    macrodata_to_equilibrium(lev);

    compute_derived(lev);
}

void LBM::initialize_f(const int lev)
{
    BL_PROFILE("LBM::initialize_f()");

    m_ic_op->initialize(lev, geom[lev].data());
    for (int i = 0; i < m_n_components; ++i) {
        m_component_ic_ops[i]->initialize_lattice(lev, geom[lev].data(), m_component_lattices[i][lev]);
    }

    fill_f_inside_eb(lev);

    // Zero out inside EB for additional components
    auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
    for (int c = 0; c < m_n_components; ++c) {
        auto const& f_arrs = m_component_lattices[c][lev].arrays();
        amrex::ParallelFor(
            m_component_lattices[c][lev], m_component_lattices[c][lev].nGrowVect(), constants::N_MICRO_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
                if (is_fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) == 0) {
                    f_arrs[nbx](i, j, k, q) = 0.0;
                }
            });
    }
    amrex::Gpu::synchronize();

    m_f[lev].FillBoundary(Geom(lev).periodicity());
    m_g[lev].FillBoundary(Geom(lev).periodicity());
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
    }
}

void LBM::initialize_moving_body_shape(int lev)
{
    if (m_using_voxel_body) return;

    amrex::ParmParse pp("eb2");
    std::string geom_type;
    pp.query("geom_type", geom_type);
    
    std::string stl_file;
    pp.query("stl_file", stl_file);
    
    amrex::ParmParse ppvc("voxel_cracks");
    std::string vc_file;
    ppvc.query("crack_file", vc_file);
    int use_voxel_cracks = 0;
    pp.query("use_voxel_cracks", use_voxel_cracks);

    bool is_file_based = (geom_type == "stl") || 
                         (!stl_file.empty()) || 
                         (use_voxel_cracks != 0) || 
                         (!vc_file.empty());

    if (is_file_based) {
        amrex::Print() << "Initializing moving body reference from current fluid field..." << std::endl;
        
        const amrex::Geometry& geom = Geom(lev);
        const amrex::Box& domain = geom.Domain();
        const int nx = domain.length(0);
        const int ny = domain.length(1);
        const int nz = domain.length(2);
        size_t num_cells = static_cast<size_t>(nx) * ny * nz;
        
        // Gather m_is_fluid[lev]
        amrex::BoxArray ba_full(domain);
        amrex::Vector<int> pmap(1, amrex::ParallelDescriptor::MyProc());
        amrex::DistributionMapping dm_local(pmap);
        amrex::iMultiFab local_imf(ba_full, dm_local, 1, 0);
        
        local_imf.ParallelCopy(m_is_fluid[lev]);
        
        m_body_voxel_data.resize(num_cells);
        auto* voxel_ptr = m_body_voxel_data.data();
        
        for (amrex::MFIter mfi(local_imf); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            auto const& fab_arr = local_imf.array(mfi);
            
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                size_t idx = k * (nx * ny) + j * nx + i;
                voxel_ptr[idx] = static_cast<uint16_t>(fab_arr(i, j, k, lbm::constants::IS_FLUID_IDX));
            });
        }
        amrex::Gpu::synchronize();
        
        m_using_voxel_body = true;
        
        // Set metadata
        m_body_voxel_dims = domain.length();
        m_body_voxel_origin = geom.ProbLoArray();
        m_body_voxel_dx = geom.CellSizeArray();
        
        amrex::Vector<amrex::Real> center(3, 0.0);
        pp.queryarr("stl_center", center);
        
        if (center[0] == 0.0 && center[1] == 0.0 && center[2] == 0.0) {
             m_body_initial_center[0] = m_body_center[0];
             m_body_initial_center[1] = m_body_center[1];
             m_body_initial_center[2] = m_body_center[2];
        } else {
             m_body_initial_center[0] = center[0];
             m_body_initial_center[1] = center[1];
             m_body_initial_center[2] = center[2];
        }
    }
}

void LBM::init_stationary_body(int lev)
{
    BL_PROFILE("LBM::init_stationary_body()");
    
    amrex::ParmParse pp("eb2");
    std::string stl_file;
    std::string crack_file;
    
    m_stationary_mask[lev].setVal(1); // Default to Fluid (1)
    
    bool has_stl = pp.query("stationary_stl_file", stl_file);
    bool has_crack = pp.query("stationary_crack_file", crack_file);
    
    if (has_stl) {
        m_has_stationary_body = true;
        amrex::Print() << "Loading stationary STL: " << stl_file << std::endl;
        
        amrex::Real scale = 1.0;
        int reverse_normal = 0;
        amrex::Array<amrex::Real, 3> center = {0.0, 0.0, 0.0};
        pp.query("stationary_stl_scale", scale);
        pp.query("stationary_stl_reverse_normal", reverse_normal);
        pp.query("stationary_stl_center", center);

        amrex::STLtools stlobj;
        stlobj.read_stl_file(stl_file, scale, center, reverse_normal);

        amrex::MultiFab marker(
            m_stationary_mask[lev].boxArray(), m_stationary_mask[lev].DistributionMap(), 1,
            m_stationary_mask[lev].nGrow());

        const amrex::Real outside_value = 1.0; // Fluid
        const amrex::Real inside_value = 0.0;  // Solid
        marker.setVal(1.0);
        stlobj.fill(
            marker, marker.nGrowVect(), Geom(lev), outside_value, inside_value);
        amrex::Gpu::synchronize();

        auto const& marker_arrs = marker.const_arrays();
        auto const& mask_arrs = m_stationary_mask[lev].arrays();
        amrex::ParallelFor(
            m_stationary_mask[lev], m_stationary_mask[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Combine with existing mask (intersection of fluids -> min)
                // 0=Solid, 1=Fluid. min(1, 0) = 0 (Solid).
                int val = static_cast<int>(marker_arrs[nbx](i, j, k, 0));
                mask_arrs[nbx](i, j, k) = amrex::min(mask_arrs[nbx](i, j, k), val);
            });
        amrex::Gpu::synchronize();
    }
    
    if (has_crack) {
        m_has_stationary_body = true;
        const auto& geom = Geom(lev);
        const amrex::Box& domain = geom.Domain();
        const int nx = domain.length(0);
        const int ny = domain.length(1);
        const int nz = domain.length(2);
        
        amrex::Print() << "Loading stationary crack file: " << crack_file << std::endl;
        
        std::vector<uint16_t> crack_data = read_crack_file(crack_file, nx, ny, nz);
        
        amrex::Gpu::DeviceVector<uint16_t> d_crack_data(crack_data.size());
        amrex::Gpu::copyAsync(
            amrex::Gpu::hostToDevice, crack_data.begin(), crack_data.end(),
            d_crack_data.begin());
        amrex::Gpu::synchronize();

        auto const* crack_ptr = d_crack_data.data();
        
        for (amrex::MFIter mfi(m_stationary_mask[lev]); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            auto const& mask_arr = m_stationary_mask[lev].array(mfi);

            amrex::ParallelFor(
                box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    int file_index = k * (nx * ny) + j * nx + i;
                    // File: 0=Fluid, 1=Solid
                    // Mask: 1=Fluid, 0=Solid
                    int val = (crack_ptr[file_index] == 0) ? 1 : 0;
                    mask_arr(i, j, k) = amrex::min(mask_arr(i, j, k), val);
                });
        }
        amrex::Gpu::synchronize();
    }
    
    // Also check for stationary parser function (handled in reconstruct_body_sdf, but we set flag here)
    std::string stationary_parser_function;
    if (pp.query("stationary_parser_function", stationary_parser_function)) {
        m_has_stationary_body = true;
    }
}

void LBM::initialize_is_fluid(const int lev)
{
    BL_PROFILE("LBM::initialize_is_fluid()");
    const auto* factory =
        static_cast<amrex::EBFArrayBoxFactory*>(m_factory[lev].get());
    auto const& flags = factory->getMultiEBCellFlagFab();
    auto const& flag_arrs = flags.const_arrays();
    m_is_fluid[lev].setVal(0.0);
    auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
    amrex::ParallelFor(
        m_is_fluid[lev], m_is_fluid[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            is_fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) =
                !(flag_arrs[nbx](i, j, k).isRegular() ||
                  flag_arrs[nbx](i, j, k).isSingleValued())
                    ? 0
                    : 1;
                        // ensure new 4th component is initialized to 0
                        //is_fluid_arrs[nbx](i, j, k, 3) = 0;
        });

    initialize_from_stl(Geom(lev), m_is_fluid[lev]);

    // If body is moving, reconstruct the SDF at t=0 to ensure correct initial position
    if (m_body_is_moving) {
        initialize_moving_body_shape(lev);
        reconstruct_body_sdf(lev, 0.0);
        // Update is_fluid from the reconstructed fraction
        update_is_fluid_from_fraction_and_mark(lev, m_is_fluid_fraction_threshold);
    }

    if (m_has_stationary_body) {
        auto const& stationary_mask_arrs = m_stationary_mask[lev].const_arrays();
        auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Merge stationary mask: if stationary mask is 0 (solid), is_fluid becomes 0
                // is_fluid = min(is_fluid, stationary_mask)
                is_fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) = amrex::min(
                    is_fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX),
                    stationary_mask_arrs[nbx](i, j, k)
                );
            });
        amrex::Gpu::synchronize();
    }

    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());

    // Compute the boundary cells
    amrex::ParallelFor(
        m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const auto if_arr = is_fluid_arrs[nbx];

            bool all_covered = true;
            const amrex::IntVect nn(1);
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
                const auto dimvec = amrex::IntVect::TheDimensionVector(idir);
                for (int n = 1; n <= nn[idir]; n++) {
                    all_covered &= (if_arr(iv - n * dimvec, lbm::constants::IS_FLUID_IDX) == 0) &&
                                   (if_arr(iv + n * dimvec, lbm::constants::IS_FLUID_IDX) == 0);
                }
            }

            if ((all_covered) || (if_arr(iv, lbm::constants::IS_FLUID_IDX) == 1)) {
                if_arr(iv, lbm::constants::EB_BOUNDARY_IDX) = 0;
            } else {
                if_arr(iv, lbm::constants::EB_BOUNDARY_IDX) = 1;
            }
        });

    // Compute the boundary cells on the fluid side
    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    amrex::ParallelFor(
        m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const auto if_arr = is_fluid_arrs[nbx];

            bool all_covered = true;
            for (int idir = 0; idir < constants::N_MICRO_STATES; idir++) {
                const auto& dimvec = evs[idir];
                all_covered &= (if_arr(iv - dimvec, lbm::constants::IS_FLUID_IDX) == 1);
            }

            if ((all_covered) || (if_arr(iv, lbm::constants::IS_FLUID_IDX) == 0)) {
                if_arr(iv, lbm::constants::IS_FLUID_SIDE_IDX) = 0;
            } else {
                if_arr(iv, lbm::constants::IS_FLUID_SIDE_IDX) = 1;
            }
        });

    // Compute the boundary cells of the fluid side boundary
    amrex::ParallelFor(
        m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            const auto if_arr = is_fluid_arrs[nbx];

            // mark cells that are fluid, not already marked as side boundary
            // (component 2), but that see at least one neighbor with comp 2
            constexpr int IS_FLUID = lbm::constants::IS_FLUID_IDX;
            constexpr int IS_FLUID_SIDE = lbm::constants::IS_FLUID_SIDE_IDX;
            constexpr int IS_FLUID_SIDE_BOUNDARY =
                lbm::constants::IS_FLUID_SIDE_BOUNDARY_IDX;

            bool sees_side = false;
            for (int idir = 0; idir < constants::N_MICRO_STATES; ++idir) {
                const auto& dimvec = evs[idir];
                if (if_arr(iv - dimvec, IS_FLUID_SIDE) == 1) {
                    sees_side = true;
                    break;
                }
            }

            if ((if_arr(iv, IS_FLUID) == 1) && (if_arr(iv, IS_FLUID_SIDE) == 0) &&
                sees_side) {
                if_arr(iv, IS_FLUID_SIDE_BOUNDARY) = 1;
            } else {
                if_arr(iv, IS_FLUID_SIDE_BOUNDARY) = 0;
            }
        });
  

    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());
}

void LBM::update_is_fluid_from_fraction_and_mark(const int lev, amrex::Real threshold)
{
    BL_PROFILE("LBM::update_is_fluid_from_fraction_and_mark()");

    if (threshold < 0.0) threshold = m_is_fluid_fraction_threshold;

    // Step 1: threshold fractional field into integer mask component 0
    {
        auto const& frac_arrs = m_is_fluid_fraction[lev].const_arrays();
        auto const& isf_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real val = frac_arrs[nbx](i, j, k, 0);
                isf_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) =
                    (val >= threshold) ? 1 : 0;
            });
        amrex::Gpu::synchronize();
    }

    // After modifying the integer mask, recompute the boundary markers
    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());

    // Compute EB_BOUNDARY similar to initialize_is_fluid
    {
        auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto if_arr = is_fluid_arrs[nbx];

                bool all_covered = true;
                const amrex::IntVect nn(1);
                for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
                    const auto dimvec = amrex::IntVect::TheDimensionVector(idir);
                    for (int n = 1; n <= nn[idir]; n++) {
                        all_covered &= (if_arr(iv - n * dimvec, lbm::constants::IS_FLUID_IDX) == 0) &&
                                       (if_arr(iv + n * dimvec, lbm::constants::IS_FLUID_IDX) == 0);
                    }
                }

                if ((all_covered) || (if_arr(iv, lbm::constants::IS_FLUID_IDX) == 1)) {
                    if_arr(iv, lbm::constants::EB_BOUNDARY_IDX) = 0;
                } else {
                    if_arr(iv, lbm::constants::EB_BOUNDARY_IDX) = 1;
                }
            });
        amrex::Gpu::synchronize();
    }

    // Compute IS_FLUID_SIDE
    {
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto if_arr = is_fluid_arrs[nbx];

                bool all_covered = true;
                for (int idir = 0; idir < constants::N_MICRO_STATES; idir++) {
                    const auto& dimvec = evs[idir];
                    all_covered &= (if_arr(iv - dimvec, lbm::constants::IS_FLUID_IDX) == 1);
                }

                if ((all_covered) || (if_arr(iv, lbm::constants::IS_FLUID_IDX) == 0)) {
                    if_arr(iv, lbm::constants::IS_FLUID_SIDE_IDX) = 0;
                } else {
                    if_arr(iv, lbm::constants::IS_FLUID_SIDE_IDX) = 1;
                }
            });
        amrex::Gpu::synchronize();
    }

    // Compute IS_FLUID_SIDE_BOUNDARY
    {
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto if_arr = is_fluid_arrs[nbx];

                constexpr int IS_FLUID = lbm::constants::IS_FLUID_IDX;
                constexpr int IS_FLUID_SIDE = lbm::constants::IS_FLUID_SIDE_IDX;
                constexpr int IS_FLUID_SIDE_BOUNDARY =
                    lbm::constants::IS_FLUID_SIDE_BOUNDARY_IDX;

                bool sees_side = false;
                for (int idir = 0; idir < constants::N_MICRO_STATES; ++idir) {
                    const auto& dimvec = evs[idir];
                    if (if_arr(iv - dimvec, IS_FLUID_SIDE) == 1) {
                        sees_side = true;
                        break;
                    }
                }

                if ((if_arr(iv, IS_FLUID) == 1) && (if_arr(iv, IS_FLUID_SIDE) == 0) &&
                    sees_side) {
                    if_arr(iv, IS_FLUID_SIDE_BOUNDARY) = 1;
                } else {
                    if_arr(iv, IS_FLUID_SIDE_BOUNDARY) = 0;
                }
            });
        amrex::Gpu::synchronize();
    }

    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());
}

void LBM::refill_and_spill(const int lev, amrex::Real threshold)
{
    BL_PROFILE("LBM::refill_and_spill()");

    if (threshold < 0.0) threshold = m_is_fluid_fraction_threshold;

    // Step 1: Fill boundary cells for all data we'll need
    m_is_fluid_fraction[lev].FillBoundary(Geom(lev).periodicity());
    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());
    m_cell_type[lev].FillBoundary(Geom(lev).periodicity());
    m_phi_fslbm[lev].FillBoundary(Geom(lev).periodicity());
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
    }
    m_g[lev].FillBoundary(Geom(lev).periodicity());

    // Step 2: Save old fluid mask AND boundary layers BEFORE updating
    amrex::iMultiFab old_is_fluid(
        m_is_fluid[lev].boxArray(), m_is_fluid[lev].DistributionMap(), 1, 1);
    amrex::iMultiFab old_fluid_side_boundary(
        m_is_fluid[lev].boxArray(), m_is_fluid[lev].DistributionMap(), 1, 1);
    
    amrex::iMultiFab::Copy(old_is_fluid, m_is_fluid[lev], 
                          lbm::constants::IS_FLUID_IDX, 0, 1, 0);
    amrex::iMultiFab::Copy(old_fluid_side_boundary, m_is_fluid[lev], 
                          lbm::constants::IS_FLUID_SIDE_BOUNDARY_IDX, 0, 1, 0);
    old_is_fluid.FillBoundary(Geom(lev).periodicity());
    old_fluid_side_boundary.FillBoundary(Geom(lev).periodicity());

    // Step 3: Update fluid mask based on new fractional values
    update_is_fluid_from_fraction_and_mark(lev, threshold);

    // Step 4a: Identify cells that changed state
    amrex::iMultiFab newly_fluid(
        m_is_fluid[lev].boxArray(), m_is_fluid[lev].DistributionMap(), 1, 0);
    amrex::iMultiFab newly_solid(
        m_is_fluid[lev].boxArray(), m_is_fluid[lev].DistributionMap(), 1, 0);

    {
        auto const& old_arrs = old_is_fluid.const_arrays();
        auto const& new_arrs = m_is_fluid[lev].const_arrays();
        auto const& newly_fluid_arrs = newly_fluid.arrays();
        auto const& newly_solid_arrs = newly_solid.arrays();
        
        amrex::ParallelFor(newly_fluid, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            int old_val = old_arrs[nbx](i, j, k, 0);
            int new_val = new_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX);
            newly_fluid_arrs[nbx](i, j, k, 0) = (old_val == 0 && new_val == 1) ? 1 : 0;
            newly_solid_arrs[nbx](i, j, k, 0) = (old_val == 1 && new_val == 0) ? 1 : 0;
        });
    }
    amrex::Gpu::synchronize();

    // Step 4b: SPILL - Distribute mass/energy from newly solid cells to OLD outer boundary layer
    // Use stencil weights (proportional to velocity) for distribution
    
    // Create temporary MultiFabs to accumulate spilled mass (to handle ghost cell updates correctly)
    amrex::MultiFab spill_f(m_f[lev].boxArray(), m_f[lev].DistributionMap(), constants::N_MICRO_STATES, m_f[lev].nGrow());
    amrex::MultiFab spill_g(m_g[lev].boxArray(), m_g[lev].DistributionMap(), constants::N_MICRO_STATES, m_g[lev].nGrow());
    spill_f.setVal(0.0);
    spill_g.setVal(0.0);

    {
        auto const& newly_solid_arrs = newly_solid.const_arrays();
        auto const& frac_arrs = m_is_fluid_fraction[lev].const_arrays();
        auto const& old_boundary_arrs = old_fluid_side_boundary.const_arrays();
        auto const& curr_fluid_arrs = m_is_fluid[lev].const_arrays();
        auto const& f_arrs = m_f[lev].arrays();
        auto const& g_arrs = m_g[lev].arrays();
        auto const& spill_f_arrs = spill_f.arrays();
        auto const& spill_g_arrs = spill_g.arrays();
        // FSLBM: don't spill cells that FSLBM classifies as free-surface/liquid/gas.
        // These cells' IS_FLUID can oscillate each step near solid walls due to the
        // tanh-smoothed SDF, but their f distributions are managed by FSLBM, not
        // by the body-motion refill/spill.
        auto const& ct_arrs_sp = m_cell_type[lev].const_arrays();
        
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        const auto& weights = stencil.weights;
        
        amrex::ParallelFor(m_f[lev], amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Only process cells that became solid AND are covered by EB
                if (newly_solid_arrs[nbx](i, j, k, 0) != 1) return;
                
                // Skip cells whose f is zero by definition and have nothing to spill:
                //   CELL_GAS — above the free surface; f is always zero.
                // CELL_INTERFACE and CELL_LIQUID cells must be spilled: if the
                // impeller blade sweeps into a free-surface or bulk-liquid cell,
                // the f content must be redistributed to fluid neighbors instead
                // of being silently discarded.
                const int ct_cell = ct_arrs_sp[nbx](i, j, k, 0);
                if (ct_cell == lbm::constants::CELL_GAS) { return; }
                
                amrex::Real frac = frac_arrs[nbx](i, j, k, 0);
                if (frac >= 1.0) return; // Not an EB cell
                
                // Get bounds for safety
                const auto& farr = f_arrs[nbx];
                const auto lo = amrex::lbound(farr);
                const auto hi = amrex::ubound(farr);
                
                // First pass: compute sum of weights for OLD IS_FLUID_SIDE_BOUNDARY neighbors
                amrex::Real weight_sum = 0.0;
                for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                    int ni = i + evs[nq][0];
                    int nj = j + evs[nq][1];
                    int nk = k + evs[nq][2];
                    
                    // Check bounds
                    if (ni < lo.x || ni > hi.x || 
                        nj < lo.y || nj > hi.y || 
                        nk < lo.z || nk > hi.z) continue;
                    
                    // Add weight if neighbor was on the OLD outer boundary layer
                    if (old_boundary_arrs[nbx](ni, nj, nk, 0) == 1) {
                        weight_sum += weights[nq];
                    }
                }
                
                // Fallback: if no old-boundary neighbor, widen to any currently-fluid neighbor
                bool use_fallback = (weight_sum == 0.0);
                if (use_fallback) {
                    for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                        int ni = i + evs[nq][0];
                        int nj = j + evs[nq][1];
                        int nk = k + evs[nq][2];
                        if (ni < lo.x || ni > hi.x ||
                            nj < lo.y || nj > hi.y ||
                            nk < lo.z || nk > hi.z) continue;
                        if (curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                            weight_sum += weights[nq];
                        }
                    }
                }

                // If still no fluid neighbor (cell fully buried in solid), mass is truly lost
                if (weight_sum == 0.0) {
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        f_arrs[nbx](i, j, k, q) = 0.0;
                        g_arrs[nbx](i, j, k, q) = 0.0;
                    }
                    return;
                }
                
                // Second pass: distribute to neighbors using normalized weights
                for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                    int ni = i + evs[nq][0];
                    int nj = j + evs[nq][1];
                    int nk = k + evs[nq][2];
                    
                    // Check bounds
                    if (ni < lo.x || ni > hi.x || 
                        nj < lo.y || nj > hi.y || 
                        nk < lo.z || nk > hi.z) continue;
                    
                    // Primary: old boundary neighbor; fallback: any fluid neighbor
                    bool valid = use_fallback
                        ? (curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1)
                        : (old_boundary_arrs[nbx](ni, nj, nk, 0) == 1);

                    if (valid) {
                        amrex::Real w = weights[nq] / weight_sum;
                        for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                            amrex::Gpu::Atomic::AddNoRet(&spill_f_arrs[nbx](ni, nj, nk, q), f_arrs[nbx](i, j, k, q) * w);
                            amrex::Gpu::Atomic::AddNoRet(&spill_g_arrs[nbx](ni, nj, nk, q), g_arrs[nbx](i, j, k, q) * w);
                        }
                    }
                }
                
                // Zero out the newly solid cell after distribution
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    f_arrs[nbx](i, j, k, q) = 0.0;
                    g_arrs[nbx](i, j, k, q) = 0.0;
                }
            });
    }
    amrex::Gpu::synchronize();

    // Sum spilled mass from ghost cells to valid cells
    spill_f.SumBoundary(Geom(lev).periodicity());
    spill_g.SumBoundary(Geom(lev).periodicity());

    // Add spilled mass to the main fluid arrays
    amrex::MultiFab::Add(m_f[lev], spill_f, 0, 0, constants::N_MICRO_STATES, 0);
    amrex::MultiFab::Add(m_g[lev], spill_g, 0, 0, constants::N_MICRO_STATES, 0);

    // Sync boundaries so everyone sees the updated mass
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
    }
    m_g[lev].FillBoundary(Geom(lev).periodicity());

    // Spill for components
    for (int c = 0; c < m_n_components; ++c) {
        spill_f.setVal(0.0); // Reuse spill_f buffer

        auto const& newly_solid_arrs = newly_solid.const_arrays();
        auto const& frac_arrs = m_is_fluid_fraction[lev].const_arrays();
        auto const& old_boundary_arrs = old_fluid_side_boundary.const_arrays();
        auto const& curr_fluid_arrs = m_is_fluid[lev].const_arrays();
        auto const& f_comp_arrs = m_component_lattices[c][lev].arrays();
        auto const& spill_comp_arrs = spill_f.arrays();

        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        const auto& weights = stencil.weights;

        amrex::ParallelFor(
            m_component_lattices[c][lev], amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Only process cells that became solid AND are covered by EB
                if (newly_solid_arrs[nbx](i, j, k, 0) != 1) return;

                amrex::Real frac = frac_arrs[nbx](i, j, k, 0);
                if (frac >= 1.0) return; // Not an EB cell

                // Get bounds for safety
                const auto& farr = f_comp_arrs[nbx];
                const auto lo = amrex::lbound(farr);
                const auto hi = amrex::ubound(farr);

                // First pass: compute sum of weights for OLD
                // IS_FLUID_SIDE_BOUNDARY neighbors that are STILL FLUID
                amrex::Real weight_sum = 0.0;
                for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                    int ni = i + evs[nq][0];
                    int nj = j + evs[nq][1];
                    int nk = k + evs[nq][2];

                    // Check bounds
                    if (ni < lo.x || ni > hi.x || nj < lo.y || nj > hi.y ||
                        nk < lo.z || nk > hi.z)
                        continue;

                    // Add weight if neighbor was on the OLD outer boundary
                    // layer AND is still fluid in the current state
                    if (old_boundary_arrs[nbx](ni, nj, nk, 0) == 1 &&
                        curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                        weight_sum += weights[nq];
                    }
                }

                // Fallback: if no old-boundary+still-fluid neighbor found,
                // widen search to ANY currently-fluid neighbor (captures
                // blade leading-edge / corner cells that lose all their
                // old boundary neighbors in the same step).
                bool use_fallback = (weight_sum == 0.0);
                if (use_fallback) {
                    for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                        int ni = i + evs[nq][0];
                        int nj = j + evs[nq][1];
                        int nk = k + evs[nq][2];
                        if (ni < lo.x || ni > hi.x || nj < lo.y || nj > hi.y ||
                            nk < lo.z || nk > hi.z) continue;
                        if (curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                            weight_sum += weights[nq];
                        }
                    }
                }

                // If still no fluid neighbor (cell fully buried in solid), mass is lost
                if (weight_sum == 0.0) {
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        f_comp_arrs[nbx](i, j, k, q) = 0.0;
                    }
                    return;
                }

                // Second pass: distribute to neighbors using normalized weights
                for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                    int ni = i + evs[nq][0];
                    int nj = j + evs[nq][1];
                    int nk = k + evs[nq][2];

                    // Check bounds
                    if (ni < lo.x || ni > hi.x || nj < lo.y || nj > hi.y ||
                        nk < lo.z || nk > hi.z)
                        continue;

                    // Primary: old boundary + still fluid; fallback: any fluid
                    bool valid = use_fallback
                        ? (curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1)
                        : (old_boundary_arrs[nbx](ni, nj, nk, 0) == 1 &&
                           curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1);

                    if (valid) {
                        amrex::Real w = weights[nq] / weight_sum;
                        for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                            amrex::Gpu::Atomic::AddNoRet(
                                &spill_comp_arrs[nbx](ni, nj, nk, q),
                                f_comp_arrs[nbx](i, j, k, q) * w);
                        }
                    }
                }

                // Zero out the newly solid cell after distribution
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    f_comp_arrs[nbx](i, j, k, q) = 0.0;
                }
            });
        amrex::Gpu::synchronize();

        spill_f.SumBoundary(Geom(lev).periodicity());
        amrex::MultiFab::Add(
            m_component_lattices[c][lev], spill_f, 0, 0,
            constants::N_MICRO_STATES, 0);
        m_component_lattices[c][lev].FillBoundary(Geom(lev).periodicity());
    }

    // Check if there are any newly fluid cells - if not, skip refill
    amrex::Long num_newly_fluid = newly_fluid.sum(0);
    if (num_newly_fluid == 0) {
        // No cells transitioned to fluid - nothing to refill
        m_f[lev].FillBoundary(Geom(lev).periodicity());
        for (int i = 0; i < m_n_components; ++i) {
            m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
        }
        m_g[lev].FillBoundary(Geom(lev).periodicity());
        m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());
        return;
    }

    // Step 5: Refill newly fluid cells by finding donor in the normal direction
    // Normal is computed from averaged evs of persistent fluid neighbors
    
    amrex::iMultiFab donor_recipient_count(
        m_is_fluid[lev].boxArray(), m_is_fluid[lev].DistributionMap(), 1, 1);
    donor_recipient_count.setVal(0);

    // Common arrays/flags shared by Step 5 and Step 6 refill lambdas
    auto const& newly_fluid_arrs_outer = newly_fluid.const_arrays();
    auto const& f_arrs_outer = m_f[lev].arrays();
    auto const& g_arrs_outer = m_g[lev].arrays();
    auto const& old_fluid_arrs_outer = old_is_fluid.const_arrays();
    auto const& curr_fluid_arrs_outer = m_is_fluid[lev].const_arrays();
    auto const& ct_arrs_refill = m_cell_type[lev].const_arrays();
    const bool is_free_surface_refill = m_free_surface;
    // phi arrays used to classify body-vacated CELL_SOLID cells: average neighbor
    // phi > 0.5 → cell will become CELL_LIQUID (let refill seed it from neighbors);
    // average phi <= 0.5 → cell will become CELL_GAS (Case B zeros f/g, skip here).
    auto const& phi_fslbm_arrs_ref = m_phi_fslbm[lev].const_arrays();
    auto const& donor_count_arrs_outer = donor_recipient_count.arrays();

    {
        auto const& newly_fluid_arrs = newly_fluid_arrs_outer;
        auto const& f_arrs = f_arrs_outer;
        auto const& g_arrs = g_arrs_outer;
        auto const& old_fluid_arrs = old_fluid_arrs_outer;
        auto const& curr_fluid_arrs = curr_fluid_arrs_outer;
        auto const& donor_count_arrs = donor_recipient_count.arrays();
        auto const& ct_arrs_ref = ct_arrs_refill;
        auto const& phi_arrs_ref = phi_fslbm_arrs_ref;
        const bool is_free_surface = is_free_surface_refill;
        
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;

        amrex::ParallelFor(m_f[lev], amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            
            // Only process newly fluid cells
            if (newly_fluid_arrs[nbx](i, j, k, 0) != 1) return;

            // Skip cells managed exclusively by FSLBM:
            //   CELL_GAS       — f=0 by definition; no liquid neighbors to donate from.
            //   CELL_INTERFACE — managed by FSLBM ABB/mass-flux.
            //   CELL_SOLID     — body just vacated this cell. Use neighbor-phi average
            //     to predict the outcome of Case B in fslbm_advance_surface:
            //       avg_phi > 0.5  → will become CELL_LIQUID: fall through to
            //                        neighbor-interpolation refill.
            //       avg_phi ≤ 0.5  → will become CELL_GAS: skip (Case B zeros f/g).
            // CELL_LIQUID cells (bulk liquid adjacent to impeller) always refilled.
            if (is_free_surface) {
                const int ct_val = ct_arrs_ref[nbx](i, j, k, 0);
                if (ct_val == lbm::constants::CELL_GAS ||
                    ct_val == lbm::constants::CELL_INTERFACE) { return; }
                if (ct_val == lbm::constants::CELL_SOLID) {
                    // Average phi of non-solid neighbors to predict Case B outcome.
                    const auto& phi_arr = phi_arrs_ref[nbx];
                    const auto lo = amrex::lbound(phi_arr);
                    const auto hi = amrex::ubound(phi_arr);
                    amrex::Real phi_sum = amrex::Real(0.0);
                    int phi_cnt = 0;
                    for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                        const int ni = i + evs[nq][0];
                        const int nj = j + evs[nq][1];
                        const int nk = k + evs[nq][2];
                        if (ni < lo.x || ni > hi.x ||
                            nj < lo.y || nj > hi.y ||
                            nk < lo.z || nk > hi.z) { continue; }
                        if (ct_arrs_ref[nbx](ni, nj, nk, 0) == lbm::constants::CELL_SOLID) { continue; }
                        phi_sum += phi_arr(ni, nj, nk, 0);
                        ++phi_cnt;
                    }
                    const amrex::Real avg_phi = (phi_cnt > 0)
                        ? (phi_sum / phi_cnt) : amrex::Real(1.0);
                    if (avg_phi <= amrex::Real(0.5)) { return; }  // will become CELL_GAS
                    // avg_phi > 0.5: will become CELL_LIQUID — fall through to refill
                }
            }
            
            // Get bounds for safety
            const auto& farr = f_arrs[nbx];
            const auto lo = amrex::lbound(farr);
            const auto hi = amrex::ubound(farr);
            
            // Step 1 & 2: Sum evs vectors for neighbors that were fluid AND still are fluid
            amrex::Real normal_x = 0.0;
            amrex::Real normal_y = 0.0;
            amrex::Real normal_z = 0.0;
            int num_persistent = 0;
            
            for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                int ni = i + evs[nq][0];
                int nj = j + evs[nq][1];
                int nk = k + evs[nq][2];
                
                // Check bounds
                if (ni < lo.x || ni > hi.x || 
                    nj < lo.y || nj > hi.y || 
                    nk < lo.z || nk > hi.z) continue;
                
                // Check if neighbor was fluid BEFORE and is still fluid NOW
                if (old_fluid_arrs[nbx](ni, nj, nk, 0) == 1 && 
                    curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                    normal_x += evs[nq][0];
                    normal_y += evs[nq][1];
                    normal_z += evs[nq][2];
                    num_persistent++;
                }
            }
            
            // If no persistent neighbors, zero out and exit.
            // This is the validated behaviour from the single-phase case.
            // A copy-from-any-neighbor fallback without a corresponding donor
            // deduction (Step 7 only covers the normal donor path) creates mass
            // every step for blade cells that are simultaneously uncovered.
            if (num_persistent == 0) {
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    f_arrs[nbx](i, j, k, q) = 0.0;
                    g_arrs[nbx](i, j, k, q) = 0.0;
                }
                return;
            }

            // Step 3: Normalize the normal vector
            amrex::Real norm = std::sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
            if (norm == 0.0) {
                // Normal is zero - fall back to first persistent neighbor
                for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                    int ni = i + evs[nq][0];
                    int nj = j + evs[nq][1];
                    int nk = k + evs[nq][2];
                    
                    if (ni < lo.x || ni > hi.x || 
                        nj < lo.y || nj > hi.y || 
                        nk < lo.z || nk > hi.z) continue;
                    
                    if (old_fluid_arrs[nbx](ni, nj, nk, 0) == 1 && 
                        curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                        for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                            f_arrs[nbx](i, j, k, q) = f_arrs[nbx](ni, nj, nk, q);
                            g_arrs[nbx](i, j, k, q) = g_arrs[nbx](ni, nj, nk, q);
                        }
                        return;
                    }
                }
            }
            
            normal_x /= norm;
            normal_y /= norm;
            normal_z /= norm;
            
            // Step 4: Find neighbor with maximum dot product with normal
            amrex::Real max_dot = -1e10;
            int donor_i = -1;
            int donor_j = -1;
            int donor_k = -1;
            
            for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                int ni = i + evs[nq][0];
                int nj = j + evs[nq][1];
                int nk = k + evs[nq][2];
                
                // Check bounds
                if (ni < lo.x || ni > hi.x || 
                    nj < lo.y || nj > hi.y || 
                    nk < lo.z || nk > hi.z) continue;
                
                // Check if neighbor was fluid BEFORE and is still fluid NOW
                if (old_fluid_arrs[nbx](ni, nj, nk, 0) == 1 && 
                    curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                    
                    // Compute dot product
                    amrex::Real dot = evs[nq][0]*normal_x + evs[nq][1]*normal_y + evs[nq][2]*normal_z;
                    
                    if (dot > max_dot) {
                        max_dot = dot;
                        donor_i = ni;
                        donor_j = nj;
                        donor_k = nk;
                    }
                }
            }
            
            // Step 5: Increment donor recipient count
            if (donor_i >= 0) {
                // Atomically increment the count for this donor
                amrex::Gpu::Atomic::Add(&donor_count_arrs[nbx](donor_i, donor_j, donor_k, 0), 1);
            }
        });
    amrex::Gpu::synchronize();
    
    // Synchronize donor counts across ghost cells
    donor_recipient_count.SumBoundary(Geom(lev).periodicity());
    donor_recipient_count.FillBoundary(Geom(lev).periodicity());
    
    // Step 6: Second pass - Transfer ONLY q=0 component from donors to newly-fluid cells
    // Recipients get mass/energy at rest (no momentum/flux), avoiding discontinuities
    // Conservation: donor gives ALL of its q=0 to be shared among N recipients
    {
        auto const& newly_fluid_arrs = newly_fluid_arrs_outer;
        auto const& f_arrs = f_arrs_outer;
        auto const& g_arrs = g_arrs_outer;
        auto const& old_fluid_arrs = old_fluid_arrs_outer;
        auto const& curr_fluid_arrs = curr_fluid_arrs_outer;
        auto const& donor_count_arrs = donor_count_arrs_outer;
        auto const& ct_arrs_ref = ct_arrs_refill;
        const bool   is_free_surface = is_free_surface_refill;
        const stencil::Stencil stencil6;
        const auto& evs = stencil6.evs;
    amrex::ParallelFor(m_f[lev], amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            
            // Only process newly fluid cells
            if (newly_fluid_arrs[nbx](i, j, k, 0) != 1) return;

            // Skip FSLBM-managed cells (same guard as first pass):
            // CELL_GAS, CELL_INTERFACE, CELL_SOLID skipped; CELL_LIQUID handled
            // by body-motion refill as in the validated single-phase case.
            if (is_free_surface) {
                const int ct_val = ct_arrs_ref[nbx](i, j, k, 0);
                if (ct_val == lbm::constants::CELL_GAS      ||
                    ct_val == lbm::constants::CELL_INTERFACE ||
                    ct_val == lbm::constants::CELL_SOLID) { return; }
            }
            
            // Get bounds for safety
            const auto& farr = f_arrs[nbx];
            const auto lo = amrex::lbound(farr);
            const auto hi = amrex::ubound(farr);
            
            // Recompute the normal direction (same as first pass)
            amrex::Real normal_x = 0.0;
            amrex::Real normal_y = 0.0;
            amrex::Real normal_z = 0.0;
            int num_persistent = 0;
            
            for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                int ni = i + evs[nq][0];
                int nj = j + evs[nq][1];
                int nk = k + evs[nq][2];
                
                if (ni < lo.x || ni > hi.x || 
                    nj < lo.y || nj > hi.y || 
                    nk < lo.z || nk > hi.z) continue;
                
                if (old_fluid_arrs[nbx](ni, nj, nk, 0) == 1 && 
                    curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                    normal_x += evs[nq][0];
                    normal_y += evs[nq][1];
                    normal_z += evs[nq][2];
                    num_persistent++;
                }
            }
            
            if (num_persistent == 0) {
                // No persistent neighbors - initialize to zero (same as Step 5)
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    f_arrs[nbx](i, j, k, q) = 0.0;
                    g_arrs[nbx](i, j, k, q) = 0.0;
                }
                return;
            }
            
            // Normalize the normal vector
            amrex::Real norm = std::sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
            if (norm == 0.0) {
                // Fallback: use first persistent neighbor (with /N scaling)
                for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                    int ni = i + evs[nq][0];
                    int nj = j + evs[nq][1];
                    int nk = k + evs[nq][2];
                    
                    if (ni < lo.x || ni > hi.x || 
                        nj < lo.y || nj > hi.y || 
                        nk < lo.z || nk > hi.z) continue;
                    
                    if (old_fluid_arrs[nbx](ni, nj, nk, 0) == 1 && 
                        curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                        
                        int n_recipients = donor_count_arrs[nbx](ni, nj, nk, 0);
                        amrex::Real scale = 1.0 / amrex::Real(n_recipients);
                        f_arrs[nbx](i, j, k, 0) = f_arrs[nbx](ni, nj, nk, 0) * scale;
                        g_arrs[nbx](i, j, k, 0) = g_arrs[nbx](ni, nj, nk, 0) * scale;
                        for (int q = 1; q < constants::N_MICRO_STATES; ++q) {
                            f_arrs[nbx](i, j, k, q) = 0.0;
                            g_arrs[nbx](i, j, k, q) = 0.0;
                        }
                        return;
                    }
                }
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    f_arrs[nbx](i, j, k, q) = 0.0;
                    g_arrs[nbx](i, j, k, q) = 0.0;
                }
                return;
            }
            
            // Find neighbor with maximum dot product with normal
            amrex::Real max_dot = -1e10;
            int donor_i = -1;
            int donor_j = -1;
            int donor_k = -1;
            
            for (int nq = 1; nq < constants::N_MICRO_STATES; ++nq) {
                int ni = i + evs[nq][0];
                int nj = j + evs[nq][1];
                int nk = k + evs[nq][2];
                
                if (ni < lo.x || ni > hi.x || 
                    nj < lo.y || nj > hi.y || 
                    nk < lo.z || nk > hi.z) continue;
                
                if (old_fluid_arrs[nbx](ni, nj, nk, 0) == 1 && 
                    curr_fluid_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                    
                    amrex::Real dot = evs[nq][0]*normal_x + evs[nq][1]*normal_y + evs[nq][2]*normal_z;
                    
                    if (dot > max_dot) {
                        max_dot = dot;
                        donor_i = ni;
                        donor_j = nj;
                        donor_k = nk;
                    }
                }
            }
            
            // Transfer ONLY q=0 component with proper conservation
            if (donor_i >= 0) {
                int n_recipients = donor_count_arrs[nbx](donor_i, donor_j, donor_k, 0);
                // Scale = 1/N to conserve mass (donor gives away everything)
                amrex::Real scale = 1.0 / amrex::Real(n_recipients);
                
                // Recipient gets ONLY q=0 component (mass/energy at rest)
                f_arrs[nbx](i, j, k, 0) = f_arrs[nbx](donor_i, donor_j, donor_k, 0) * scale;
                g_arrs[nbx](i, j, k, 0) = g_arrs[nbx](donor_i, donor_j, donor_k, 0) * scale;
                
                // All other components are zero (no momentum or heat flux)
                for (int q = 1; q < constants::N_MICRO_STATES; ++q) {
                    f_arrs[nbx](i, j, k, q) = 0.0;
                    g_arrs[nbx](i, j, k, q) = 0.0;
                }
            } else {
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    f_arrs[nbx](i, j, k, q) = 0.0;
                    g_arrs[nbx](i, j, k, q) = 0.0;
                }
            }
        });
    amrex::Gpu::synchronize();
    } // end Step 6 block
    
    // Step 7: Third pass - Reduce donor's q=0 component to conserve mass
    // Donors give away ALL of their q=0 component
    {
        auto const& donor_count_arrs = donor_count_arrs_outer;
        auto const& curr_fluid_arrs  = curr_fluid_arrs_outer;
        auto const& f_arrs = f_arrs_outer;
        auto const& g_arrs = g_arrs_outer;
    amrex::ParallelFor(m_f[lev], amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            
            // Check if this cell is a donor (has recipients)
            int n_recipients = donor_count_arrs[nbx](i, j, k, 0);
            if (n_recipients == 0) return;
            
            // Check if this cell is still fluid (donors must be fluid)
            if (curr_fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) != 1) return;
            
            // Donor gives away ALL of q=0
            f_arrs[nbx](i, j, k, 0) = 0.0;
            g_arrs[nbx](i, j, k, 0) = 0.0;
        });
    amrex::Gpu::synchronize();
    } // end Step 7 block

    // Refill for components: DISABLED
    // For passive scalars/tracers, refilling is not physically meaningful.
    // When cells become fluid (e.g., as a moving body rotates away), they should
    // remain at zero or their initial state, not be "refilled" from neighbors.
    // Spilling (when cells become solid) is still performed above to conserve mass.
    //
    // Explicitly zero out component lattices in newly fluid cells
    for (int c = 0; c < m_n_components; ++c) {
        auto const& f_comp_arrs = m_component_lattices[c][lev].arrays();
        auto const& newly_fluid_arrs = newly_fluid.const_arrays();
        
        amrex::ParallelFor(
            m_component_lattices[c][lev], amrex::IntVect(0), constants::N_MICRO_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
                if (newly_fluid_arrs[nbx](i, j, k, 0) == 1) {
                    f_comp_arrs[nbx](i, j, k, q) = 0.0;
                }
            });
        amrex::Gpu::synchronize();
    }

    } // End of refill block

    // Step 8: Reset populations in ALL solid cells (after spill and refill are complete)
    // This ensures no residual populations remain inside the solid body.
    // EXCEPTION: CELL_INTERFACE cells (FSLBM free-surface cells) must retain their
    // f distributions even if the body SDF fraction temporarily classifies them as
    // solid. Those cells are active free-surface cells and zeroing their f creates
    // unphysical all-zero populations that cause negative rho and blow-up.
    {
        auto const& fluid_arrs = m_is_fluid[lev].const_arrays();
        auto const& f_arrs = m_f[lev].arrays();
        auto const& g_arrs = m_g[lev].arrays();
        auto const& md_arrs = m_macrodata[lev].arrays();
        auto const& d_arrs = m_derived[lev].arrays();
        auto const& ct_arrs = m_cell_type[lev].const_arrays();
        
        amrex::ParallelFor(m_f[lev], m_f[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Zero out populations in solid cells — but NOT FSLBM interface cells,
                // which are active free-surface cells with valid f distributions.
                if (fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) == 0
                    && ct_arrs[nbx](i, j, k, 0) != lbm::constants::CELL_INTERFACE) {
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        f_arrs[nbx](i, j, k, q) = 0.0;
                        g_arrs[nbx](i, j, k, q) = 0.0;
                    }
                    
                    // Reset macrodata (check bounds as it might have fewer ghosts)
                    if (md_arrs[nbx].contains(i, j, k)) {
                        for (int n = 0; n < constants::N_MACRO_STATES; ++n) {
                            md_arrs[nbx](i, j, k, n) = 0.0;
                        }
                    }

                    // Reset derived data
                    if (d_arrs[nbx].contains(i, j, k)) {
                        for (int n = 0; n < constants::N_DERIVED; ++n) {
                            d_arrs[nbx](i, j, k, n) = 0.0;
                        }
                    }
                }
            });
    }
    amrex::Gpu::synchronize();

    // Reset components
    for (int c = 0; c < m_n_components; ++c) {
        auto const& fluid_arrs = m_is_fluid[lev].const_arrays();
        auto const& f_comp_arrs = m_component_lattices[c][lev].arrays();

        amrex::ParallelFor(
            m_component_lattices[c][lev],
            m_component_lattices[c][lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Zero out populations in all solid cells
                if (fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) ==
                    0) {
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        f_comp_arrs[nbx](i, j, k, q) = 0.0;
                    }
                }
            });
        amrex::Gpu::synchronize();
    }
    
    // Step 9: Fill boundary cells for updated data
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
    }
    m_g[lev].FillBoundary(Geom(lev).periodicity());
    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());
}

void LBM::reconstruct_body_sdf(const int lev, amrex::Real time)
{
    BL_PROFILE("LBM::reconstruct_body_sdf()");

    if (!m_body_is_moving) return;

    const auto& geom = Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    
    // Update body position and orientation
    amrex::Real dt = time - m_ts_old[lev];
    // If m_ts_old is uninitialized (LOW_NUM), assume dt = 0 (initialization step)
    if (m_ts_old[lev] < -1.0e20) {
        dt = 0.0;
    }

    const amrex::Real vx = m_body_velocity[0];
    const amrex::Real vy = m_body_velocity[1];
    const amrex::Real vz = m_body_velocity[2];
    
    // Update rotation angle (simple Euler integration)
    const amrex::Real omega_mag = std::sqrt(
        m_body_angular_velocity[0] * m_body_angular_velocity[0] +
        m_body_angular_velocity[1] * m_body_angular_velocity[1] +
        m_body_angular_velocity[2] * m_body_angular_velocity[2]);
    
    if (omega_mag > 1e-12) {
        m_body_rotation_angle += omega_mag * dt;
        if (m_print_int > 0 && m_isteps[0] % m_print_int == 0) {
            amrex::Print() << "Updating rotation: dt=" << dt 
                           << " omega=" << omega_mag 
                           << " angle=" << m_body_rotation_angle << std::endl;
        }
    }
    
    // Current body center (for translation)
    const amrex::Real cx = m_body_center[0] + vx * time;
    const amrex::Real cy = m_body_center[1] + vy * time;
    const amrex::Real cz = m_body_center[2] + vz * time;
    
    // Rotation axis (normalized)
    amrex::Real axis_x = 0.0, axis_y = 0.0, axis_z = 1.0;
    if (omega_mag > 1e-12) {
        axis_x = m_body_angular_velocity[0] / omega_mag;
        axis_y = m_body_angular_velocity[1] / omega_mag;
        axis_z = m_body_angular_velocity[2] / omega_mag;
    }
    
    const amrex::Real theta = m_body_rotation_angle;
    const amrex::Real cos_theta = std::cos(theta);
    const amrex::Real sin_theta = std::sin(theta);
    
    // Capture geometry parameters from ParmParse for SDF evaluation
    amrex::Real cyl_radius = 0.1;
    int cyl_direction = 2;
    bool cyl_has_fluid_inside = false;
    
    // Parser support
    amrex::Parser parser;
    amrex::ParserExecutor<3> parser_exe;
    bool use_parser = (m_body_geom_type == "parser");

    // Stationary Parser support
    amrex::Parser parser_stat;
    amrex::ParserExecutor<3> parser_stat_exe;
    bool use_stationary_parser = false;
    {
        amrex::ParmParse pp("eb2");
        std::string stationary_parser_function;
        if (pp.query("stationary_parser_function", stationary_parser_function)) {
            use_stationary_parser = true;
            parser_stat.define(stationary_parser_function);
            parser_stat.registerVariables({"x", "y", "z"});
            parser_stat_exe = parser_stat.compile<3>();
        }
    }

    if (m_isteps[lev] == 0) {
        amrex::Print() << "reconstruct_body_sdf: geom_type='" << m_body_geom_type 
                       << "', use_parser=" << use_parser 
                       << ", use_stationary_parser=" << use_stationary_parser << std::endl;
    }

    if (use_parser) {
        amrex::ParmParse pp("eb2");
        std::string parser_function;
        pp.get("parser_function", parser_function);
        parser.define(parser_function);
        parser.registerVariables({"x", "y", "z"});
        parser_exe = parser.compile<3>();
    } else if (m_body_geom_type == "rotated_cylinder" || m_body_geom_type == "cylinder") {
        amrex::ParmParse pp("eb2");
        pp.query("cylinder_radius", cyl_radius);
        pp.query("cylinder_direction", cyl_direction);
        pp.query("cylinder_has_fluid_inside", cyl_has_fluid_inside);
    }
    
    // Reconstruct fractional field from SDF
    auto const& frac_arrs = m_is_fluid_fraction[lev].arrays();

    // Capture voxel data
    const bool using_voxel_body = m_using_voxel_body;
    const uint16_t* voxel_ptr = using_voxel_body ? m_body_voxel_data.data() : nullptr;
    const amrex::IntVect voxel_dims = m_body_voxel_dims;
    const auto voxel_origin = m_body_voxel_origin;
    const auto voxel_dx = m_body_voxel_dx;
    const auto initial_center = m_body_initial_center;

    // Capture stationary mask
    const bool has_stationary_body = m_has_stationary_body;
    auto const& stat_mask_arrs = m_stationary_mask[lev].const_arrays();
    
    amrex::ParallelFor(
        m_is_fluid_fraction[lev], m_is_fluid_fraction[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // World coordinates of cell center
            const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
            const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
            
            // Translate to body center
            amrex::Real xb = x - cx;
            amrex::Real yb = y - cy;
            amrex::Real zb = z - cz;
            
            // Apply inverse rotation (rotate point in opposite direction)
            // Rodrigues' rotation formula: v_rot = v*cos(θ) + (k×v)*sin(θ) + k(k·v)(1-cos(θ))
            const amrex::Real dot = axis_x * xb + axis_y * yb + axis_z * zb;
            const amrex::Real cross_x = axis_y * zb - axis_z * yb;
            const amrex::Real cross_y = axis_z * xb - axis_x * zb;
            const amrex::Real cross_z = axis_x * yb - axis_y * xb;
            
            // Rotate by -theta (inverse rotation)
            const amrex::Real xr = xb * cos_theta - cross_x * sin_theta + axis_x * dot * (1.0 - cos_theta);
            const amrex::Real yr = yb * cos_theta - cross_y * sin_theta + axis_y * dot * (1.0 - cos_theta);
            const amrex::Real zr = zb * cos_theta - cross_z * sin_theta + axis_z * dot * (1.0 - cos_theta);
            
            // Compute signed distance based on geometry type
            amrex::Real sdf = 0.0;
            
            if (using_voxel_body) {
                 amrex::Real x_init = xr + initial_center[0];
                 amrex::Real y_init = yr + initial_center[1];
                 amrex::Real z_init = zr + initial_center[2];
                 
                 int i_idx = static_cast<int>(std::floor((x_init - voxel_origin[0]) / voxel_dx[0]));
                 int j_idx = static_cast<int>(std::floor((y_init - voxel_origin[1]) / voxel_dx[1]));
                 int k_idx = static_cast<int>(std::floor((z_init - voxel_origin[2]) / voxel_dx[2]));
                 
                 bool is_fluid = true; 
                 
                 if (i_idx >= 0 && i_idx < voxel_dims[0] &&
                     j_idx >= 0 && j_idx < voxel_dims[1] &&
                     k_idx >= 0 && k_idx < voxel_dims[2]) {
                     
                     size_t idx = k_idx * (voxel_dims[0] * voxel_dims[1]) + 
                                  j_idx * voxel_dims[0] + i_idx;
                     
                     // 1=fluid, 0=solid
                     is_fluid = (voxel_ptr[idx] != 0);
                 }
                 
                 // Large value for sharp interface
                 // For solid bodies (default), is_fluid=true means we are outside the body.
                 // The downstream logic expects sdf < 0 for fluid (outside) and sdf > 0 for solid (inside).
                 sdf = is_fluid ? -1.0 : 1.0;

            } else if (use_parser) {
                // Evaluate parser function
                // The parser function defines the shape in the local body frame (centered at 0,0,0).
                // We pass the local coordinates (xr, yr, zr) directly.
                sdf = -parser_exe(xr, yr, zr);
            } else if (cyl_direction == 0) {
                // X-aligned cylinder
                sdf = cyl_radius - std::sqrt(yr * yr + zr * zr);
            } else if (cyl_direction == 1) {
                // Y-aligned cylinder
                sdf = cyl_radius - std::sqrt(xr * xr + zr * zr);
            } else {
                // Z-aligned cylinder (default)
                sdf = cyl_radius - std::sqrt(xr * xr + yr * yr);
            }

            if (use_stationary_parser) {
                // Evaluate stationary parser function in lab frame
                // Note: parser function is expected to be SDF (negative inside, positive outside)
                // But we use -parser() convention here to match the above logic where sdf > 0 is solid.
                // So if user provides standard SDF (neg inside), -SDF is pos inside (solid).
                amrex::Real sdf_stat = -parser_stat_exe(x, y, z);
                sdf = amrex::max(sdf, sdf_stat);
            }

            if (has_stationary_body) {
                // Check mask (STL/CSV)
                // Mask: 1=Fluid, 0=Solid
                // If 0, force sdf to be positive (Solid)
                int is_fluid_stat = stat_mask_arrs[nbx](i, j, k);
                if (is_fluid_stat == 0) {
                    sdf = amrex::max(sdf, 1.0);
                }
            }
            
            // Convert SDF to fractional field
            // Positive SDF = outside (fluid), negative = inside (solid)
            // Use smooth transition with tanh for better numerical behavior
            const amrex::Real interface_width = 1.5 * dx[0];  // ~1.5 cells
            amrex::Real phi;
            
            if (cyl_has_fluid_inside) {
                // Fluid inside, solid outside
                phi = 0.5 * (1.0 + std::tanh(sdf / interface_width));
            } else {
                // Solid inside, fluid outside
                phi = 0.5 * (1.0 + std::tanh(-sdf / interface_width));
            }
            
            // Clamp to [0, 1]
            phi = amrex::max(0.0, amrex::min(1.0, phi));
            
            frac_arrs[nbx](i, j, k, 0) = phi;
        });
    
    amrex::Gpu::synchronize();
    m_is_fluid_fraction[lev].FillBoundary(Geom(lev).periodicity());
}

void LBM::initialize_mask(const int lev)
{
    BL_PROFILE("LBM::initialize_mask()");
    m_mask[lev].setVal(0.0);

    if (lev < finest_level) {
        const amrex::iMultiFab mask = makeFineMask(
            boxArray(lev), DistributionMap(lev), boxArray(lev + 1),
            refRatio(lev));
        amrex::iMultiFab::Copy(
            m_mask[lev], mask, 0, 0, m_mask[lev].nComp(), m_mask[lev].nGrow());
    }
}

void LBM::fill_f_inside_eb(const int lev)
{
    BL_PROFILE("LBM::fill_f_inside_eb()");

    auto const& f_arrs = m_f[lev].arrays();
    auto const& g_arrs = m_g[lev].arrays();

    auto const& is_fluid_arrs = m_is_fluid[lev].arrays();

    amrex::ParallelFor(
        m_f[lev], m_f[lev].nGrowVect(), constants::N_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
            if (is_fluid_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) == 0) {

                f_arrs[nbx](i, j, k, q) = 0.0;
                g_arrs[nbx](i, j, k, q) = 0.0;
            }
        });

    amrex::Gpu::synchronize();
}

// Remake an existing level using provided BoxArray and DistributionMapping
// and fill with existing fine and coarse data.
void LBM::RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("LBM::RemakeLevel()");

    if (Verbose() > 0) {
        amrex::Print() << "Remaking level " << lev << std::endl;
    }

    m_factory[lev] = amrex::makeEBFabFactory(
        Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);
    amrex::MultiFab new_f(
        ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    amrex::Vector<amrex::MultiFab> new_component_lattices(m_n_components);
    for (int i = 0; i < m_n_components; ++i) {
        new_component_lattices[i].define(
            ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
            *(m_factory[lev]));
        new_component_lattices[i].setVal(0.0);
    }
    amrex::MultiFab new_g(
        ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
        *(m_factory[lev]));

    m_fillpatch_op->fillpatch(lev, time, new_f);
    for (int i = 0; i < m_n_components; ++i) {
        m_component_fillpatch_ops[i]->fillpatch(lev, time, new_component_lattices[i]);
    }

    m_fillpatch_g_op->fillpatch(lev, time, new_g);

    std::swap(new_f, m_f[lev]);
    for (int i = 0; i < m_n_components; ++i) {
         std::swap(new_component_lattices[i], m_component_lattices[i][lev]);
    }
    std::swap(new_g, m_g[lev]);

    m_macrodata[lev].define(
        ba, dm, constants::N_MACRO_STATES, m_macrodata_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_is_fluid[lev].define(ba, dm, constants::N_IS_FLUID, m_f[lev].nGrow());
    m_is_fluid_fraction[lev].define(ba, dm, 1, m_is_fluid[lev].nGrow());
    m_eq[lev].define(
        ba, dm, constants::N_MICRO_STATES, m_eq_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_eq_g[lev].define(
        ba, dm, constants::N_MICRO_STATES, m_eq_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_derived[lev].define(
        ba, dm, constants::N_DERIVED, m_derived_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
    m_mask[lev].define(ba, dm, 1, 0);
    m_cell_type[lev].define(ba, dm, 1, m_f_nghost);
    m_cell_type[lev].setVal(constants::CELL_LIQUID);
    m_phi_fslbm[lev].define(ba, dm, 1, m_f_nghost);
    m_phi_fslbm[lev].setVal(amrex::Real(1.0));

    initialize_is_fluid(lev);
    initialize_mask(lev);
    fill_f_inside_eb(lev);
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    for (int i = 0; i < m_n_components; ++i) {
        m_component_lattices[i][lev].FillBoundary(Geom(lev).periodicity());
    }
    m_g[lev].FillBoundary(Geom(lev).periodicity());
    m_macrodata[lev].setVal(0.0);
    m_eq[lev].setVal(0.0);
    m_eq_g[lev].setVal(0.0);
    m_derived[lev].setVal(0.0);

    f_to_macrodata(lev);

    compute_q_corrections(lev);

    macrodata_to_equilibrium(lev);

    compute_derived(lev);

    m_ts_new[lev] = time;
    m_ts_old[lev] = time - constants::SMALL_NUM;
}

// Delete level data
void LBM::ClearLevel(int lev)
{
    BL_PROFILE("LBM::ClearLevel()");
    m_macrodata[lev].clear();
    m_f[lev].clear();
    m_g[lev].clear();
    m_eq[lev].clear();
    m_eq_g[lev].clear();
    m_derived[lev].clear();
    m_is_fluid[lev].clear();
    m_is_fluid_fraction[lev].clear();
    m_plt_mf[lev].clear();
    m_mask[lev].clear();
    m_cell_type[lev].clear();
    m_phi_fslbm[lev].clear();
}

// Set the user defined BC functions
void LBM::set_bcs()
{

    BL_PROFILE("LBM::set_bcs()");
    const bool is_an_energy_lattice(true);
    m_component_fillpatch_ops.resize(m_n_components);

    if (m_velocity_bc_type == "noop") {

        using VelBCOp = bc::BCOpCreator<bc::NoOp>;

        m_fillpatch_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect()), m_f);

        for (int i = 0; i < m_n_components; ++i) {
            m_component_fillpatch_ops[i] = std::make_unique<FillPatchOps<VelBCOp>>(
                geom, refRatio(), m_bcs,
                VelBCOp(m_mesh_speed, m_bc_type, m_component_lattices[i][0].nGrowVect()), m_component_lattices[i]);
        }

        m_fillpatch_g_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(
                m_mesh_speed, m_bc_type, m_g[0].nGrowVect(),
                is_an_energy_lattice),
            m_g);

    } else if (m_velocity_bc_type == "constant") {

        using VelBCOp = bc::BCOpCreator<bc::Constant>;

        m_fillpatch_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect(), "velocity_bc_constant"), m_f);

        for (int i = 0; i < m_n_components; ++i) {
            std::string prefix = "velocity_bc_constant_component_" + std::to_string(i);
            m_component_fillpatch_ops[i] = std::make_unique<FillPatchOps<VelBCOp>>(
                geom, refRatio(), m_bcs,
                VelBCOp(m_mesh_speed, m_bc_type, m_component_lattices[i][0].nGrowVect(), prefix), m_component_lattices[i]);
        }

        m_fillpatch_g_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(
                m_mesh_speed, m_bc_type, m_g[0].nGrowVect(),
                "velocity_bc_constant", is_an_energy_lattice),
            m_g);

    } else if (m_velocity_bc_type == "channel") {

        using VelBCOp = bc::BCOpCreator<bc::Channel>;

        m_fillpatch_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect(), "velocity_bc_channel"), m_f);

        for (int i = 0; i < m_n_components; ++i) {
            std::string prefix = "velocity_bc_channel_component_" + std::to_string(i);
            m_component_fillpatch_ops[i] = std::make_unique<FillPatchOps<VelBCOp>>(
                geom, refRatio(), m_bcs,
                VelBCOp(m_mesh_speed, m_bc_type, m_component_lattices[i][0].nGrowVect(), prefix), m_component_lattices[i]);
        }

        m_fillpatch_g_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(
                m_mesh_speed, m_bc_type, m_g[0].nGrowVect(),
                "velocity_bc_channel", is_an_energy_lattice),
            m_g);

    } else if (m_velocity_bc_type == "parabolic") {

        using VelBCOp = bc::BCOpCreator<bc::Parabolic>;

        m_fillpatch_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect(), "velocity_bc_parabolic"), m_f);

        for (int i = 0; i < m_n_components; ++i) {
            std::string prefix = "velocity_bc_parabolic_component_" + std::to_string(i);
            m_component_fillpatch_ops[i] = std::make_unique<FillPatchOps<VelBCOp>>(
                geom, refRatio(), m_bcs,
                VelBCOp(m_mesh_speed, m_bc_type, m_component_lattices[i][0].nGrowVect(), prefix), m_component_lattices[i]);
        }

        m_fillpatch_g_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(
                m_mesh_speed, m_bc_type, m_g[0].nGrowVect(),
                "velocity_bc_parabolic", is_an_energy_lattice),
            m_g);

    } else {
        amrex::Abort("LBM::set_bcs(): Unknown velocity BC");
    }
}

void LBM::set_ics()
{
    BL_PROFILE("LBM::set_ics()");
    if (m_ic_type == "constant") {
        m_ic_op = std::make_unique<ic::Initializer<ic::Constant>>(
            m_mesh_speed, ic::Constant(ic::Constant()), m_f, m_g);
    } else if (m_ic_type == "taylorgreen") {
        m_ic_op = std::make_unique<ic::Initializer<ic::TaylorGreen>>(
            m_mesh_speed, ic::TaylorGreen(ic::TaylorGreen()), m_f, m_g);
    } else if (m_ic_type == "viscosity_test") {
        m_ic_op = std::make_unique<ic::Initializer<ic::ViscosityTest>>(
            m_mesh_speed, ic::ViscosityTest(ic::ViscosityTest()), m_f, m_g);
    } else if (m_ic_type == "thermaldiffusivity_test") {
        m_ic_op = std::make_unique<ic::Initializer<ic::ThermalDiffusivityTest>>(
            m_mesh_speed,
            ic::ThermalDiffusivityTest(ic::ThermalDiffusivityTest()), m_f, m_g);
    } else if (m_ic_type == "sod") {
        m_ic_op = std::make_unique<ic::Initializer<ic::SodTest>>(
            m_mesh_speed, ic::SodTest(ic::SodTest()), m_f, m_g);
    } else {
        amrex::Abort(
            "LBM::set_ics(): User must specify a valid initial condition");
    }

    m_component_ic_ops.resize(m_n_components);
    for (int i = 0; i < m_n_components; ++i) {
        std::string prefix = "ic_constant_component_" + std::to_string(i);
        if (m_ic_type == "constant") {
            m_component_ic_ops[i] = std::make_unique<ic::Initializer<ic::Constant>>(
                m_mesh_speed, ic::Constant(prefix), m_component_lattices[i]);
        } else {
            // Fallback or error if other IC types are not supported for components yet
            // For now, assume constant IC for components if main IC is constant
             m_component_ic_ops[i] = std::make_unique<ic::Initializer<ic::Constant>>(
                m_mesh_speed, ic::Constant(prefix), m_component_lattices[i]);
        }
    }
}

// Check if a field exists
bool LBM::check_field_existence(const std::string& name)
{
    BL_PROFILE("LBM::check_field_existence()");

    {
        const auto vnames = {
            m_macrodata_varnames, m_microdata_varnames, m_microdata_g_varnames,
            m_deriveddata_varnames, m_idata_varnames, m_fracdata_varnames};

        return std::any_of(vnames.begin(), vnames.end(), [=](const auto& vn) {
            return get_field_component(name, vn) != -1;
        });
    }
}

// Get field component
int LBM::get_field_component(
    const std::string& name, const amrex::Vector<std::string>& varnames)
{
    BL_PROFILE("LBM::get_field_component()");
    const auto itr = std::find(varnames.begin(), varnames.end(), name);
    if (itr != varnames.cend()) {
        return static_cast<int>(std::distance(varnames.begin(), itr));
    }
    return -1;
}

// get a field based on a variable name
std::unique_ptr<amrex::MultiFab>
LBM::get_field(const std::string& name, const int lev, const int ngrow)
{
    BL_PROFILE("LBM::get_field()");

    if (!check_field_existence(name)) {
        amrex::Abort("LBM::get_field(): this field was not found: " + name);
    }

    const int nc = 1;
    std::unique_ptr<amrex::MultiFab> mf = std::make_unique<amrex::MultiFab>(
        boxArray(lev), DistributionMap(lev), nc, ngrow);

    const int srccomp_mad = get_field_component(name, m_macrodata_varnames);
    if (srccomp_mad != -1) {
        amrex::MultiFab::Copy(*mf, m_macrodata[lev], srccomp_mad, 0, nc, ngrow);
    }
    const int srccomp_mid = get_field_component(name, m_microdata_varnames);
    if (srccomp_mid != -1) {
        amrex::MultiFab::Copy(*mf, m_f[lev], srccomp_mid, 0, nc, ngrow);
    }
    const int srccomp_g_mid = get_field_component(name, m_microdata_g_varnames);
    if (srccomp_g_mid != -1) {
        amrex::MultiFab::Copy(*mf, m_g[lev], srccomp_g_mid, 0, nc, ngrow);
    }
    const int srccomp_mdd = get_field_component(name, m_deriveddata_varnames);
    if (srccomp_mdd != -1) {
        amrex::MultiFab::Copy(*mf, m_derived[lev], srccomp_mdd, 0, nc, ngrow);
    }
    const int srccomp_id = get_field_component(name, m_idata_varnames);
    if (srccomp_id != -1) {
        auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
        auto const& mf_arrs = mf->arrays();
        amrex::ParallelFor(
            *mf, mf->nGrowVect(), m_is_fluid[lev].nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                mf_arrs[nbx](i, j, k, n) = is_fluid_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();
    }

    const int srccomp_frac = get_field_component(name, m_fracdata_varnames);
    if (srccomp_frac != -1) {
        auto const& frac_arrs = m_is_fluid_fraction[lev].const_arrays();
        auto const& mf_arrs = mf->arrays();
        amrex::ParallelFor(
            *mf, mf->nGrowVect(), 1,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                mf_arrs[nbx](i, j, k, n) = frac_arrs[nbx](i, j, k, 0);
            });
        amrex::Gpu::synchronize();
    }
    


    amrex::Vector<amrex::BCRec> bcs(nc);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (auto& bc : bcs) {
            bc.setLo(idim, amrex::BCType::foextrap);
            bc.setHi(idim, amrex::BCType::foextrap);
        }
    }
    amrex::FillDomainBoundary(*mf, Geom(lev), bcs);

    return mf;
}


// set covered coarse cells to be the average of overlying fine cells
void LBM::average_down(amrex::IntVect crse_ng)
{
    BL_PROFILE("LBM::average_down()");
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        average_down_to(lev, crse_ng);
    }
}

// more flexible version of AverageDown() that lets you average down across
// multiple levels
void LBM::average_down_to(int crse_lev, amrex::IntVect crse_ng)
{
    BL_PROFILE("LBM::average_down_to()");

    average_down_with_ghosts(
        m_f[crse_lev + 1], m_f[crse_lev], Geom(crse_lev), crse_ng,
        refRatio(crse_lev));

    average_down_with_ghosts(
        m_g[crse_lev + 1], m_g[crse_lev], Geom(crse_lev), crse_ng,
        refRatio(crse_lev));

    amrex::Gpu::synchronize();
}

void LBM::sanity_check_f(const int lev)
{
    BL_PROFILE("LBM::sanity_check_f()");

    for (int q = 0; q < constants::N_MICRO_STATES; q++) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_f[lev].min(q) >= 0.0, "Negative number found in f");
    }
}

// tag cells for refinement
void LBM::ErrorEst(
    int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    BL_PROFILE("LBM::ErrorEst()");

    for (const auto& m_err_tag : m_err_tags) {
        std::unique_ptr<amrex::MultiFab> mf;
        if (!m_err_tag.Field().empty()) {
            mf = get_field(m_err_tag.Field(), lev, m_err_tag.NGrow());
        }
        m_err_tag(
            tags, mf.get(), amrex::TagBox::CLEAR, amrex::TagBox::SET, time, lev,
            Geom(lev));
    }
}

amrex::Vector<std::string> LBM::plot_file_var_names() const
{
    return m_lbm_varnames;
}

std::string LBM::plot_file_name(const int step) const
{
    return amrex::Concatenate(m_plot_file, step, m_file_name_digits);
}

std::string LBM::chk_file_name(const int step) const
{
    return amrex::Concatenate(m_chk_file, step, m_file_name_digits);
}

// put together an array of multifabs for writing
amrex::Vector<const amrex::MultiFab*> LBM::plot_file_mf()
{
    amrex::Vector<const amrex::MultiFab*> r;
    for (int lev = 0; lev <= finest_level; ++lev) {

        m_plt_mf[lev].define(
            boxArray(lev), DistributionMap(lev),
            static_cast<int>(plot_file_var_names().size()), 0);
        int cnt = 0;

        amrex::MultiFab::Copy(
            m_plt_mf[lev], m_macrodata[lev], 0, cnt, m_macrodata[lev].nComp(),
            0);
        cnt += m_macrodata[lev].nComp();

        if (m_save_streaming) {
            amrex::MultiFab::Copy(
                m_plt_mf[lev], m_f[lev], 0, cnt, m_f[lev].nComp(), 0);
            cnt += m_f[lev].nComp();
        }

        if (m_save_streaming) {
            amrex::MultiFab::Copy(
                m_plt_mf[lev], m_g[lev], 0, cnt, m_g[lev].nComp(), 0);
            cnt += m_g[lev].nComp();
        }

        if (m_save_derived) {
            amrex::MultiFab::Copy(
                m_plt_mf[lev], m_derived[lev], 0, cnt, m_derived[lev].nComp(),
                0);
            cnt += m_derived[lev].nComp();
        }
        auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
        auto const& plt_mf_arrs = m_plt_mf[lev].arrays();
        amrex::ParallelFor(
            m_plt_mf[lev], m_plt_mf[lev].nGrowVect(), m_is_fluid[lev].nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                plt_mf_arrs[nbx](i, j, k, n + cnt) =
                    is_fluid_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();
        cnt += m_is_fluid[lev].nComp();
        // copy fractional field (1 component)
        auto const& frac_arrs = m_is_fluid_fraction[lev].const_arrays();
        amrex::ParallelFor(
            m_plt_mf[lev], m_plt_mf[lev].nGrowVect(), 1,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                plt_mf_arrs[nbx](i, j, k, n + cnt) = frac_arrs[nbx](i, j, k, 0);
            });
        amrex::Gpu::synchronize();
        cnt += 1;

        auto const& md_arrs = m_macrodata[lev].const_arrays();
        for (int c = 0; c < m_n_components; ++c) {
            auto const& f_comp_arrs = m_component_lattices[c][lev].const_arrays();
            amrex::ParallelFor(
                m_plt_mf[lev], m_plt_mf[lev].nGrowVect(),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                    // Only report Y_k in fluid cells; gas and solid cells
                    // may hold stale populations that never contribute to
                    // transport and should appear as zero in the plotfile.
                    if (is_fluid_arrs[nbx](i, j, k,
                            lbm::constants::IS_FLUID_IDX) != 1) {
                        plt_mf_arrs[nbx](i, j, k, cnt) = 0.0;
                        return;
                    }
                    amrex::Real rho_comp = 0.0;
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        rho_comp += f_comp_arrs[nbx](i, j, k, q);
                    }
                    amrex::Real rho_total =
                        md_arrs[nbx](i, j, k, constants::RHO_IDX);
                    amrex::Real Y_k =
                        (rho_total > 0.0) ? (rho_comp / rho_total) : 0.0;
                    plt_mf_arrs[nbx](i, j, k, cnt) = Y_k;
                });
            amrex::Gpu::synchronize();
            cnt += 1;
        }

        r.push_back(&m_plt_mf[lev]);
    }
    return r;
}

void LBM::write_plot_file()
{
    BL_PROFILE("LBM::write_plot_file()");
    const std::string& plotfilename = plot_file_name(m_isteps[0]);
    const auto& mf = plot_file_mf();
    const auto& varnames = plot_file_var_names();

    amrex::Print() << "Writing plot file " << plotfilename << " at time "
                   << m_ts_new[0] << std::endl;

    amrex::WriteMultiLevelPlotfile(
        plotfilename, finest_level + 1, mf, varnames, Geom(), m_ts_new[0],
        m_isteps, refRatio());

    // Write Lagrangian bubble particles into the same plotfile directory so
    // ParaView's "AMReX/BoxLib Particles Reader" can load them alongside the
    // mesh fields.  The subdirectory will be plt00000/Bubbles/.
    if (m_enable_bubbles) {
        // Before writing, convert rdata to more meaningful output units:
        //   diameter : SI [m]  → LB cells  (scale factor 1.0 in ParaView)
        //   n_o2     : mol     → C_g = n_O2/V_b [mol/m³]  (starts ~44.6, drops to 0)
        const amrex::Real inv_dx  = 1.0 / m_bubble_params.dx_phys;
        const amrex::Real pi_over_6 = amrex::Math::pi<amrex::Real>() / 6.0;
        auto& container = m_bubbles.container();
        for (int lev = 0; lev <= container.finestLevel(); ++lev) {
            for (auto& kv : container.GetParticles(lev)) {
                for (auto& p : kv.second.GetArrayOfStructs()()) {
                    if (!p.id().is_valid()) continue;
                    const amrex::Real d  = p.rdata(lbm::BubbleIdx::DIAMETER);  // SI [m]
                    const amrex::Real Vb = pi_over_6 * d * d * d;              // m³
                    // n_o2 → C_g [mol/m³]
                    p.rdata(lbm::BubbleIdx::N_O2) = (Vb > 0.0)
                        ? p.rdata(lbm::BubbleIdx::N_O2) / Vb : 0.0;
                    // diameter → LB cells
                    p.rdata(lbm::BubbleIdx::DIAMETER) *= inv_dx;
                }
            }
        }
        container.WritePlotFile(
            plotfilename, "Bubbles",
            {"vx", "vy", "vz", "diameter", "C_g_mol_m3", "ax", "ay", "az"});
        // Write the simulation time into plt*/Bubbles/time so that the ParaView
        // AMReX Grid Reader reports the same time for both the fluid fields and
        // the bubble particles.  The AMReX particle sub-header (plt*/Bubbles/Header)
        // contains no time stamp; without this sidecar file ParaView infers the
        // particle time from the directory name (step number as an integer) while
        // the fluid reader reads the stored float time from plt*/Header — producing
        // two incoherent time axes and doubling the apparent step count.
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream tfile(plotfilename + "/Bubbles/time");
            tfile << std::setprecision(17) << m_ts_new[0] << '\n';
        }
        // Restore original rdata (diameter → SI, C_g → n_o2)
        const amrex::Real dx = m_bubble_params.dx_phys;
        for (int lev = 0; lev <= container.finestLevel(); ++lev) {
            for (auto& kv : container.GetParticles(lev)) {
                for (auto& p : kv.second.GetArrayOfStructs()()) {
                    if (!p.id().is_valid()) continue;
                    // diameter: LB cells → SI [m]
                    p.rdata(lbm::BubbleIdx::DIAMETER) *= dx;
                    const amrex::Real d  = p.rdata(lbm::BubbleIdx::DIAMETER);
                    const amrex::Real Vb = pi_over_6 * d * d * d;
                    // C_g → n_o2 [mol]
                    p.rdata(lbm::BubbleIdx::N_O2) *= Vb;
                }
            }
        }
    }
}

void LBM::write_checkpoint_file() const
{
    BL_PROFILE("LBM::write_checkpoint_file()");
    const auto& varnames = m_microdata_varnames;
    const auto& varnames_g = m_microdata_g_varnames;
    const auto& varnames_frac = m_fracdata_varnames;

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g.,
    // finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data
    // at each level of refinement

    const std::string& checkpointname = chk_file_name(m_isteps[0]);

    amrex::Print() << "Writing checkpoint file " << checkpointname
                   << " at time " << m_ts_new[0] << std::endl;

    const int nlevels = finest_level + 1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then
    // build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
    if (amrex::ParallelDescriptor::IOProcessor()) {

        const std::string header_file_name(checkpointname + "/Header");
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream header_file;
        header_file.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        header_file.open(
            header_file_name.c_str(),
            std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if (!header_file.good()) {
            amrex::FileOpenFailed(header_file_name);
        }

        header_file.precision(17);

        // write out title line
        header_file << "Checkpoint file for LBM\n";

        // write out finest_level
        header_file << finest_level << "\n";

        // write out array of istep
        for (int m_istep : m_isteps) {
            header_file << m_istep << " ";
        }
        header_file << "\n";

        // write out array of dt
        for (double m_dt : m_dts) {
            header_file << m_dt << " ";
        }
        header_file << "\n";

        // write out array of t_new
        for (double i : m_ts_new) {
            header_file << i << " ";
        }
        header_file << "\n";

        // write the BoxArray at each level
        for (int lev = 0; lev <= finest_level; ++lev) {
            boxArray(lev).writeOn(header_file);
            header_file << '\n';
        }
    }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Write(
            m_f[lev], amrex::MultiFabFileFullPrefix(
                          lev, checkpointname, "Level_", varnames[0]));
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Write(
            m_g[lev], amrex::MultiFabFileFullPrefix(
                          lev, checkpointname, "Level_", varnames_g[0]));
    }

    for (int c = 0; c < m_n_components; ++c) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::string mf_name = "f_comp_" + std::to_string(c);
            amrex::VisMF::Write(
                m_component_lattices[c][lev],
                amrex::MultiFabFileFullPrefix(
                    lev, checkpointname, "Level_", mf_name));
        }
    }

    // write fractional is_fluid field
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Write(
            m_is_fluid_fraction[lev], amrex::MultiFabFileFullPrefix(
                                          lev, checkpointname, "Level_",
                                          varnames_frac[0]));
    }
}

void LBM::read_checkpoint_file()
{
    BL_PROFILE("LBM::read_checkpoint_file()");
    const auto& varnames = m_microdata_varnames;
    const auto& varnames_g = m_microdata_g_varnames;
    const auto& varnames_frac = m_fracdata_varnames;

    amrex::Print() << "Restarting from checkpoint file " << m_restart_chkfile
                   << std::endl;

    // Header
    const std::string file(m_restart_chkfile + "/Header");

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

    amrex::Vector<char> file_char_ptr;
    amrex::ParallelDescriptor::ReadAndBcastFile(file, file_char_ptr);
    std::string file_char_ptr_string(file_char_ptr.dataPtr());
    std::istringstream is(file_char_ptr_string, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    goto_next_line(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_isteps[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_dts[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_ts_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        amrex::BoxArray ba;
        ba.readFrom(is);
        goto_next_line(is);

        // create a distribution mapping
        amrex::DistributionMapping dm{ba, amrex::ParallelDescriptor::NProcs()};

        // set BoxArray grids and DistributionMapping dmap in
        // AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFabs
        const int ncomp = static_cast<int>(varnames.size());
        AMREX_ASSERT(ncomp == constants::N_MICRO_STATES);
        m_factory[lev] = amrex::makeEBFabFactory(
            Geom(lev), ba, dm, {5, 5, 5}, amrex::EBSupport::basic);
        m_f[lev].define(
            ba, dm, ncomp, m_f_nghost, amrex::MFInfo(), *(m_factory[lev]));
        m_g[lev].define(
            ba, dm, ncomp, m_f_nghost, amrex::MFInfo(), *(m_factory[lev]));
        for (int i = 0; i < m_n_components; ++i) {
            m_component_lattices[i][lev].define(
                ba, dm, ncomp, m_f_nghost, amrex::MFInfo(), *(m_factory[lev]));
        }
        m_macrodata[lev].define(
            ba, dm, constants::N_MACRO_STATES, m_macrodata_nghost,
            amrex::MFInfo(), *(m_factory[lev]));
        m_is_fluid[lev].define(ba, dm, constants::N_IS_FLUID, m_f[lev].nGrow());
        m_eq[lev].define(
            ba, dm, constants::N_MICRO_STATES, m_eq_nghost, amrex::MFInfo(),
            *(m_factory[lev]));
        m_eq_g[lev].define(
            ba, dm, constants::N_MICRO_STATES, m_eq_nghost, amrex::MFInfo(),
            *(m_factory[lev]));
        m_derived[lev].define(
            ba, dm, constants::N_DERIVED, m_derived_nghost, amrex::MFInfo(),
            *(m_factory[lev]));
        m_mask[lev].define(ba, dm, 1, 0);
        // define fractional field storage
        m_is_fluid_fraction[lev].define(ba, dm, 1, m_is_fluid[lev].nGrow());
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
            m_f[lev], amrex::MultiFabFileFullPrefix(
                          lev, m_restart_chkfile, "Level_", varnames[0]));
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
            m_g[lev], amrex::MultiFabFileFullPrefix(
                          lev, m_restart_chkfile, "Level_", varnames_g[0]));
    }

    for (int c = 0; c < m_n_components; ++c) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::string mf_name = "f_comp_" + std::to_string(c);
            amrex::VisMF::Read(
                m_component_lattices[c][lev],
                amrex::MultiFabFileFullPrefix(
                    lev, m_restart_chkfile, "Level_", mf_name));
        }
    }

    // read fractional is_fluid field
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
            m_is_fluid_fraction[lev], amrex::MultiFabFileFullPrefix(
                                          lev, m_restart_chkfile, "Level_",
                                          varnames_frac[0]));
    }

    // Populate the other data
    for (int lev = 0; lev <= finest_level; ++lev) {
        initialize_is_fluid(lev);
        initialize_mask(lev);
        fill_f_inside_eb(lev);
        m_f[lev].FillBoundary(Geom(lev).periodicity());
        m_g[lev].FillBoundary(Geom(lev).periodicity());
        m_macrodata[lev].setVal(0.0);
        m_eq[lev].setVal(0.0);
        m_eq_g[lev].setVal(0.0);
        m_derived[lev].setVal(0.0);

        f_to_macrodata(lev);

        compute_q_corrections(lev);

        macrodata_to_equilibrium(lev);

        compute_derived(lev);
    }
}

// utility to skip to next line in Header
void LBM::goto_next_line(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}

void LBM::open_forces_file(const bool initialize)
{
    BL_PROFILE("LBM::open_forces_file()");
    if (m_compute_forces) {
        if ((file_exists(m_forces_file)) && (!initialize)) {
            m_forces_stream.open(m_forces_file, std::ios::app);
        } else {
            m_forces_stream.open(m_forces_file, std::ios::out);
            m_forces_stream << std::setw(constants::DATWIDTH)
                            << "          time";
            AMREX_D_DECL(
                m_forces_stream << std::setw(constants::DATWIDTH)
                                << "          fx",
                m_forces_stream << std::setw(constants::DATWIDTH)
                                << "          fy",
                m_forces_stream << std::setw(constants::DATWIDTH)
                                << "          fz");
            m_forces_stream << std::endl;
        }
    }
}

void LBM::close_forces_file()
{
    BL_PROFILE("LBM::close_forces_file()");
    if (m_forces_stream) {
        m_forces_stream.close();
    }
}

void LBM::output_forces_file(const amrex::Vector<amrex::Real>& forces)
{
    BL_PROFILE("LBM::output_forces_file()");
    if (m_compute_forces) {
        m_forces_stream << std::setw(constants::DATWIDTH)
                        << std::setprecision(constants::DATPRECISION)
                        << m_ts_new[0];
        for (const auto& val : forces) {
            m_forces_stream << std::setw(constants::DATWIDTH)
                            << std::setprecision(constants::DATPRECISION)
                            << val;
        }
        m_forces_stream << std::endl;
    }
}

// ---------------------------------------------------------------
// Two-step catalytic reaction: S + C <--(k_f/k_r)--> I --k_p--> P + C
//
// Operator-split source terms applied after each BGK collision step.
// The change in local density for each species over one time step is:
//
//   R_fwd = k_f * rho_S * rho_C    (bimolecular)
//   R_rev = k_r * rho_I             (unimolecular)
//   R_prd = k_p * rho_I             (unimolecular)
//
//   Delta_rho_S = R_rev - R_fwd
//   Delta_rho_C = (R_rev + R_prd) - R_fwd
//   Delta_rho_I = R_fwd - (R_rev + R_prd)
//   Delta_rho_P = R_prd
//
// Each increment is spread uniformly over all N_MICRO_STATES populations
// (isotropic source => no spurious momentum injection).
// ---------------------------------------------------------------
void LBM::apply_reaction_source_terms(const int lev)
{
    BL_PROFILE("LBM::apply_reaction_source_terms()");

    // Need exactly four components: S(0), C(1), I(2), P(3)
    AMREX_ASSERT(m_n_components >= 4);

    const amrex::Real k_f  = m_rxn_k_forward;
    const amrex::Real k_r  = m_rxn_k_reverse;
    const amrex::Real k_p  = m_rxn_k_product;

    // Thermodynamic constants needed for the equilibrium distribution
    const amrex::Real specific_gas_constant = m_R_u / m_m_bar;
    const amrex::Real l_mesh_speed          = m_mesh_speed;

    // Stencil: lattice velocities and weights
    const stencil::Stencil stencil;
    const auto& evs    = stencil.evs;
    const auto& weight = stencil.weights;

    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& md_arrs       = m_macrodata[lev].const_arrays();

    // Non-const arrays for the four reactive components
    auto const& fS_arrs = m_component_lattices[0][lev].arrays();
    auto const& fC_arrs = m_component_lattices[1][lev].arrays();
    auto const& fI_arrs = m_component_lattices[2][lev].arrays();
    auto const& fP_arrs = m_component_lattices[3][lev].arrays();

    amrex::ParallelFor(
        m_component_lattices[0][lev],
        amrex::IntVect(0),   // no ghost cells — source terms only on interior
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept
        {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) != 1) {
                return;
            }

            // Compute local species densities by summing over populations
            amrex::Real rho_S = 0.0, rho_C = 0.0, rho_I = 0.0, rho_P = 0.0;
            for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                rho_S += fS_arrs[nbx](iv, q);
                rho_C += fC_arrs[nbx](iv, q);
                rho_I += fI_arrs[nbx](iv, q);
                rho_P += fP_arrs[nbx](iv, q);
            }

            // Guard against numerical noise producing negative densities
            rho_S = amrex::max(rho_S, 0.0);
            rho_C = amrex::max(rho_C, 0.0);
            rho_I = amrex::max(rho_I, 0.0);

            // Reaction rates
            const amrex::Real R_fwd = k_f * rho_S * rho_C;
            const amrex::Real R_rev = k_r * rho_I;
            const amrex::Real R_prd = k_p * rho_I;

            // Density increments
            const amrex::Real d_S = R_rev - R_fwd;
            const amrex::Real d_C = (R_rev + R_prd) - R_fwd;
            const amrex::Real d_I = R_fwd - (R_rev + R_prd);
            const amrex::Real d_P = R_prd;

            // Local fluid velocity and temperature from macrodata
            const amrex::RealVect vel = {AMREX_D_DECL(
                md_arrs[nbx](iv, constants::VELX_IDX),
                md_arrs[nbx](iv, constants::VELY_IDX),
                md_arrs[nbx](iv, constants::VELZ_IDX))};
            const amrex::Real temperature = md_arrs[nbx](iv, constants::TEMPERATURE_IDX);

            // Equilibrium stress components (purely kinetic — no viscous correction)
            const amrex::Real Rg_T  = specific_gas_constant * temperature;
            const amrex::Real pxx_eq = vel[0] * vel[0] + Rg_T;
            const amrex::Real pyy_eq = vel[1] * vel[1] + Rg_T;
            const amrex::Real pzz_eq = AMREX_D_PICK(0.0, 0.0, vel[2] * vel[2] + Rg_T);

            // Distribute density increments using the local equilibrium shape.
            // set_extended_equilibrium_value is linear in rho, so the unit
            // equilibrium f_eq(rho=1, vel, T) gives the correct per-q weight.
            // Adding d_X * f_eq_unit to each population preserves the net
            // momentum added by the reaction (zero here since d_X is scalar).
            for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                const amrex::Real wt = weight[q];
                const auto& ev       = evs[q];
                const amrex::Real f_eq_unit = set_extended_equilibrium_value(
                    1.0, vel, pxx_eq, pyy_eq, pzz_eq, l_mesh_speed, wt, ev);
                fS_arrs[nbx](iv, q) += d_S * f_eq_unit;
                fC_arrs[nbx](iv, q) += d_C * f_eq_unit;
                fI_arrs[nbx](iv, q) += d_I * f_eq_unit;
                fP_arrs[nbx](iv, q) += d_P * f_eq_unit;
            }
        });
    amrex::Gpu::synchronize();
}

// ---------------------------------------------------------------
// Timed catalyst injection
// On the first call where m_isteps[0] >= m_cat_inject_step, fill all
// fluid cells whose physical centres fall inside the injection box with
// a uniform equilibrium-like population for component 1 (catalyst C).
// The flag m_cat_inject_done prevents any repeated application.
// ---------------------------------------------------------------
void LBM::apply_timed_catalyst_injection(const int lev)
{
    BL_PROFILE("LBM::apply_timed_catalyst_injection()");

    if (m_cat_inject_done || m_cat_inject_step < 0) { return; }
    if (m_isteps[0] < m_cat_inject_step)             { return; }
    if (m_n_components < 2)                           { return; }

    m_cat_inject_done = true;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "\n[Reaction] Injecting catalyst (component 1) at step "
                       << m_isteps[0] << "  (LB time = " << m_ts_new[lev] << ")\n";
    }

    const amrex::Real rho_inject = m_cat_inject_density;
    const amrex::Real pop_val    = rho_inject / constants::N_MICRO_STATES;

    const amrex::RealVect box_lo = m_cat_inject_box_lo;
    const amrex::RealVect box_hi = m_cat_inject_box_hi;

    const auto prob_lo = Geom(lev).ProbLoArray();
    const auto dx      = Geom(lev).CellSizeArray();

    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& fC_arrs       = m_component_lattices[1][lev].arrays();

    amrex::ParallelFor(
        m_component_lattices[1][lev],
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept
        {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) != 1) {
                return;
            }
            // Physical cell centre
            const amrex::Real cx = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real cy = prob_lo[1] + (j + 0.5) * dx[1];
#if AMREX_SPACEDIM == 3
            const amrex::Real cz = prob_lo[2] + (k + 0.5) * dx[2];
#else
            const amrex::Real cz = 0.0;
            (void)(cz); // suppress unused-variable warning in 2D builds
#endif

            const bool inside =
                (cx >= box_lo[0] && cx <= box_hi[0]) &&
                (cy >= box_lo[1] && cy <= box_hi[1])
#if AMREX_SPACEDIM == 3
                && (cz >= box_lo[2] && cz <= box_hi[2])
#endif
                ;

            if (inside) {
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    fC_arrs[nbx](iv, q) = pop_val;
                }
            }
        });
    amrex::Gpu::synchronize();

    // Ensure ghost cells are consistent after the injection
    m_component_lattices[1][lev].FillBoundary(Geom(lev).periodicity());
}

// ---------------------------------------------------------------
// Species statistics file (species_stats.csv)
// Columns: step, LBtime, mean_rho_S, mean_rho_C, mean_rho_I, mean_rho_P,
//          cat_cycle_sum (= rho_C + rho_I, should be const after injection),
//          sub_cycle_sum (= rho_S + rho_I + rho_P, should be const)
// ---------------------------------------------------------------
void LBM::open_species_stats_file(const bool initialize)
{
    BL_PROFILE("LBM::open_species_stats_file()");
    if (!(m_enable_reactions && m_n_components >= 4)) { return; }

    const bool file_already_exists = file_exists(m_species_stats_file);
    if (file_already_exists && !initialize) {
        m_species_stats_stream.open(m_species_stats_file, std::ios::app);
    } else {
        m_species_stats_stream.open(m_species_stats_file, std::ios::out);
        // Header
        m_species_stats_stream
            << std::setw(12) << "step"
            << std::setw(constants::DATWIDTH) << "LBtime"
            << std::setw(constants::DATWIDTH) << "mean_rho_S"
            << std::setw(constants::DATWIDTH) << "mean_rho_C"
            << std::setw(constants::DATWIDTH) << "mean_rho_I"
            << std::setw(constants::DATWIDTH) << "mean_rho_P"
            << std::setw(constants::DATWIDTH) << "cat_cycle_sum"
            << std::setw(constants::DATWIDTH) << "sub_cycle_sum"
            << "\n";
    }
}

void LBM::close_species_stats_file()
{
    BL_PROFILE("LBM::close_species_stats_file()");
    if (m_species_stats_stream.is_open()) {
        m_species_stats_stream.close();
    }
}

void LBM::write_species_stats()
{
    BL_PROFILE("LBM::write_species_stats()");
    if (!m_species_stats_stream.is_open()) { return; }

    // Volume-average each species density over all fluid cells on level 0.
    // We reduce four sums (one per component) plus a fluid-cell counter.
    amrex::Real sum_S = 0.0, sum_C = 0.0, sum_I = 0.0, sum_P = 0.0;
    amrex::Real n_fluid = 0.0;

    auto const& is_fluid_arrs = m_is_fluid[0].const_arrays();

    // Helper lambda: reduce total density of one component lattice.
    auto reduce_component = [&](int c) -> amrex::Real {
        amrex::Real total = 0.0;
        auto const& fc_arrs = m_component_lattices[c][0].const_arrays();
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(
            m_component_lattices[c][0], amrex::IntVect(0),
            reduce_data,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) -> ReduceTuple
            {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) != 1) {
                    return {0.0};
                }
                amrex::Real rho = 0.0;
                const auto fc = fc_arrs[nbx];
                for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                    rho += fc(iv, q);
                }
                return {rho};
            });
        ReduceTuple hv = reduce_data.value(reduce_op);
        total = amrex::get<0>(hv);
        amrex::ParallelDescriptor::ReduceRealSum(total);
        return total;
    };

    // Count fluid cells (only needs to be done once; reuse for all species)
    {
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(
            m_component_lattices[0][0], amrex::IntVect(0), reduce_data,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) -> ReduceTuple
            {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                return {amrex::Real(is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) == 1 ? 1.0 : 0.0)};
            });
        ReduceTuple hv = reduce_data.value(reduce_op);
        n_fluid = amrex::get<0>(hv);
        amrex::ParallelDescriptor::ReduceRealSum(n_fluid);
    }

    if (n_fluid <= 0.5) { return; }

    sum_S = reduce_component(0);
    sum_C = reduce_component(1);
    sum_I = reduce_component(2);
    sum_P = reduce_component(3);

    const amrex::Real inv_n = 1.0 / n_fluid;
    const amrex::Real mean_S = sum_S * inv_n;
    const amrex::Real mean_C = sum_C * inv_n;
    const amrex::Real mean_I = sum_I * inv_n;
    const amrex::Real mean_P = sum_P * inv_n;

    // Conservation diagnostics
    // cat_cycle_sum  = mean_C + mean_I  (constant after catalyst injection)
    // sub_cycle_sum  = mean_S + mean_I + mean_P  (constant everywhere)
    const amrex::Real cat_cycle = mean_C + mean_I;
    const amrex::Real sub_cycle = mean_S + mean_I + mean_P;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        m_species_stats_stream
            << std::setw(12) << m_isteps[0]
            << std::setw(constants::DATWIDTH)
            << std::setprecision(constants::DATPRECISION) << m_ts_new[0]
            << std::setw(constants::DATWIDTH)
            << std::setprecision(constants::DATPRECISION) << mean_S
            << std::setw(constants::DATWIDTH)
            << std::setprecision(constants::DATPRECISION) << mean_C
            << std::setw(constants::DATWIDTH)
            << std::setprecision(constants::DATPRECISION) << mean_I
            << std::setw(constants::DATWIDTH)
            << std::setprecision(constants::DATPRECISION) << mean_P
            << std::setw(constants::DATWIDTH)
            << std::setprecision(constants::DATPRECISION) << cat_cycle
            << std::setw(constants::DATWIDTH)
            << std::setprecision(constants::DATPRECISION) << sub_cycle
            << "\n";
        m_species_stats_stream.flush();
    }
}

#if 0  // FSLBM (Körner 2005) replaces advance_phi; kept for reference only.
// ============================================================================
// advance_phi — one step of the Chiu & Lin (2011) conservative phase-field
// equation for the liquid-gas free surface.
//
// Governing equation (Chiu & Lin Eq. 18):
//   dPhi/dt + div(u * Phi) = gamma * div[ grad(Phi) - Phi*(1-Phi)/eps * nhat ]
// where
//   nhat  = grad(Phi) / |grad(Phi)|   (unit normal, regularised)
//   gamma = gamma_coeff * |u|          (interface compression coefficient)
//   eps   = 1.5 * dx                   (interface half-width, set at init)
//
// Algorithm (per step):
//   1. Compute cell-centred compression vector C = Phi*(1-Phi)/eps * nhat
//      using central diffs from 1-ghost-cell halo of m_is_fluid_fraction.
//   2. FillBoundary(C) so its divergence at interior cells has valid neighbours.
//   3. Forward-Euler update:
//        Phi_new = Phi + dt * [ -upwind_div(u,Phi)
//                               + gamma * (Laplacian(Phi) - div(C)) ]
//   4. Mass redistribution (Sec. 2.5): clip to [0,1], compute
//      G = M0 - M, add G/N_interface to transition cells (0.001 < Phi < 0.999).
//   5. Copy Phi_new -> m_is_fluid_fraction; call refill_and_spill to update
//      the LBM fluid/gas domain (identical path to moving solid body).
// ============================================================================
void LBM::advance_phi(const int lev)
{
    BL_PROFILE("LBM::advance_phi()");

    // Zero-gradient (Neumann) fill for domain boundary ghost cells.
    // FillBoundary(periodicity()) only fills periodic faces; on a
    // non-periodic domain all boundary ghost cells remain at their
    // allocation value (0), causing spurious gradients in the stencil.
    // This lambda copies the nearest interior cell into each ghost cell
    // that lies outside the domain — equivalent to a 90° contact-angle BC.
    auto fill_neumann_bc = [](amrex::MultiFab& mf, const amrex::Geometry& geom) {
        const amrex::Box& domain = geom.Domain();
        const int ncomp = mf.nComp();
        auto const& arrs = mf.arrays();
        amrex::ParallelFor(mf, mf.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const int ii = amrex::max(domain.smallEnd(0), amrex::min(domain.bigEnd(0), i));
                const int jj = amrex::max(domain.smallEnd(1), amrex::min(domain.bigEnd(1), j));
                const int kk = amrex::max(domain.smallEnd(2), amrex::min(domain.bigEnd(2), k));
                if (ii != i || jj != j || kk != k) {
                    for (int c = 0; c < ncomp; ++c) {
                        arrs[nbx](i, j, k, c) = arrs[nbx](ii, jj, kk, c);
                    }
                }
            });
        amrex::Gpu::synchronize();
    };

    const auto& geom   = Geom(lev);
    const auto  dx_arr = geom.CellSizeArray();
    const amrex::Real inv2dx  = 0.5 / dx_arr[0];
    const amrex::Real inv2dy  = 0.5 / dx_arr[1];
    const amrex::Real inv2dz  = 0.5 / dx_arr[2];
    const amrex::Real inv_dx2 = 1.0 / (dx_arr[0] * dx_arr[0]);
    const amrex::Real inv_dy2 = 1.0 / (dx_arr[1] * dx_arr[1]);
    const amrex::Real inv_dz2 = 1.0 / (dx_arr[2] * dx_arr[2]);
    const amrex::Real dx0     = dx_arr[0];
    const amrex::Real dy0     = dx_arr[1];
    const amrex::Real dz0     = dx_arr[2];

    const amrex::Real eps_phi    = m_free_surface_eps;
    const amrex::Real gamma_coef = m_phi_gamma_coeff;
    const amrex::Real dt         = m_dts[lev];
    const amrex::Real nhat_reg   = 1.0e-8;   // |grad_phi| regulariser

    // ------------------------------------------------------------------
    // Ghost-fill Phi so the stencil has valid halo values.
    // ------------------------------------------------------------------
    m_is_fluid_fraction[lev].FillBoundary(geom.periodicity());
    fill_neumann_bc(m_is_fluid_fraction[lev], geom);

    auto const& phi_arrs = m_is_fluid_fraction[lev].const_arrays();
    auto const& md_arrs  = m_macrodata[lev].const_arrays();

    // ------------------------------------------------------------------
    // Step 1: compression vector  C = Phi*(1-Phi)/eps * grad(Phi)/|grad(Phi)|
    // Computed at every cell (interior + 1-ghost layer).
    // ------------------------------------------------------------------
    amrex::MultiFab phi_comp(
        m_is_fluid_fraction[lev].boxArray(),
        m_is_fluid_fraction[lev].DistributionMap(),
        AMREX_SPACEDIM, 1);
    phi_comp.setVal(0.0);

    {
        auto const& c_arrs = phi_comp.arrays();
        auto const& stat_c = m_stationary_mask[lev].const_arrays();
        amrex::ParallelFor(
            phi_comp, amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Solid cells: compression vector = 0 (no interface here)
                if (stat_c[nbx](i, j, k) == 0) {
                    c_arrs[nbx](i, j, k, 0) = 0.0;
                    c_arrs[nbx](i, j, k, 1) = 0.0;
                    c_arrs[nbx](i, j, k, 2) = 0.0;
                    return;
                }
                const amrex::Real phi = phi_arrs[nbx](i, j, k, 0);

                // Neumann BC at solid faces: substitute current cell's phi for
                // any solid neighbor to eliminate spurious contact-line gradients.
                auto phi_nbr = [&](int ii, int jj, int kk) -> amrex::Real {
                    return (stat_c[nbx](ii, jj, kk) == 0) ? phi
                           : phi_arrs[nbx](ii, jj, kk, 0);
                };
                const amrex::Real gpx = (phi_nbr(i+1,j,k) - phi_nbr(i-1,j,k)) * inv2dx;
                const amrex::Real gpy = (phi_nbr(i,j+1,k) - phi_nbr(i,j-1,k)) * inv2dy;
                const amrex::Real gpz = (phi_nbr(i,j,k+1) - phi_nbr(i,j,k-1)) * inv2dz;

                const amrex::Real mag     = std::sqrt(gpx*gpx + gpy*gpy + gpz*gpz);
                const amrex::Real inv_mag = 1.0 / amrex::max(mag, nhat_reg);
                const amrex::Real coeff   = phi * (1.0 - phi) / eps_phi;

                c_arrs[nbx](i, j, k, 0) = coeff * gpx * inv_mag;
                c_arrs[nbx](i, j, k, 1) = coeff * gpy * inv_mag;
                c_arrs[nbx](i, j, k, 2) = coeff * gpz * inv_mag;
            });
        amrex::Gpu::synchronize();
    }

    phi_comp.FillBoundary(geom.periodicity());
    fill_neumann_bc(phi_comp, geom);

    // ------------------------------------------------------------------
    // Step 2+3: forward-Euler update
    // Gamma = gamma_coef * |u|_max  (Chiu & Lin §2.1: "c̄ = |u_max|")
    // Use the global maximum liquid-cell speed so the reinitialization term
    // stays active everywhere — using the local speed would switch it off at
    // the interface where LBM speeds drop toward zero due to bounce-back.
    // ------------------------------------------------------------------
    amrex::Real umag_max = 0.0;
    {
        auto const& md_r     = m_macrodata[lev].const_arrays();
        auto const& isf_r    = m_is_fluid[lev].const_arrays();
        umag_max = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real>{},
            m_macrodata[lev], amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real> {
                if (isf_r[nbx](i, j, k, lbm::constants::IS_FLUID_IDX) != 1)
                    return {amrex::Real(0.0)};
                const amrex::Real ux = md_r[nbx](i, j, k, constants::VELX_IDX);
                const amrex::Real uy = md_r[nbx](i, j, k, constants::VELY_IDX);
                const amrex::Real uz = md_r[nbx](i, j, k, constants::VELZ_IDX);
                return {std::sqrt(ux*ux + uy*uy + uz*uz)};
            });
        amrex::ParallelDescriptor::ReduceRealMax(umag_max);
        // Floor: if u=0 everywhere (t=0), use mesh speed so interface is stable.
        umag_max = amrex::max(umag_max, m_mesh_speed * amrex::Real(1.0e-3));
    }

    amrex::MultiFab phi_new(
        m_is_fluid_fraction[lev].boxArray(),
        m_is_fluid_fraction[lev].DistributionMap(),
        1, 0);

    {
        auto const& c_arrs    = phi_comp.const_arrays();
        auto const& pn_arrs   = phi_new.arrays();
        auto const& stat_arrs = m_stationary_mask[lev].const_arrays();
        auto const& frac_arrs_ro = m_is_fluid_fraction[lev].const_arrays();

        amrex::ParallelFor(
            phi_new,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Solid cells (wall, baffles, impeller) must not be updated by the
                // phase-field PDE — their Phi is permanently 0 (solid) and the
                // large gradient at their fluid-facing faces is NOT a free surface.
                // Updating them causes the interface to blow up at every wall boundary.
                if (stat_arrs[nbx](i, j, k) == 0) {
                    pn_arrs[nbx](i, j, k, 0) = 0.0;
                    return;
                }
                // Also skip cells that are inside the moving body (phi already = 0)
                const amrex::Real phi_cur = frac_arrs_ro[nbx](i, j, k, 0);
                if (phi_cur <= 0.0) {
                    pn_arrs[nbx](i, j, k, 0) = 0.0;
                    return;
                }

                const amrex::Real phi = phi_cur;

                // Conservative 1st-order upwind advection  -div(u*Phi)
                // Uses face-centred velocities interpolated as cell averages.
                // This form is exact regardless of whether div(u)=0, so it
                // remains correct when bubble forces / reaction sources make
                // the velocity field locally non-solenoidal.
                //
                // Face flux (Godunov upwind):
                //   F_{i+1/2} = max(ux_{i+1/2}, 0)*Phi_i
                //              + min(ux_{i+1/2}, 0)*Phi_{i+1}
                //   ux_{i+1/2} = 0.5*(ux_i + ux_{i+1})
                //   div(u*Phi) = (F_{i+1/2} - F_{i-1/2}) / dx  + ...
                // Neumann BC at solid faces for phi and velocity:
                // - phi: use current cell's phi (zero gradient across solid face)
                // - velocity: use 0 (no-slip at solid wall)
                auto phi_nbr2 = [&](int ii, int jj, int kk) -> amrex::Real {
                    return (stat_arrs[nbx](ii, jj, kk) == 0) ? phi
                           : phi_arrs[nbx](ii, jj, kk, 0);
                };
                auto vel_nbr = [&](int ii, int jj, int kk, int comp) -> amrex::Real {
                    return (stat_arrs[nbx](ii, jj, kk) == 0) ? amrex::Real(0.0)
                           : md_arrs[nbx](ii, jj, kk, comp);
                };

                const amrex::Real ux  = md_arrs[nbx](i, j, k, constants::VELX_IDX);
                const amrex::Real uy  = md_arrs[nbx](i, j, k, constants::VELY_IDX);
                const amrex::Real uz  = md_arrs[nbx](i, j, k, constants::VELZ_IDX);
                const amrex::Real uxm = vel_nbr(i-1, j, k, constants::VELX_IDX);
                const amrex::Real uxp = vel_nbr(i+1, j, k, constants::VELX_IDX);
                const amrex::Real uym = vel_nbr(i, j-1, k, constants::VELY_IDX);
                const amrex::Real uyp = vel_nbr(i, j+1, k, constants::VELY_IDX);
                const amrex::Real uzm = vel_nbr(i, j, k-1, constants::VELZ_IDX);
                const amrex::Real uzp = vel_nbr(i, j, k+1, constants::VELZ_IDX);

                const amrex::Real phim  = phi_nbr2(i-1, j, k);
                const amrex::Real phip  = phi_nbr2(i+1, j, k);
                const amrex::Real phijm = phi_nbr2(i, j-1, k);
                const amrex::Real phijp = phi_nbr2(i, j+1, k);
                const amrex::Real phikm = phi_nbr2(i, j, k-1);
                const amrex::Real phikp = phi_nbr2(i, j, k+1);

                // Face-centre velocities (arithmetic average)
                const amrex::Real ux_hi = 0.5 * (ux + uxp);
                const amrex::Real ux_lo = 0.5 * (uxm + ux);
                const amrex::Real uy_hi = 0.5 * (uy + uyp);
                const amrex::Real uy_lo = 0.5 * (uym + uy);
                const amrex::Real uz_hi = 0.5 * (uz + uzp);
                const amrex::Real uz_lo = 0.5 * (uzm + uz);

                // Upwind face fluxes
                const amrex::Real Fx_hi = amrex::max(ux_hi, 0.0)*phi  + amrex::min(ux_hi, 0.0)*phip;
                const amrex::Real Fx_lo = amrex::max(ux_lo, 0.0)*phim + amrex::min(ux_lo, 0.0)*phi;
                const amrex::Real Fy_hi = amrex::max(uy_hi, 0.0)*phi  + amrex::min(uy_hi, 0.0)*phijp;
                const amrex::Real Fy_lo = amrex::max(uy_lo, 0.0)*phijm+ amrex::min(uy_lo, 0.0)*phi;
                const amrex::Real Fz_hi = amrex::max(uz_hi, 0.0)*phi  + amrex::min(uz_hi, 0.0)*phikp;
                const amrex::Real Fz_lo = amrex::max(uz_lo, 0.0)*phikm+ amrex::min(uz_lo, 0.0)*phi;

                // Conservative flux divergence
                const amrex::Real adv_x = (Fx_hi - Fx_lo) / dx0;
                const amrex::Real adv_y = (Fy_hi - Fy_lo) / dy0;
                const amrex::Real adv_z = (Fz_hi - Fz_lo) / dz0;

                // Laplacian(Phi) for diffusion term — Neumann at solid faces
                const amrex::Real lap_phi =
                    (phi_nbr2(i+1,j,k) - 2.0*phi + phi_nbr2(i-1,j,k)) * inv_dx2
                  + (phi_nbr2(i,j+1,k) - 2.0*phi + phi_nbr2(i,j-1,k)) * inv_dy2
                  + (phi_nbr2(i,j,k+1) - 2.0*phi + phi_nbr2(i,j,k-1)) * inv_dz2;

                // div(C) for compression term — Neumann at solid faces:
                // solid cells have C=0 (set in Step 1), but using that value
                // for a fluid cell adjacent to a solid face creates a spurious
                // divergence  (0 - C_fluid) / 2dx.  Substitute the current
                // cell's C-vector for any solid neighbor (zero-gradient BC).
                auto c_nbr = [&](int ii, int jj, int kk, int comp) -> amrex::Real {
                    return (stat_arrs[nbx](ii, jj, kk) == 0)
                           ? c_arrs[nbx](i,  j,  k,  comp)
                           : c_arrs[nbx](ii, jj, kk, comp);
                };
                const amrex::Real div_c =
                    (c_nbr(i+1,j,k,0) - c_nbr(i-1,j,k,0)) * inv2dx
                  + (c_nbr(i,j+1,k,1) - c_nbr(i,j-1,k,1)) * inv2dy
                  + (c_nbr(i,j,k+1,2) - c_nbr(i,j,k-1,2)) * inv2dz;

                const amrex::Real gamma = gamma_coef * umag_max;

                const amrex::Real rhs = -(adv_x + adv_y + adv_z)
                                        + gamma * (lap_phi - div_c);

                pn_arrs[nbx](i, j, k, 0) = phi + dt * rhs;
            });
        amrex::Gpu::synchronize();
    }

    // ------------------------------------------------------------------
    // Step 4: mass redistribution (Chiu & Lin Sec. 2.5)
    // ------------------------------------------------------------------
    constexpr amrex::Real phi_lo = 0.001;
    constexpr amrex::Real phi_hi = 0.999;
    const amrex::Real M0 = m_phi_M0;

    // (i) clip
    {
        auto const& pn = phi_new.arrays();
        amrex::ParallelFor(phi_new,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                pn[nbx](i,j,k,0) = amrex::max(0.0, amrex::min(1.0, pn[nbx](i,j,k,0)));
            });
        amrex::Gpu::synchronize();
    }

    // (ii–iii) total mass M and interface cell count N_G
    amrex::Real M_cur = 0.0;
    long        N_G   = 0;
    {
        auto const& pn = phi_new.const_arrays();
        auto [sum_phi, cnt] = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
            amrex::TypeList<amrex::Real, long>{},
            phi_new, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real, long> {
                const amrex::Real p = pn[nbx](i,j,k,0);
                return {p, (p > phi_lo && p < phi_hi) ? 1L : 0L};
            });
        amrex::ParallelDescriptor::ReduceRealSum(sum_phi);
        amrex::ParallelDescriptor::ReduceLongSum(cnt);
        M_cur = sum_phi;
        N_G   = cnt;
    }

    // (iv) distribute residual G = M0 - M over interface cells
    if (N_G > 0) {
        const amrex::Real G_per_cell =
            (M0 - M_cur) / static_cast<amrex::Real>(N_G);
        auto const& pn = phi_new.arrays();
        amrex::ParallelFor(phi_new,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real p = pn[nbx](i,j,k,0);
                if (p > phi_lo && p < phi_hi) {
                    pn[nbx](i,j,k,0) =
                        amrex::max(0.0, amrex::min(1.0, p + G_per_cell));
                }
            });
        amrex::Gpu::synchronize();
    }

    // ------------------------------------------------------------------
    // Step 5: commit and update LBM fluid domain
    // ------------------------------------------------------------------
    amrex::MultiFab::Copy(m_is_fluid_fraction[lev], phi_new, 0, 0, 1, 0);

    // Re-enforce solid cells to Phi=0.  The phase-field PDE skips them (above),
    // but re-apply here for safety after any interpolation / copy operations.
    {
        auto const& frac = m_is_fluid_fraction[lev].arrays();
        auto const& stat = m_stationary_mask[lev].const_arrays();
        amrex::ParallelFor(m_is_fluid_fraction[lev],
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                if (stat[nbx](i, j, k) == 0) {
                    frac[nbx](i, j, k, 0) = 0.0;
                }
            });
        amrex::Gpu::synchronize();
    }

    // ------------------------------------------------------------------
    // Step 6: update IS_FLUID from new Phi, then initialize f/g for
    // cells that transitioned between gas and liquid.
    //
    // This is intentionally separate from the solid-body refill_and_spill
    // in advance().  For the gas-liquid interface:
    //   - newly-gas  (liquid→gas): zero f,g — mass leaves the LBM domain.
    //   - newly-liquid(gas→liquid): set f,g = equilibrium(rho_nbr, u=0, T_nbr)
    //     with rho and T taken from the nearest persistently-liquid neighbor
    //     and u=0.  This avoids importing spurious momentum (Körner 2005).
    // ------------------------------------------------------------------

    // Ghost-fill before IS_FLUID update so boundary cells are coherent.
    m_is_fluid_fraction[lev].FillBoundary(geom.periodicity());
    m_f[lev].FillBoundary(geom.periodicity());
    m_g[lev].FillBoundary(geom.periodicity());
    for (int ci = 0; ci < m_n_components; ++ci) {
        m_component_lattices[ci][lev].FillBoundary(geom.periodicity());
    }

    // Snapshot old IS_FLUID *before* the threshold update.
    amrex::iMultiFab old_is_fluid_fs(
        m_is_fluid[lev].boxArray(), m_is_fluid[lev].DistributionMap(), 1, 1);
    amrex::iMultiFab::Copy(old_is_fluid_fs, m_is_fluid[lev],
                           lbm::constants::IS_FLUID_IDX, 0, 1, 0);
    old_is_fluid_fs.FillBoundary(geom.periodicity());

    // Recompute IS_FLUID mask from new Phi.
    update_is_fluid_from_fraction_and_mark(lev, m_is_fluid_fraction_threshold);

    {
        const stencil::Stencil stencil;
        const auto& evs     = stencil.evs;
        const auto& weights = stencil.weights;
        const amrex::Real theta0 = stencil::Stencil::THETA0;

        auto const& old_arrs  = old_is_fluid_fs.const_arrays();
        auto const& new_arrs  = m_is_fluid[lev].const_arrays();
        auto const& f_arrs    = m_f[lev].arrays();
        auto const& g_arrs    = m_g[lev].arrays();
        auto const& md_arrs   = m_macrodata[lev].const_arrays();

        const amrex::Real l_mesh_speed  = m_mesh_speed;
        const amrex::Real Rg            = m_R_u / m_m_bar;
        const amrex::Real cv            = Rg / (m_adiabaticExponent - 1.0);
        const amrex::Real T_ref         = m_initialTemperature;
        const amrex::Real adiabaticExp  = m_adiabaticExponent;

        amrex::ParallelFor(m_f[lev], amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const int old_if = old_arrs[nbx](i, j, k, 0);
                const int new_if = new_arrs[nbx](i, j, k, lbm::constants::IS_FLUID_IDX);

                // Newly gas: liquid left the domain — zero distributions.
                if (old_if == 1 && new_if == 0) {
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        f_arrs[nbx](i, j, k, q) = 0.0;
                        g_arrs[nbx](i, j, k, q) = 0.0;
                    }
                    return;
                }

                // Newly liquid: initialize to equilibrium at rest.
                if (old_if == 0 && new_if == 1) {
                    const auto& farr = f_arrs[nbx];
                    const auto lo = amrex::lbound(farr);
                    const auto hi = amrex::ubound(farr);

                    // Find rho and T from a persistently-liquid neighbor.
                    amrex::Real rho_init = 0.0;
                    amrex::Real T_init   = T_ref;
                    bool found = false;

                    // First pass: look for a neighbor that was AND still is fluid.
                    for (int nq = 1; nq < constants::N_MICRO_STATES && !found; ++nq) {
                        const int ni = i + evs[nq][0];
                        const int nj = j + evs[nq][1];
                        const int nk = k + evs[nq][2];
                        if (ni < lo.x || ni > hi.x ||
                            nj < lo.y || nj > hi.y ||
                            nk < lo.z || nk > hi.z) continue;
                        if (old_arrs[nbx](ni, nj, nk, 0) == 1 &&
                            new_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                            rho_init = md_arrs[nbx](ni, nj, nk, constants::RHO_IDX);
                            T_init   = md_arrs[nbx](ni, nj, nk, constants::TEMPERATURE_IDX);
                            found    = true;
                        }
                    }
                    // Fallback: any currently-fluid neighbor.
                    for (int nq = 1; nq < constants::N_MICRO_STATES && !found; ++nq) {
                        const int ni = i + evs[nq][0];
                        const int nj = j + evs[nq][1];
                        const int nk = k + evs[nq][2];
                        if (ni < lo.x || ni > hi.x ||
                            nj < lo.y || nj > hi.y ||
                            nk < lo.z || nk > hi.z) continue;
                        if (new_arrs[nbx](ni, nj, nk, lbm::constants::IS_FLUID_IDX) == 1) {
                            rho_init = md_arrs[nbx](ni, nj, nk, constants::RHO_IDX);
                            T_init   = md_arrs[nbx](ni, nj, nk, constants::TEMPERATURE_IDX);
                            found    = true;
                        }
                    }

                    if (!found || rho_init <= 0.0) return; // truly isolated — leave as zero

                    // f = f_eq(rho_init, u=0, T_init), u=0 so pxx=pyy=pzz=Rg*T
                    const amrex::RealVect zero_vel = {AMREX_D_DECL(0.0, 0.0, 0.0)};
                    const amrex::Real cs2 = Rg * T_init;

                    // g = g_eq(two_rho_e, u=0, T_init) using exact IC recipe.
                    const amrex::Real two_rho_e = get_energy(T_init, rho_init,
                        0.0, 0.0, 0.0, cv);

                    amrex::Real rxx_eq(0.0), ryy_eq(0.0), rzz_eq(0.0),
                                rxy_eq(0.0), rxz_eq(0.0), ryz_eq(0.0);
                    amrex::RealVect heat_flux = {AMREX_D_DECL(0.0, 0.0, 0.0)};
                    get_equilibrium_moments(rho_init, zero_vel, two_rho_e, cv,
                        Rg, heat_flux, rxx_eq, ryy_eq, rzz_eq,
                        rxy_eq, rxz_eq, ryz_eq);
                    amrex::GpuArray<amrex::Real, 6> flux_of_hf = {
                        rxx_eq, ryy_eq, rzz_eq, rxy_eq, rxz_eq, ryz_eq};

                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        const auto& ev = evs[q];
                        const amrex::Real wt = weights[q];
                        f_arrs[nbx](i, j, k, q) = set_equilibrium_value(
                            rho_init, zero_vel, cs2, l_mesh_speed, wt, ev);
                        g_arrs[nbx](i, j, k, q) =
                            set_extended_grad_expansion_generic(
                                two_rho_e, heat_flux, flux_of_hf,
                                l_mesh_speed, wt, ev, theta0, zero_vel, 1.0);
                    }
                }
            });
        amrex::Gpu::synchronize();
    }

    // Component lattices at the free surface:
    //   newly-gas   → zero (dissolved species exits with the liquid).
    //   newly-liquid → zero (fresh liquid has no dissolved species yet).
    // Both cases map to the same action: zero on any IS_FLUID state change.
    for (int ci = 0; ci < m_n_components; ++ci) {
        auto const& comp_arrs   = m_component_lattices[ci][lev].arrays();
        auto const& old_arrs_c  = old_is_fluid_fs.const_arrays();
        auto const& new_arrs_c  = m_is_fluid[lev].const_arrays();
        amrex::ParallelFor(m_component_lattices[ci][lev], amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const int old_if = old_arrs_c[nbx](i, j, k, 0);
                const int new_if = new_arrs_c[nbx](i, j, k, lbm::constants::IS_FLUID_IDX);
                if (old_if != new_if) {
                    for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                        comp_arrs[nbx](i, j, k, q) = 0.0;
                    }
                }
            });
        amrex::Gpu::synchronize();
    }

    m_is_fluid[lev].FillBoundary(geom.periodicity());
    m_f[lev].FillBoundary(geom.periodicity());
    m_g[lev].FillBoundary(geom.periodicity());
    for (int ci = 0; ci < m_n_components; ++ci) {
        m_component_lattices[ci][lev].FillBoundary(geom.periodicity());
    }
}
#endif  // advance_phi disabled — FSLBM active

} // namespace lbm

// ============================================================================
// Bubble body force: correct forcing for the entropic thermal D3Q27 model.
//
// Goal: add momentum F*dt to the fluid while leaving density unchanged.
//
// The wrong approach: He-Luo  delta_f = w_q * (e_q.F) / cs^2  assumes cs^2
// = 1/3 (standard isothermal LBM) and modifies only the e_q-linear moment.
//
// The correct approach for this model (cs^2 = Ru/m_bar * T, cell-local):
//   delta_f_q = f_eq(rho, u + delta_u, T) - f_eq(rho, u, T)
//   where delta_u = F * dt / rho
//
// By construction of the equilibrium:
//   sum_q  delta_f_q          = 0              (density unchanged)
//   sum_q  e_qa * delta_f_q   = rho * delta_u  = F * dt  (correct momentum)
//   Higher moments (stress, energy) shift self-consistently.
//
// The extended stress tensor pxx = ux^2 + r_temperature + dt*omega_corr*D_CORR_X
// (where r_temperature = P/rho = (R/m_bar)*T) is recomputed with the shifted
// velocity; the SGS correction term is the same for both states.
// ============================================================================
namespace lbm {

void LBM::apply_bubble_body_force(int lev, const amrex::MultiFab& force_mf)
{
    BL_PROFILE("LBM::apply_bubble_body_force()");

    const stencil::Stencil stencil;
    const auto& evs    = stencil.evs;
    const auto& weight = stencil.weights;

    const amrex::Real l_mesh_speed   = m_mesh_speed;
    const amrex::Real spec_gas_const = m_R_u / m_m_bar;  // (R/m_bar) = P/(rho*T)
    const amrex::Real nu             = m_nu;
    const amrex::Real dt             = m_dts[lev];

    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& force_arrs    = force_mf.const_arrays();
    auto const& f_arrs        = m_f[lev].arrays();
    auto const& md_arrs       = m_macrodata[lev].const_arrays();
    auto const& d_arrs        = m_derived[lev].const_arrays();

    amrex::ParallelFor(
        m_f[lev], amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) != 1) {
                return;
            }

            const amrex::Real Fx = force_arrs[nbx](iv, 0);
            const amrex::Real Fy = force_arrs[nbx](iv, 1);
            const amrex::Real Fz = force_arrs[nbx](iv, 2);
            if (Fx == 0.0 && Fy == 0.0 && Fz == 0.0) { return; }

            const auto md_arr = md_arrs[nbx];
            const auto d_arr  = d_arrs[nbx];

            const amrex::Real rho = md_arr(iv, constants::RHO_IDX);
            if (rho < 1.0e-12) { return; }

            const amrex::Real ux = md_arr(iv, constants::VELX_IDX);
            const amrex::Real uy = md_arr(iv, constants::VELY_IDX);
            const amrex::Real uz = md_arr(iv, constants::VELZ_IDX);

            const amrex::Real temperature = md_arr(iv, constants::TEMPERATURE_IDX);
            // r_temperature = (R/m_bar)*T = P/rho — the isotropic stress entry
            // used in the equilibrium function (pxx = ux^2 + r_temperature, etc.)
            // Note: this is NOT cs^2; the acoustic cs^2 = gamma*(R/m_bar)*T.
            // The momentum equilibrium only sees P/rho, not gamma.
            const amrex::Real r_temperature = spec_gas_const * temperature;

            // SGS correction coefficient (same formula as macrodata_to_equilibrium)
            const amrex::Real omega      = 1.0 / (nu / (r_temperature * dt) + 0.5);
            const amrex::Real omega_corr = (2.0 - omega) / (2.0 * omega * rho);

            // Extended stress tensor at current state
            const amrex::Real pxx_0 = ux*ux + r_temperature
                + dt * omega_corr * d_arr(iv, constants::D_Q_CORR_X_IDX);
            const amrex::Real pyy_0 = uy*uy + r_temperature
                + dt * omega_corr * d_arr(iv, constants::D_Q_CORR_Y_IDX);
            const amrex::Real pzz_0 = uz*uz + r_temperature
                + dt * omega_corr * d_arr(iv, constants::D_Q_CORR_Z_IDX);

            // Shifted velocity: delta_u = F * dt / rho
            const amrex::Real inv_rho = 1.0 / rho;
            const amrex::Real ux1 = ux + Fx * dt * inv_rho;
            const amrex::Real uy1 = uy + Fy * dt * inv_rho;
            const amrex::Real uz1 = uz + Fz * dt * inv_rho;

            // Safety guard: skip if the velocity perturbation is unphysically
            // large (> 0.5 LB/step).  This should never trigger in a correct run;
            // if it does it indicates a unit-conversion bug in the force deposition
            // and is better treated as a NaN-guard than as a silent clamp.
            const amrex::Real du2 = (ux1-ux)*(ux1-ux) + (uy1-uy)*(uy1-uy) + (uz1-uz)*(uz1-uz);
            if (du2 > 0.25) { return; }  // |δu| > 0.5 — skip this cell

            // Extended stress tensor at shifted state
            // (kinematic u^2 changes; r_temperature and SGS D_CORR unchanged)
            const amrex::Real pxx_1 = ux1*ux1 + r_temperature
                + dt * omega_corr * d_arr(iv, constants::D_Q_CORR_X_IDX);
            const amrex::Real pyy_1 = uy1*uy1 + r_temperature
                + dt * omega_corr * d_arr(iv, constants::D_Q_CORR_Y_IDX);
            const amrex::Real pzz_1 = uz1*uz1 + r_temperature
                + dt * omega_corr * d_arr(iv, constants::D_Q_CORR_Z_IDX);

            const amrex::RealVect vel0 = {AMREX_D_DECL(ux,  uy,  uz )};
            const amrex::RealVect vel1 = {AMREX_D_DECL(ux1, uy1, uz1)};

            for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                const auto& ev = evs[q];
                const amrex::Real wt = weight[q];
                const amrex::Real feq0 = set_extended_equilibrium_value(
                    rho, vel0, pxx_0, pyy_0, pzz_0, l_mesh_speed, wt, ev);
                const amrex::Real feq1 = set_extended_equilibrium_value(
                    rho, vel1, pxx_1, pyy_1, pzz_1, l_mesh_speed, wt, ev);
                f_arrs[nbx](iv, q) += feq1 - feq0;
            }
        });
    amrex::Gpu::synchronize();
}

// ============================================================================
// Bubble O2 source term
// Convert mol/(m³·s) → LB_rho/step, distribute using local equilibrium shape.
//
// The density increment d_rho is deposited as:
//   delta_f_q = d_rho * f_eq(rho=1, vel, T)
// consistent with apply_reaction_source_terms.  This adds dissolved O2 mass
// in an equilibrium state advecting at the local fluid velocity, which avoids
// spurious non-equilibrium stress contributions at the next collision step.
// Uniform 1/N_MICRO_STATES weighting would be equivalent only at u=0 and T=T0.
// ============================================================================
void LBM::apply_bubble_o2_source(int lev, const amrex::MultiFab& o2_src_mf)
{
    BL_PROFILE("LBM::apply_bubble_o2_source()");

    if (m_n_components < 1) { return; }

    // Conversion: [mol/(m³·s)] * dt_phys [s] / C_ref [mol/m³ per LB_rho] = [LB_rho/step]
    const amrex::Real conv = m_bubble_params.dt_phys / m_bubble_o2_C_ref;

    const amrex::Real specific_gas_constant = m_R_u / m_m_bar;
    const amrex::Real l_mesh_speed          = m_mesh_speed;

    const stencil::Stencil stencil;
    const auto& evs    = stencil.evs;
    const auto& weight = stencil.weights;

    auto const& is_fluid_arrs = m_is_fluid[lev].const_arrays();
    auto const& src_arrs      = o2_src_mf.const_arrays();
    auto const& md_arrs       = m_macrodata[lev].const_arrays();
    auto const& fO2_arrs      = m_component_lattices[0][lev].arrays();

    amrex::ParallelFor(
        m_component_lattices[0][lev], amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, lbm::constants::IS_FLUID_IDX) != 1) {
                return;
            }

            const amrex::Real d_rho = src_arrs[nbx](iv, 0) * conv;
            if (d_rho == 0.0) { return; }

            // Local fluid velocity and r_temperature = (R/m_bar)*T = P/rho
            const amrex::RealVect vel = {AMREX_D_DECL(
                md_arrs[nbx](iv, constants::VELX_IDX),
                md_arrs[nbx](iv, constants::VELY_IDX),
                md_arrs[nbx](iv, constants::VELZ_IDX))};
            const amrex::Real r_temperature =
                specific_gas_constant * md_arrs[nbx](iv, constants::TEMPERATURE_IDX);

            // Equilibrium stress entries (no viscous correction — pure kinematic)
            const amrex::Real pxx_eq = vel[0]*vel[0] + r_temperature;
            const amrex::Real pyy_eq = vel[1]*vel[1] + r_temperature;
            const amrex::Real pzz_eq = AMREX_D_PICK(0.0, 0.0, vel[2]*vel[2] + r_temperature);

            for (int q = 0; q < constants::N_MICRO_STATES; ++q) {
                const amrex::Real f_eq_unit = set_extended_equilibrium_value(
                    1.0, vel, pxx_eq, pyy_eq, pzz_eq, l_mesh_speed, weight[q], evs[q]);
                fO2_arrs[nbx](iv, q) += d_rho * f_eq_unit;
            }
        });
    amrex::Gpu::synchronize();
}

// ============================================================================
// FSLBM helper: derive IS_FLUID_IDX and all boundary markers directly from
// m_cell_type, bypassing m_is_fluid_fraction entirely.
//
// This is the FSLBM replacement for update_is_fluid_from_fraction_and_mark.
// We do NOT use m_is_fluid_fraction because:
//   - reconstruct_body_sdf writes impeller SDF values into it each step
//   - refill_and_spill applies threshold=0.5 which would clobber FSLBM
//     interface cells (phi<0.5 -> IS_FLUID=0 -> body_sync Case A -> CELL_SOLID)
// Instead, derive IS_FLUID directly from m_cell_type:
//   CELL_LIQUID or CELL_INTERFACE -> IS_FLUID=1
//   CELL_GAS   or CELL_SOLID      -> IS_FLUID=0
// ============================================================================
void LBM::fslbm_sync_isfluid_markers(const int lev)
{
    BL_PROFILE("LBM::fslbm_sync_isfluid_markers()");
    using namespace lbm::constants;

    // Step 1: derive IS_FLUID_IDX from m_cell_type
    {
        auto const& ct_arrs  = m_cell_type[lev].const_arrays();
        auto const& isf_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const int ct = ct_arrs[nbx](i, j, k, 0);
                isf_arrs[nbx](i, j, k, IS_FLUID_IDX) =
                    (ct == CELL_LIQUID || ct == CELL_INTERFACE) ? 1 : 0;
            });
        amrex::Gpu::synchronize();
    }

    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());

    // Steps 2-4: recompute EB_BOUNDARY, IS_FLUID_SIDE, IS_FLUID_SIDE_BOUNDARY
    // — identical logic to update_is_fluid_from_fraction_and_mark.
    {
        auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto if_arr = is_fluid_arrs[nbx];
                bool all_covered = true;
                const amrex::IntVect nn(1);
                for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
                    const auto dimvec = amrex::IntVect::TheDimensionVector(idir);
                    for (int n = 1; n <= nn[idir]; n++) {
                        all_covered &=
                            (if_arr(iv - n*dimvec, IS_FLUID_IDX) == 0) &&
                            (if_arr(iv + n*dimvec, IS_FLUID_IDX) == 0);
                    }
                }
                if (all_covered || if_arr(iv, IS_FLUID_IDX) == 1)
                    if_arr(iv, EB_BOUNDARY_IDX) = 0;
                else
                    if_arr(iv, EB_BOUNDARY_IDX) = 1;
            });
        amrex::Gpu::synchronize();
    }
    // IS_FLUID_SIDE: fluid cells adjacent to MOVING SOLID (impeller) only.
    //
    // We check m_cell_type == CELL_SOLID AND stat_mask == 1 (moving body in
    // stationary-mask convention: 1=not-stationary = moving or open fluid).
    //
    // Exclusions:
    //   - CELL_GAS (free surface): IS_FLUID=0 from fslbm_sync, BUT
    //     stat_mask=1 (default).  If GAS were treated as a solid here,
    //     f_to_macrodata would apply the impeller velocity to surface cells.
    //   - CELL_SOLID from baffle: stat_mask=0.  Baffle-adjacent cells should
    //     NOT receive body velocity; the baffle is stationary.
    //   - CELL_SOLID from impeller: stat_mask=1 → IS_FLUID_SIDE=1.  These
    //     cells genuinely need the no-slip impeller BC.  ✓
    {
        const stencil::Stencil stencil_s;
        const auto& evs_s = stencil_s.evs;
        auto const& ct_arrs2        = m_cell_type[lev].const_arrays();
        auto const& stat_mask_arrs2 = m_stationary_mask[lev].const_arrays();
        auto const& is_fluid_arrs2  = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto if_arr   = is_fluid_arrs2[nbx];
                const auto ct_arr   = ct_arrs2[nbx];
                const auto sm_arr   = stat_mask_arrs2[nbx];
                if (if_arr(iv, IS_FLUID_IDX) == 0) {
                    if_arr(iv, IS_FLUID_SIDE_IDX) = 0;
                    return;
                }
                // IS_FLUID_SIDE=1 only if adjacent to MOVING solid (impeller).
                //   CELL_SOLID + stat_mask==1  =>  moving solid (impeller)
                //   CELL_SOLID + stat_mask==0  =>  stationary solid (baffle) — excluded
                //   CELL_GAS   (any stat_mask) =>  gas, not a wall — excluded
                bool sees_moving_solid = false;
                for (int q = 0; q < N_MICRO_STATES; q++) {
                    const auto ivn = iv - evs_s[q];
                    if (ct_arr(ivn, 0) == CELL_SOLID && sm_arr(ivn) == 1) {
                        sees_moving_solid = true; break;
                    }
                }
                if_arr(iv, IS_FLUID_SIDE_IDX) = sees_moving_solid ? 1 : 0;
            });
        amrex::Gpu::synchronize();
    }
    {
        const stencil::Stencil stencil;
        const auto& evs = stencil.evs;
        auto const& is_fluid_arrs = m_is_fluid[lev].arrays();
        amrex::ParallelFor(
            m_is_fluid[lev], m_is_fluid[lev].nGrowVect() - 1,
            [=] AMREX_GPU_DEVICE(
                int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto if_arr = is_fluid_arrs[nbx];
                bool sees_side = false;
                for (int q = 0; q < N_MICRO_STATES; ++q)
                    if (if_arr(iv - evs[q], IS_FLUID_SIDE_IDX) == 1) {
                        sees_side = true; break;
                    }
                if (if_arr(iv, IS_FLUID_IDX) == 1 &&
                    if_arr(iv, IS_FLUID_SIDE_IDX) == 0 && sees_side)
                    if_arr(iv, IS_FLUID_SIDE_BOUNDARY_IDX) = 1;
                else
                    if_arr(iv, IS_FLUID_SIDE_BOUNDARY_IDX) = 0;
            });
        amrex::Gpu::synchronize();
    }
    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());
}

// ============================================================================
// FSLBM: Initialize cell type field and fill level (φ) from a sharp interface
// at z = m_free_surface_z (physical units).
//
// Cell classifications (Körner 2005):
//   CELL_SOLID     : IS_FLUID_IDX == 0 (EB or moving body)
//   CELL_GAS       : z_cell > z_surf + 0.5*dz
//   CELL_LIQUID    : z_cell < z_surf - 0.5*dz
//   CELL_INTERFACE : |z_cell - z_surf| <= 0.5*dz   (one-cell-thick band)
//
// m_is_fluid_fraction stores φ (0 for gas, 1 for liquid, linear for interface).
// update_is_fluid_from_fraction_and_mark syncs the integer mask so the gas
// headspace is treated as solid by the LBM flux loops.
// ============================================================================
void LBM::fslbm_init_cell_type(const int lev)
{
    BL_PROFILE("LBM::fslbm_init_cell_type()");

    const auto& geom_l  = Geom(lev);
    const auto  dx      = geom_l.CellSizeArray();
    const auto  plo     = geom_l.ProbLoArray();
    const amrex::Real z_surf = m_free_surface_z;
    const amrex::Real dz     = dx[2];

    using namespace lbm::constants;
    const amrex::Real l_phi_lo = FSLBM_PHI_LO;
    const amrex::Real l_phi_hi = FSLBM_PHI_HI;

    auto const& if_arrs   = m_is_fluid[lev].const_arrays();
    auto const& phi_arrs  = m_phi_fslbm[lev].arrays();
    auto const& ct_arrs   = m_cell_type[lev].arrays();

    amrex::ParallelFor(
        m_cell_type[lev], m_cell_type[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // EB / moving-body cell — no PDFs
            if (if_arrs[nbx](i, j, k, IS_FLUID_IDX) == 0) {
                ct_arrs[nbx](i, j, k, 0) = CELL_SOLID;
                phi_arrs[nbx](i, j, k, 0) = amrex::Real(0.0);
                return;
            }
            const amrex::Real z_cell = plo[2] + (k + amrex::Real(0.5)) * dz;
            // Use a ±1·dz interface band so that z_surf falling exactly on a cell
            // face still produces interface cells with φ ∈ (PHI_LO, PHI_HI).
            // Linear interpolation: φ = (z_surf − z_cell)/(2·dz) + 0.5
            //   z_cell = z_surf − dz  →  φ = 0.5 + 0.5 = 1.0  (clamped to PHI_HI)
            //   z_cell = z_surf       →  φ = 0.5
            //   z_cell = z_surf + dz  →  φ = 0.5 − 0.5 = 0.0  (clamped to PHI_LO)
            if (z_cell >= z_surf + dz) {
                ct_arrs[nbx](i, j, k, 0)  = CELL_GAS;
                phi_arrs[nbx](i, j, k, 0) = amrex::Real(0.0);
            } else if (z_cell <= z_surf - dz) {
                ct_arrs[nbx](i, j, k, 0)  = CELL_LIQUID;
                phi_arrs[nbx](i, j, k, 0) = amrex::Real(1.0);
            } else {
                // Two-cell interface band — linear fill, clamped to (PHI_LO, PHI_HI)
                // so neither boundary cell is immediately converted in Step 5.
                ct_arrs[nbx](i, j, k, 0) = CELL_INTERFACE;
                const amrex::Real phi_lin =
                    (z_surf - z_cell) / (amrex::Real(2.0) * dz) + amrex::Real(0.5);
                phi_arrs[nbx](i, j, k, 0) =
                    amrex::max(l_phi_lo, amrex::min(l_phi_hi, phi_lin));
            }
        });
    amrex::Gpu::synchronize();

    m_phi_fslbm[lev].FillBoundary(geom_l.periodicity());

    // Derive IS_FLUID markers directly from m_cell_type (never from m_phi_fslbm /
    // m_is_fluid_fraction, which may hold moving-body SDF values after
    // reconstruct_body_sdf runs).  This avoids the 0.5 threshold clobbering
    // FSLBM interface cells.
    fslbm_sync_isfluid_markers(lev);

    m_cell_type[lev].FillBoundary(geom_l.periodicity());

    amrex::Print() << "FSLBM cell types initialized at lev=" << lev
                   << "  z_surf=" << z_surf << " m\n";
}

// ============================================================================
// FSLBM: One free-surface time step (Körner 2005, Schwarzmeier 2023 JCP).
//
// This function REPLACES both advance_phi(lev) and stream(lev, m_f).
//
// Algorithm (Phase 1 — no excess-mass redistribution or spawning):
//  Step 1: Push streaming of m_f with Anti-Bounce-Back (ABB) at interface-gas
//          boundaries, yielding f_star.
//  Step 2: Compute per-interface-cell mass flux Δm (pull scheme, pre-stream f).
//  Step 3: Copy f_star -> m_f and FillBoundary.
//  Step 4: Update fill level phi <- phi + Δm/rho, clamp to [0,1].
//  Step 5: Convert cells: phi<1e-4 -> GAS (zero f), phi>1-1e-4 -> LIQUID.
//  Step 6: Sync integer is_fluid mask and FillBoundary everything.
// ============================================================================
void LBM::fslbm_advance_surface(const int lev)
{
    BL_PROFILE("LBM::fslbm_advance_surface()");

    using namespace lbm::constants;

    const stencil::Stencil stencil;
    const auto& evs          = stencil.evs;
    const auto& bounce_dirs  = stencil.bounce_dirs;
    const auto& weights      = stencil.weights;
    const amrex::Real l_mesh_speed    = m_mesh_speed;
    const amrex::Real l_R_u           = m_R_u;
    const amrex::Real l_m_bar         = m_m_bar;
    // Reference density from ic_constant.density — used for ABB ρ₀, seeding
    // targets (feq), and repair threshold so all FSLBM numerics scale with
    // the actual initial bulk density rather than being hard-coded to 1.
    const amrex::Real l_fslbm_rho_ref = m_fslbm_rho_ref;
    // Reference pressure-tensor diagonal for seeding: pdiag = (R/m_bar) * T_ref.
    // This is what collide() uses as p_by_rho = spec_gas_const * temperature,
    // so seeded cells produce feq consistent with the rest of the solver.
    // Using mesh_speed² instead would seed at T = gamma * T_ref (~30 × T_ref
    // for this case), causing a transient temperature spike in newly activated cells.
    const amrex::Real l_fslbm_pdiag_ref = (m_R_u / m_m_bar) * m_initialTemperature;
    // Thermal constants needed to build the correct g equilibrium when seeding
    // newly-activated / repaired cells.  g_eq is NOT the same as f_eq — its
    // zeroth moment is 2ρe = 2ρ·cv·T, not ρ.
    const amrex::Real l_Rg     = m_R_u / m_m_bar;
    const amrex::Real l_cv     = l_Rg / (m_adiabaticExponent - amrex::Real(1.0));
    const amrex::Real l_T_ref  = m_initialTemperature;
    const amrex::Real l_theta0 = stencil::Stencil::THETA0;
    // Surface tension + Laplace-pressure density correction for ABB BC.
    // Δρ_G = -2*sigma*kappa / (Rg * T_interface); when sigma=0 this is zero.
    const amrex::Real l_sigma              = m_fslbm_sigma;
    // Contact angle θ_W: cos(θ) used as ghost-phi modifier for solid neighbors.
    // φ_ghost = φ_fluid + cos(θ) * |∇_tangential φ|
    // θ=90° → cos=0 → φ_ghost = φ_fluid (neutral wetting, original Körner).
    const amrex::Real l_cos_contact_angle  =
        std::cos(m_fslbm_contact_angle_deg * amrex::Real(M_PI) / amrex::Real(180.0));

    // -----------------------------------------------------------------------
    // Body-motion sync: reconcile m_cell_type / m_is_fluid_fraction with the
    // IS_FLUID_IDX written by reconstruct_body_sdf + refill_and_spill this step.
    //
    //  Case A — body swept INTO a cell (ct=LIQUID/INTERFACE, IS_FLUID=0):
    //    Reclassify → CELL_SOLID, φ = 0.
    //    Without this, the cell would still stream PDFs and Step 6's low threshold
    //    (5e-5) would re-activate it as fluid, overriding the body's IS_FLUID=0.
    //
    //  Case B — body swept OUT of a cell (ct=SOLID, IS_FLUID=1):
    //    Reclassify by averaging φ of neighboring non-solid cells:
    //      avg_phi > 0.5  →  CELL_LIQUID (φ=1); f/g already set by refill_and_spill.
    //      avg_phi ≤ 0.5  →  CELL_GAS   (φ=0); zero f/g (gas cells have f=0).
    //    This tracks the actual deformed interface, not the initial flat z_surf.
    //    Step 0 repair below catches any CELL_LIQUID cell left at low rho.
    // -----------------------------------------------------------------------
    {
        const auto& geom_sync = Geom(lev);

        // Ensure phi ghost cells are up-to-date before reading neighbor values.
        m_phi_fslbm[lev].FillBoundary(geom_sync.periodicity());
        m_cell_type[lev].FillBoundary(geom_sync.periodicity());

        auto const& ct_s   = m_cell_type[lev].arrays();
        auto const& phi_s  = m_phi_fslbm[lev].arrays();
        auto const& isf_s  = m_is_fluid[lev].const_arrays();
        auto const& f_s    = m_f[lev].arrays();
        auto const& g_s    = m_g[lev].arrays();
        const auto& l_evs_ab  = evs;

        amrex::ParallelFor(
            m_cell_type[lev], m_cell_type[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const int ct  = ct_s[nbx](i, j, k, 0);
                const int isf = isf_s[nbx](i, j, k, IS_FLUID_IDX);

                if ((ct == CELL_LIQUID || ct == CELL_INTERFACE) && isf == 0) {
                    // Case A: body has swept into this cell this step.
                    ct_s[nbx](i, j, k, 0)  = CELL_SOLID;
                    phi_s[nbx](i, j, k, 0) = amrex::Real(0.0);

                } else if (ct == CELL_SOLID && isf == 1) {
                    // Case B: body has swept out of this cell this step.
                    // Classify as CELL_LIQUID or CELL_GAS by averaging phi of
                    // neighboring non-solid cells.  This tracks the actual
                    // (possibly deformed) interface position rather than the
                    // initial flat z_surf threshold.
                    //   avg_phi > 0.5  →  cell is inside the liquid bulk  → CELL_LIQUID
                    //   avg_phi <= 0.5 →  cell is in / above the gas region → CELL_GAS
                    // Fallback: no non-solid neighbors found → CELL_LIQUID (conservative:
                    // Step 0 repair handles low-rho; FSLBM will re-classify next step).
                    const auto& phi_arr = phi_s[nbx];
                    const auto lo = amrex::lbound(phi_arr);
                    const auto hi = amrex::ubound(phi_arr);
                    amrex::Real phi_sum = amrex::Real(0.0);
                    int phi_count = 0;
                    for (int nq = 1; nq < N_MICRO_STATES; ++nq) {
                        const int ni = i + l_evs_ab[nq][0];
                        const int nj = j + l_evs_ab[nq][1];
                        const int nk = k + l_evs_ab[nq][2];
                        if (ni < lo.x || ni > hi.x ||
                            nj < lo.y || nj > hi.y ||
                            nk < lo.z || nk > hi.z) { continue; }
                        const int nct = ct_s[nbx](ni, nj, nk, 0);
                        if (nct == CELL_SOLID) { continue; }
                        phi_sum += phi_arr(ni, nj, nk, 0);
                        ++phi_count;
                    }
                    const amrex::Real avg_phi = (phi_count > 0)
                        ? (phi_sum / phi_count) : amrex::Real(1.0);
                    const int new_ct = (avg_phi > amrex::Real(0.5))
                        ? CELL_LIQUID : CELL_GAS;
                    ct_s[nbx](i, j, k, 0)  = new_ct;
                    phi_s[nbx](i, j, k, 0) = (new_ct == CELL_LIQUID)
                        ? amrex::Real(1.0) : amrex::Real(0.0);

                    if (new_ct == CELL_GAS) {
                        for (int q = 0; q < N_MICRO_STATES; ++q) {
                            f_s[nbx](i, j, k, q) = amrex::Real(0.0);
                            g_s[nbx](i, j, k, q) = amrex::Real(0.0);
                        }
                    }
                    // CELL_LIQUID: f/g already set by refill_and_spill.
                    // Step 0 repair catches any low-rho residual cases.
                }
            });
        amrex::Gpu::synchronize();
    }
    m_cell_type[lev].FillBoundary(Geom(lev).periodicity());

    // Now that m_cell_type reflects the post-body-motion state, derive IS_FLUID
    // directly from it.  This corrects the IS_FLUID that refill_and_spill set from
    // the SDF (which used threshold=0.5 and would clobber FSLBM interface cells).
    fslbm_sync_isfluid_markers(lev);

    // -----------------------------------------------------------------------
    // Step 0: Repair any CELL_INTERFACE (or CELL_LIQUID) cell whose f
    // populations are unphysical (rho < threshold).  These arise from cells
    // that start as IS_FLUID=0 (tanh SDF < 0.5 near the tank wall) so their
    // f is zero-initialized, then fslbm_sync_isfluid_markers flips IS_FLUID=1
    // without seeding f.  Negative-rho interface cells feed exponentially
    // growing dm and collapse the free surface within 4-5 steps.
    // -----------------------------------------------------------------------
    {
        const amrex::Real rho_repair_threshold = amrex::Real(0.01) * l_fslbm_rho_ref;
        auto const& ct_r   = m_cell_type[lev].const_arrays();
        auto const& f_r    = m_f[lev].arrays();
        auto const& g_r    = m_g[lev].arrays();
        const auto& l_evs_r     = evs;
        const auto& l_weights_r = weights;
        // Use pdiag = (R/m_bar)*T_ref to avoid stale/zero T from macrodata.
        // Consistent with collide's p_by_rho = spec_gas_const * temperature.
        const amrex::Real pdiag_repair = l_fslbm_pdiag_ref;
        amrex::ParallelFor(
            m_f[lev], m_f[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const int ct = ct_r[nbx](i, j, k, 0);
                if (ct != CELL_INTERFACE && ct != CELL_LIQUID) { return; }
                amrex::Real rho = amrex::Real(0.0);
                for (int q = 0; q < N_MICRO_STATES; ++q) rho += f_r[nbx](i, j, k, q);
                if (rho >= rho_repair_threshold) { return; }
                // Seed f to feq(rho_ref, u=0, T_ref) and g to the correct
                // thermal equilibrium g_eq(2ρe_ref, u=0, T_ref).
                const amrex::RealVect zero_vel(AMREX_D_DECL(0, 0, 0));
                amrex::RealVect heat_flux_seed(AMREX_D_DECL(0, 0, 0));
                const amrex::Real two_rho_e_rep = get_energy(
                    l_T_ref, l_fslbm_rho_ref, 0.0, 0.0, 0.0, l_cv);
                amrex::Real rxx_r(0), ryy_r(0), rzz_r(0),
                            rxy_r(0), rxz_r(0), ryz_r(0);
                get_equilibrium_moments(l_fslbm_rho_ref, zero_vel,
                    two_rho_e_rep, l_cv, l_Rg, heat_flux_seed,
                    rxx_r, ryy_r, rzz_r, rxy_r, rxz_r, ryz_r);
                const amrex::GpuArray<amrex::Real, 6> hf_eq_rep = {
                    rxx_r, ryy_r, rzz_r, rxy_r, rxz_r, ryz_r};
                for (int q = 0; q < N_MICRO_STATES; ++q) {
                    const amrex::Real feq_val = set_extended_equilibrium_value(
                        l_fslbm_rho_ref, zero_vel, pdiag_repair, pdiag_repair, pdiag_repair,
                        l_mesh_speed, l_weights_r[q], l_evs_r[q]);
                    f_r[nbx](i, j, k, q) = feq_val;
                    g_r[nbx](i, j, k, q) = set_extended_grad_expansion_generic(
                        two_rho_e_rep, heat_flux_seed, hf_eq_rep,
                        l_mesh_speed, l_weights_r[q], l_evs_r[q],
                        l_theta0, zero_vel, amrex::Real(1.0));
                }
            });
        amrex::Gpu::synchronize();
    }
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    m_g[lev].FillBoundary(Geom(lev).periodicity());

    // -----------------------------------------------------------------------
    // Step 0b: Overshoot repair — reset CELL_INTERFACE cells whose rho
    // exceeds a ceiling, and CELL_LIQUID cells with NaN/Inf (but NOT finite
    // over-density from legitimate spill deposits).
    //
    // WHY CELL_INTERFACE only for finite overshoot:
    //   CELL_LIQUID bulk cells adjacent to the rotating impeller legitimately
    //   receive spill mass from multiple simultaneously-swept blade cells.
    //   For example, 5 cells becoming solid in one step, all naming the same
    //   old-boundary neighbour as spill target, results in rho ≈ 5 in that
    //   cell.  This is CORRECT and TRANSIENT — FSLBM Step 1a streaming will
    //   distribute the excess over 26 neighbours within a few steps.
    //   Clamping CELL_LIQUID bulk cells to rho_ref destroyed 4 mass units per
    //   event, causing a sustained mass-loss cascade and incorrect flow.
    //
    //   CELL_INTERFACE cells at the gas-liquid boundary cannot safely hold
    //   large rho: the ABB formula and mass-flux Step 2 amplify rho → ∞ at
    //   the surface within 2–3 steps if an over-dense interface cell exists.
    //
    // CELL_LIQUID ceiling: NaN/Inf only (any finite value dissipates safely).
    // CELL_INTERFACE ceiling: 5 × rho_ref (protect ABB from amplification).
    // -----------------------------------------------------------------------
    {
        const amrex::Real rho_ceil_ifc = amrex::Real(5.0) * l_fslbm_rho_ref;
        auto const& ct_ob   = m_cell_type[lev].const_arrays();
        auto const& f_ob    = m_f[lev].arrays();
        auto const& g_ob    = m_g[lev].arrays();
        const auto& l_evs_ob     = evs;
        const auto& l_weights_ob = weights;
        const amrex::Real pdiag_ob = l_fslbm_pdiag_ref;
        amrex::ParallelFor(
            m_f[lev], m_f[lev].nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const int ct = ct_ob[nbx](i, j, k, 0);
                if (ct != CELL_INTERFACE && ct != CELL_LIQUID) { return; }
                amrex::Real rho = amrex::Real(0.0);
                for (int q = 0; q < N_MICRO_STATES; ++q) rho += f_ob[nbx](i, j, k, q);
                // For CELL_LIQUID: only clamp NaN/Inf (finite spill deposits
                // should dissipate naturally via streaming, not be clamped).
                // For CELL_INTERFACE: clamp if rho > ceiling (ABB amplification
                // risk).  !(rho <= ceil) also catches NaN.
                const bool is_interface = (ct == CELL_INTERFACE);
                const amrex::Real rho_ceil = is_interface
                    ? rho_ceil_ifc
                    : amrex::Real(1.0e30);  // only catches NaN/Inf for LIQUID
                if (!(rho > rho_ceil)) { return; }
                const amrex::RealVect zero_vel(AMREX_D_DECL(0, 0, 0));
                amrex::RealVect heat_flux_ob(AMREX_D_DECL(0, 0, 0));
                const amrex::Real two_rho_e_ob = get_energy(
                    l_T_ref, l_fslbm_rho_ref, 0.0, 0.0, 0.0, l_cv);
                amrex::Real rxx_ob(0), ryy_ob(0), rzz_ob(0),
                            rxy_ob(0), rxz_ob(0), ryz_ob(0);
                get_equilibrium_moments(l_fslbm_rho_ref, zero_vel,
                    two_rho_e_ob, l_cv, l_Rg, heat_flux_ob,
                    rxx_ob, ryy_ob, rzz_ob, rxy_ob, rxz_ob, ryz_ob);
                const amrex::GpuArray<amrex::Real, 6> hf_eq_ob = {
                    rxx_ob, ryy_ob, rzz_ob, rxy_ob, rxz_ob, ryz_ob};
                for (int q = 0; q < N_MICRO_STATES; ++q) {
                    const amrex::Real feq_val = set_extended_equilibrium_value(
                        l_fslbm_rho_ref, zero_vel, pdiag_ob, pdiag_ob, pdiag_ob,
                        l_mesh_speed, l_weights_ob[q], l_evs_ob[q]);
                    f_ob[nbx](i, j, k, q) = feq_val;
                    g_ob[nbx](i, j, k, q) = set_extended_grad_expansion_generic(
                        two_rho_e_ob, heat_flux_ob, hf_eq_ob,
                        l_mesh_speed, l_weights_ob[q], l_evs_ob[q],
                        l_theta0, zero_vel, amrex::Real(1.0));
                }
            });
        amrex::Gpu::synchronize();
    }
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    m_g[lev].FillBoundary(Geom(lev).periodicity());

    // -----------------------------------------------------------------------
    // Allocate working storage
    // f_star     : post-stream PDFs
    // mass_flux  : mass increment Δm per interface cell
    // -----------------------------------------------------------------------
    amrex::MultiFab f_star(
        boxArray(lev), DistributionMap(lev), N_MICRO_STATES,
        m_f[lev].nGrow(), amrex::MFInfo(), *(m_factory[lev]));
    f_star.setVal(amrex::Real(0.0));

    amrex::MultiFab mass_flux(
        boxArray(lev), DistributionMap(lev), 1,
        m_f[lev].nGrow(), amrex::MFInfo(), *(m_factory[lev]));
    mass_flux.setVal(amrex::Real(0.0));

    // -----------------------------------------------------------------------
    // Step 1a: Push stream from all LIQUID/INTERFACE cells to fluid neighbors.
    //          Solid neighbors → standard bounce-back.
    //          Gas neighbors   → leave f_star slot ZERO (will be filled in 1b).
    //
    //   NO ABB here on purpose.  ABB must be a separate pass to avoid a GPU
    //   data race: both the push from a face-sharing interface neighbor AND the
    //   ABB reflection of the current interface cell want to write to the SAME
    //   f_star slot (the incoming direction coming from the gas side).
    // -----------------------------------------------------------------------
    {
        auto const& fs_w    = f_star.arrays();
        auto const& f_ro    = m_f[lev].const_arrays();
        auto const& ct_arrs = m_cell_type[lev].const_arrays();

        amrex::ParallelFor(
            m_f[lev], m_f[lev].nGrowVect(), N_MICRO_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const int ct_iv = ct_arrs[nbx](iv, 0);
                if (ct_iv != CELL_LIQUID && ct_iv != CELL_INTERFACE) { return; }

                const auto& ev = evs[q];
                const int   bq = bounce_dirs[q];
                const amrex::IntVect ivn(iv + ev);

                const auto& f_arr = f_ro[nbx];
                const auto& lb    = amrex::lbound(f_arr);
                const auto& ub    = amrex::ubound(f_arr);
                const amrex::Box fbox(
                    amrex::IntVect(AMREX_D_DECL(lb.x, lb.y, lb.z)),
                    amrex::IntVect(AMREX_D_DECL(ub.x, ub.y, ub.z)));
                if (!fbox.contains(ivn)) { return; }

                const int ct_ivn = ct_arrs[nbx](ivn, 0);
                const amrex::Real f_q = f_ro[nbx](iv, q);

                if (ct_ivn == CELL_LIQUID || ct_ivn == CELL_INTERFACE) {
                    // Normal push to fluid neighbor
                    fs_w[nbx](ivn, q) = f_q;
                } else if (ct_ivn == CELL_GAS) {
                    // Leave f_star slot zero — Step 1b fills this via pull ABB.
                    // DO NOT bounce-back here: the slot fs_w[iv][bq] will be
                    // written exclusively by Step 1b (no race).
                } else {
                    // SOLID neighbor: standard bounce-back
                    fs_w[nbx](iv, bq) = f_q;
                }
            });
        amrex::Gpu::synchronize();
    }

    // -----------------------------------------------------------------------
    // Pre-Step-1b: Interface normals and curvature for gas-pressure ABB BC.
    //
    // Pass 1 — unit normal n̂ = ∇φ / |∇φ|  (1 ghost layer, FillBoundary'd).
    //   Wall correction (Donath [52]): solid neighbor → substitute current
    //   cell's φ, enforcing a 90° contact angle.  This prevents the interface
    //   normal from being dragged toward the wall when the surface touches it,
    //   which was the cause of instability when the interface hit the tank wall.
    //
    // Pass 2 — curvature  κ = −∇·n̂  (central differences, same wall correction
    //   for n̂ at solid neighbors).  Only computed for CELL_INTERFACE cells.
    //   In LB units Δx=1, so no division by dx is needed.
    // -----------------------------------------------------------------------
    amrex::MultiFab nhat_mf(boxArray(lev), DistributionMap(lev), 3, 1,
                             amrex::MFInfo(), *(m_factory[lev]));
    nhat_mf.setVal(0.0);
    {
        const amrex::Real nhat_reg = 1.0e-8;
        auto const& nh_w    = nhat_mf.arrays();
        auto const& phi_arr = m_phi_fslbm[lev].const_arrays();
        auto const& ct_nh   = m_cell_type[lev].const_arrays();
        amrex::ParallelFor(
            nhat_mf,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                if (ct_nh[nbx](i, j, k, 0) == CELL_SOLID) { return; }
                const amrex::Real phi = phi_arr[nbx](i, j, k, 0);
                // Helper: phi of neighbor using current-cell phi for SOLID (90° base)
                auto pf = [&](int ii, int jj, int kk) -> amrex::Real {
                    return (ct_nh[nbx](ii, jj, kk, 0) == CELL_SOLID)
                        ? phi : phi_arr[nbx](ii, jj, kk, 0);
                };
                // Pre-compute tangential gradients (using 90° ghost for any solid nbr)
                const amrex::Real gx0 = amrex::Real(0.5) * (pf(i+1,j,k) - pf(i-1,j,k));
                const amrex::Real gy0 = amrex::Real(0.5) * (pf(i,j+1,k) - pf(i,j-1,k));
                const amrex::Real gz0 = amrex::Real(0.5) * (pf(i,j,k+1) - pf(i,j,k-1));
                // Wall correction with contact angle θ: for solid in direction (di,dj,dk),
                // φ_ghost = φ + cos(θ) * |∇_tangential φ|.
                // At θ=90°: cos=0 → φ_ghost = φ  (neutral wetting, original Körner).
                auto phi_w = [&](int ii, int jj, int kk) -> amrex::Real {
                    if (ct_nh[nbx](ii, jj, kk, 0) != CELL_SOLID) {
                        return phi_arr[nbx](ii, jj, kk, 0);
                    }
                    int di = ii - i, dj = jj - j, dk = kk - k;
                    amrex::Real gt;
                    if      (di != 0) { gt = std::sqrt(gy0*gy0 + gz0*gz0); }
                    else if (dj != 0) { gt = std::sqrt(gx0*gx0 + gz0*gz0); }
                    else              { gt = std::sqrt(gx0*gx0 + gy0*gy0); }
                    return phi + l_cos_contact_angle * gt;
                };
                const amrex::Real gpx = amrex::Real(0.5) * (phi_w(i+1,j,k) - phi_w(i-1,j,k));
                const amrex::Real gpy = amrex::Real(0.5) * (phi_w(i,j+1,k) - phi_w(i,j-1,k));
                const amrex::Real gpz = amrex::Real(0.5) * (phi_w(i,j,k+1) - phi_w(i,j,k-1));
                const amrex::Real mag     = std::sqrt(gpx*gpx + gpy*gpy + gpz*gpz);
                const amrex::Real inv_mag = amrex::Real(1.0) / amrex::max(mag, nhat_reg);
                nh_w[nbx](i, j, k, 0) = gpx * inv_mag;
                nh_w[nbx](i, j, k, 1) = gpy * inv_mag;
                nh_w[nbx](i, j, k, 2) = gpz * inv_mag;
            });
        amrex::Gpu::synchronize();
    }
    nhat_mf.FillBoundary(Geom(lev).periodicity());

    amrex::MultiFab kappa_mf(boxArray(lev), DistributionMap(lev), 1, 1,
                              amrex::MFInfo(), *(m_factory[lev]));
    kappa_mf.setVal(0.0);
    {
        auto const& kap_w  = kappa_mf.arrays();
        auto const& nh_ro  = nhat_mf.const_arrays();
        auto const& ct_k   = m_cell_type[lev].const_arrays();
        amrex::ParallelFor(
            kappa_mf,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                if (ct_k[nbx](i, j, k, 0) != CELL_INTERFACE) { return; }
                // Wall correction: solid neighbor → use current cell's n̂ component
                auto nx_w = [&](int ii, int jj, int kk) -> amrex::Real {
                    return (ct_k[nbx](ii,jj,kk,0) == CELL_SOLID)
                        ? nh_ro[nbx](i,j,k,0) : nh_ro[nbx](ii,jj,kk,0);
                };
                auto ny_w = [&](int ii, int jj, int kk) -> amrex::Real {
                    return (ct_k[nbx](ii,jj,kk,0) == CELL_SOLID)
                        ? nh_ro[nbx](i,j,k,1) : nh_ro[nbx](ii,jj,kk,1);
                };
                auto nz_w = [&](int ii, int jj, int kk) -> amrex::Real {
                    return (ct_k[nbx](ii,jj,kk,0) == CELL_SOLID)
                        ? nh_ro[nbx](i,j,k,2) : nh_ro[nbx](ii,jj,kk,2);
                };
                const amrex::Real div_n =
                    amrex::Real(0.5) * (nx_w(i+1,j,k) - nx_w(i-1,j,k))
                  + amrex::Real(0.5) * (ny_w(i,j+1,k) - ny_w(i,j-1,k))
                  + amrex::Real(0.5) * (nz_w(i,j,k+1) - nz_w(i,j,k-1));
                kap_w[nbx](i, j, k, 0) = -div_n;  // κ = −∇·n̂
            });
        amrex::Gpu::synchronize();
    }
    kappa_mf.FillBoundary(Geom(lev).periodicity());

    // -----------------------------------------------------------------------
    // Step 1b: Pull ABB for missing-from-gas incoming populations.
    //
    // For each INTERFACE cell iv, for each direction q where the source
    // neighbor iv - ev_q is GAS, the incoming population f_star[iv][q]
    // was left zero by Step 1a (gas cell produced nothing).  Fill it via
    // Anti-Bounce-Back at atmospheric ρ=1 (Körner 2005 Eq. 7):
    //
    //   f_star[iv][q] = f_eq_bq(1, u_iv, T_iv) + f_eq_q(1, u_iv, T_iv)
    //                   - f_pre[iv][bq]
    //
    // where bq is the direction POINTING INTO gas (opposite of q).
    // Each thread writes exactly ONE slot (f_star[iv][q]) and reads only
    // f_pre[iv][bq] (pre-stream, read-only) — no write race possible.
    // -----------------------------------------------------------------------
    {
        auto const& fs_w    = f_star.arrays();
        auto const& f_ro    = m_f[lev].const_arrays();
        auto const& ct_arrs = m_cell_type[lev].const_arrays();
        auto const& md_arrs = m_macrodata[lev].const_arrays();
        auto const& kap_arr = kappa_mf.const_arrays();

        amrex::ParallelFor(
            m_f[lev], m_f[lev].nGrowVect(), N_MICRO_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int q) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                // Only interface cells need ABB
                if (ct_arrs[nbx](iv, 0) != CELL_INTERFACE) { return; }

                // q is the incoming direction (INTO iv).  The source neighbor
                // in direction -ev_q is: iv - ev_q = iv + ev_bq.
                const int   bq  = bounce_dirs[q];       // outgoing direction into gas
                const auto& ev_bq = evs[bq];
                const amrex::IntVect src(iv + ev_bq);   // neighbor in gas direction

                const auto& f_arr = f_ro[nbx];
                const auto& lb    = amrex::lbound(f_arr);
                const auto& ub    = amrex::ubound(f_arr);
                const amrex::Box fbox(
                    amrex::IntVect(AMREX_D_DECL(lb.x, lb.y, lb.z)),
                    amrex::IntVect(AMREX_D_DECL(ub.x, ub.y, ub.z)));
                if (!fbox.contains(src)) { return; }

                // Only apply if the source neighbor is GAS
                if (ct_arrs[nbx](src, 0) != CELL_GAS) { return; }

                const amrex::Real ux = md_arrs[nbx](iv, VELX_IDX);
                const amrex::Real uy = md_arrs[nbx](iv, VELY_IDX);
                const amrex::Real uz = md_arrs[nbx](iv, VELZ_IDX);
                const amrex::RealVect vel(AMREX_D_DECL(ux, uy, uz));
                const amrex::Real T_iv =
                    amrex::max(md_arrs[nbx](iv, TEMPERATURE_IDX), amrex::Real(1.0e-10));

                // Gas-side density for ABB (Schwarzmeier 2023 Eq.11-14,
                // thermal extension: ρ_G = p_G / (R_g T_interface)).
                //
                // Base (old Körner stability clamp): rho_iv = max(Σf_iv, rho_ref).
                // This guarantees feq_bq + feq_q ≥ f_pre[bq] when the impeller
                // drives ρ_iv > ρ_ref near the interface.
                //
                // Laplace correction (σ > 0): Δρ_G = -2σκ / (Rg·T_iv).
                //   κ > 0 (concave toward gas) → Δρ < 0 → lower gas pressure.
                //   κ < 0 (convex toward gas)  → Δρ > 0 → higher gas pressure.
                // When σ = 0: Δρ = 0 → rho_G = rho_iv = EXACT old behaviour,
                // regardless of local temperature T_iv.
                //
                // Wall-corrected κ (computed above) prevents instability when
                // the interface contacts a solid wall (Donath [52]).
                amrex::Real rho_iv = amrex::Real(0.0);
                for (int qq = 0; qq < N_MICRO_STATES; ++qq)
                    rho_iv += f_ro[nbx](iv, qq);
                rho_iv = amrex::max(rho_iv, l_fslbm_rho_ref);
                const amrex::Real kappa_iv = kap_arr[nbx](iv, 0);
                const amrex::Real delta_rho_laplace =
                    -amrex::Real(2.0) * l_sigma * kappa_iv / (l_Rg * T_iv);
                const amrex::Real rho_G = amrex::max(
                    rho_iv + delta_rho_laplace,
                    l_fslbm_rho_ref * amrex::Real(1.0e-3));

                const amrex::Real pxx = ux * ux + l_Rg * T_iv;
                const amrex::Real pyy = uy * uy + l_Rg * T_iv;
                const amrex::Real pzz = uz * uz + l_Rg * T_iv;
                const amrex::Real feq_q  = set_extended_equilibrium_value(
                    rho_G, vel, pxx, pyy, pzz, l_mesh_speed, weights[q],  evs[q]);
                const amrex::Real feq_bq = set_extended_equilibrium_value(
                    rho_G, vel, pxx, pyy, pzz, l_mesh_speed, weights[bq], evs[bq]);
                // Clamp to 0: ABB can in principle produce f<0 for large non-eq stress.
                fs_w[nbx](iv, q) = amrex::max(
                    amrex::Real(0.0), feq_bq + feq_q - f_ro[nbx](iv, bq));
            });
        amrex::Gpu::synchronize();
    }

    // -----------------------------------------------------------------------
    // Step 2: Mass flux per interface cell (pull scheme, pre-stream f).
    //   Δm = Σ_q  S_q * (f_pre[ivn][bq] - f_pre[iv][q])
    //   S_q = 1               if ivn is LIQUID
    //   S_q = 0.5*(phi_iv+phi_ivn) if ivn is INTERFACE
    //   S_q = 0               if ivn is GAS or SOLID
    // -----------------------------------------------------------------------
    // Diagnostic: rho min/max for LIQUID and INTERFACE cells
    {
        amrex::MultiFab rho_diag(boxArray(lev), DistributionMap(lev), 2, 0,
                                 amrex::MFInfo(), *(m_factory[lev]));
        rho_diag.setVal(amrex::Real(0.0));
        {
            auto const& rd   = rho_diag.arrays();
            auto const& f_ro = m_f[lev].const_arrays();
            auto const& ct   = m_cell_type[lev].const_arrays();
            amrex::ParallelFor(rho_diag,
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                    const int ctype = ct[nbx](i,j,k,0);
                    if (ctype != CELL_LIQUID && ctype != CELL_INTERFACE) return;
                    amrex::Real rho = amrex::Real(0.0);
                    for (int q = 0; q < N_MICRO_STATES; ++q) rho += f_ro[nbx](i,j,k,q);
                    if (ctype == CELL_LIQUID)    rd[nbx](i,j,k,0) = rho;
                    if (ctype == CELL_INTERFACE) rd[nbx](i,j,k,1) = rho;
                });
            amrex::Gpu::synchronize();
        }
        if (m_print_int > 0 && m_isteps[0] % m_print_int == 0) {
            amrex::Print() << "FSLBM rho step=" << m_isteps[0]
                           << " liq=[" << rho_diag.min(0) << "," << rho_diag.max(0) << "]"
                           << " ifc=[" << rho_diag.min(1) << "," << rho_diag.max(1) << "]\n";
            if (rho_diag.max(0) > amrex::Real(2.0)) {
                amrex::IntVect mx = rho_diag.maxIndex(0);
                amrex::Print() << "  liq rho_max cell=" << mx << " step=" << m_isteps[0] << "\n";
            }
            if (rho_diag.max(1) > amrex::Real(2.0)) {
                amrex::IntVect mx = rho_diag.maxIndex(1);
                amrex::Print() << "  ifc rho_max cell=" << mx << " step=" << m_isteps[0] << "\n";
            }
        }
    }
    {
        auto const& dm_arrs  = mass_flux.arrays();
        auto const& f_ro     = m_f[lev].const_arrays();
        auto const& ct_arrs  = m_cell_type[lev].const_arrays();
        auto const& phi_arrs = m_phi_fslbm[lev].const_arrays();

        amrex::ParallelFor(
            mass_flux, mass_flux.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (ct_arrs[nbx](iv, 0) != CELL_INTERFACE) { return; }

                const auto& f_arr = f_ro[nbx];
                const auto& lb    = amrex::lbound(f_arr);
                const auto& ub    = amrex::ubound(f_arr);
                const amrex::Box fbox(
                    amrex::IntVect(AMREX_D_DECL(lb.x, lb.y, lb.z)),
                    amrex::IntVect(AMREX_D_DECL(ub.x, ub.y, ub.z)));

                const amrex::Real phi_iv = phi_arrs[nbx](iv, 0);
                amrex::Real dm = amrex::Real(0.0);

                for (int q = 0; q < N_MICRO_STATES; ++q) {
                    const auto& ev = evs[q];
                    const int   bq = bounce_dirs[q];
                    const amrex::IntVect ivn(iv + ev);
                    if (!fbox.contains(ivn)) { continue; }

                    const int ct_ivn = ct_arrs[nbx](ivn, 0);
                    amrex::Real S_q = amrex::Real(0.0);
                    if (ct_ivn == CELL_LIQUID) {
                        S_q = amrex::Real(1.0);
                    } else if (ct_ivn == CELL_INTERFACE) {
                        S_q = amrex::Real(0.5) * (phi_iv + phi_arrs[nbx](ivn, 0));
                    }
                    if (S_q > amrex::Real(0.0)) {
                        dm += S_q * (f_ro[nbx](ivn, bq) - f_ro[nbx](iv, q));
                    }
                }
                dm_arrs[nbx](iv, 0) = dm;
            });
        amrex::Gpu::synchronize();
    }

    // -----------------------------------------------------------------------
    // Step 3: Copy f_star -> m_f and FillBoundary
    // -----------------------------------------------------------------------
    amrex::MultiFab::Copy(
        m_f[lev], f_star, 0, 0, N_MICRO_STATES, m_f[lev].nGrowVect());
    m_f[lev].FillBoundary(Geom(lev).periodicity());

    // -----------------------------------------------------------------------
    // Step 4: phi update  ->  phi^{n+1} = phi^n + Δm / rho_post
    // Use rho computed directly from the post-streaming m_f rather than
    // from m_macrodata.  m_macrodata is updated by collide() AFTER
    // fslbm_advance_surface(), so it holds values from the PREVIOUS step.
    // Cells that were IS_FLUID=0 (f=0) last step but were just repaired /
    // re-activated this step would have m_macrodata[RHO_IDX]≈0, turning
    // any nonzero dm into phi=∞ and collapsing the entire surface.
    // -----------------------------------------------------------------------
    {
        auto const& phi_arrs = m_phi_fslbm[lev].arrays();
        auto const& dm_arrs  = mass_flux.const_arrays();
        auto const& f_ro_cur = m_f[lev].const_arrays();
        auto const& ct_arrs  = m_cell_type[lev].const_arrays();

        amrex::ParallelFor(
            m_phi_fslbm[lev],
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (ct_arrs[nbx](iv, 0) != CELL_INTERFACE) { return; }
                // Compute rho from current post-streaming f (not stale macrodata).
                amrex::Real rho = amrex::Real(0.0);
                for (int q = 0; q < N_MICRO_STATES; ++q)
                    rho += f_ro_cur[nbx](iv, q);
                rho = amrex::max(rho, amrex::Real(1.0e-4) * l_fslbm_rho_ref);
                const amrex::Real phi_new =
                    phi_arrs[nbx](iv, 0) + dm_arrs[nbx](iv, 0) / rho;
                phi_arrs[nbx](iv, 0) =
                    amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), phi_new));
            });
        amrex::Gpu::synchronize();
    }

    //   phi < FSLBM_PHI_LO  ->  CELL_GAS:    zero f
    //   phi > FSLBM_PHI_HI  ->  CELL_LIQUID
    //
    // Simultaneously write a conversion flag into mass_flux (reused as scratch):
    //   +1  =>  this cell just converted to CELL_GAS   (neighbors may need spawning)
    //   -1  =>  this cell just converted to CELL_LIQUID (neighbors may need spawning)
    //    0  =>  no conversion
    // -----------------------------------------------------------------------
    mass_flux.setVal(amrex::Real(0.0));
    {
        auto const& ct_arrs   = m_cell_type[lev].arrays();
        auto const& phi_arrs  = m_phi_fslbm[lev].arrays();
        auto const& f_arrs    = m_f[lev].arrays();
        auto const& flag_arrs = mass_flux.arrays();   // reused as conversion flag

        amrex::ParallelFor(
            m_cell_type[lev],
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (ct_arrs[nbx](iv, 0) != CELL_INTERFACE) { return; }
                const amrex::Real phi = phi_arrs[nbx](iv, 0);
                if (phi < FSLBM_PHI_LO) {
                    ct_arrs[nbx](iv, 0)   = CELL_GAS;
                    phi_arrs[nbx](iv, 0)  = amrex::Real(0.0);
                    flag_arrs[nbx](iv, 0) = amrex::Real(+1.0);   // converted → GAS
                    for (int q = 0; q < N_MICRO_STATES; ++q)
                        f_arrs[nbx](iv, q) = amrex::Real(0.0);
                } else if (phi > FSLBM_PHI_HI) {
                    ct_arrs[nbx](iv, 0)   = CELL_LIQUID;
                    flag_arrs[nbx](iv, 0) = amrex::Real(-1.0);   // converted → LIQUID
                }
            });
        amrex::Gpu::synchronize();
    }

    // Zero component lattices in cells that just converted to CELL_GAS.
    // m_f was already zeroed in the kernel above; component lattices must
    // also be cleared so dissolved scalars don't persist in gas cells and
    // contaminate liquid cells if the interface retreats back over them.
    if (m_n_components > 0) {
        auto const& flag_arrs_ro = mass_flux.const_arrays();
        for (int c = 0; c < m_n_components; ++c) {
            auto const& f_comp_arrs = m_component_lattices[c][lev].arrays();
            amrex::ParallelFor(
                m_component_lattices[c][lev],
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                    if (flag_arrs_ro[nbx](i, j, k, 0) > amrex::Real(0.5)) {
                        for (int q = 0; q < N_MICRO_STATES; ++q)
                            f_comp_arrs[nbx](i, j, k, q) = amrex::Real(0.0);
                    }
                });
            amrex::Gpu::synchronize();
        }
    }

    // -----------------------------------------------------------------------
    // Step 5b: Spawn new interface cells — triggered ONLY by conversions in
    // Step 5 (Körner 2005 §3.3).
    //
    //   Neighbor of a cell that → GAS    AND is LIQUID    →  CELL_INTERFACE, φ = PHI_HI
    //   Neighbor of a cell that → LIQUID AND is GAS       →  CELL_INTERFACE, φ = PHI_LO
    //
    // Using a conversion-flag array means the spawn is idempotent and cannot
    // run on stable interface cells, preventing the spurious mass-loss that
    // occurs when liquid↔gas adjacency is checked unconditionally every step.
    // -----------------------------------------------------------------------
    mass_flux.FillBoundary(Geom(lev).periodicity());
    m_cell_type[lev].FillBoundary(Geom(lev).periodicity());
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    m_g[lev].FillBoundary(Geom(lev).periodicity());
    {
        auto const& ct_arrs   = m_cell_type[lev].arrays();
        auto const& phi_arrs  = m_phi_fslbm[lev].arrays();
        auto const& f_arrs    = m_f[lev].arrays();
        auto const& g_arrs    = m_g[lev].arrays();
        auto const& flag_arrs = mass_flux.const_arrays();

        amrex::ParallelFor(
            m_cell_type[lev],
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const int ct = ct_arrs[nbx](iv, 0);
                if (ct != CELL_LIQUID && ct != CELL_GAS) { return; }

                bool neighbor_converted_to_gas    = false;
                bool neighbor_converted_to_liquid = false;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    const auto ep = amrex::IntVect::TheDimensionVector(d);
                    for (int s : {+1, -1}) {
                        const amrex::Real fl = flag_arrs[nbx](iv + s * ep, 0);
                        if (fl > amrex::Real(0.5))  neighbor_converted_to_gas    = true;
                        if (fl < amrex::Real(-0.5)) neighbor_converted_to_liquid = true;
                    }
                }

                if (ct == CELL_LIQUID && neighbor_converted_to_gas) {
                    // Demote to interface so the surface band is maintained
                    ct_arrs[nbx](iv, 0)  = CELL_INTERFACE;
                    phi_arrs[nbx](iv, 0) = FSLBM_PHI_HI;
                    // PDFs already valid (cell was LIQUID)
                } else if (ct == CELL_GAS && neighbor_converted_to_liquid) {
                    // Promote to interface; seed PDFs from nearest fluid neighbor.
                    // Normalize to rho=1 to prevent inheriting an elevated density
                    // from a high-rho LIQUID neighbor — otherwise the spawned cell
                    // enters the next step's push streaming with rho>1, amplifying
                    // the error each conversion cycle until blow-up.
                    ct_arrs[nbx](iv, 0)  = CELL_INTERFACE;
                    phi_arrs[nbx](iv, 0) = FSLBM_PHI_LO;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        const auto ep = amrex::IntVect::TheDimensionVector(d);
                        for (int s : {+1, -1}) {
                            const auto ivn = iv + s * ep;
                            const int  ctn = ct_arrs[nbx](ivn, 0);
                            if (ctn == CELL_LIQUID || ctn == CELL_INTERFACE) {
                                // Compute neighbor rho = Σf so we can normalize
                                amrex::Real rho_ivn = amrex::Real(0.0);
                                for (int q = 0; q < N_MICRO_STATES; ++q)
                                    rho_ivn += f_arrs[nbx](ivn, q);
                                rho_ivn = amrex::max(rho_ivn, amrex::Real(1.0e-10));
                                // Normalize to l_fslbm_rho_ref (not 1) so the
                                // spawned cell enters streaming at the correct bulk density.
                                const amrex::Real inv_rho = l_fslbm_rho_ref / rho_ivn;
                                for (int q = 0; q < N_MICRO_STATES; ++q) {
                                    f_arrs[nbx](iv, q) = f_arrs[nbx](ivn, q) * inv_rho;
                                    g_arrs[nbx](iv, q) = g_arrs[nbx](ivn, q) * inv_rho;
                                }
                                return;
                            }
                        }
                    }
                }
            });
        amrex::Gpu::synchronize();
    }

    // -----------------------------------------------------------------------
    // Step 6: Sync IS_FLUID markers from updated m_cell_type and FillBoundary.
    fslbm_sync_isfluid_markers(lev);

    m_cell_type[lev].FillBoundary(Geom(lev).periodicity());
    m_phi_fslbm[lev].FillBoundary(Geom(lev).periodicity());
    m_f[lev].FillBoundary(Geom(lev).periodicity());
    m_g[lev].FillBoundary(Geom(lev).periodicity());
}

} // namespace lbm
