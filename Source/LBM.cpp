#include <memory>

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

    m_idata_varnames.push_back("is_fluid");
    m_idata_varnames.push_back("eb_boundary");
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
    m_g.resize(nlevs_max);
    m_eq.resize(nlevs_max);
    m_eq_g.resize(nlevs_max);
    m_derived.resize(nlevs_max);
    m_is_fluid.resize(nlevs_max);
    m_plt_mf.resize(nlevs_max);
    m_mask.resize(nlevs_max);

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

        if (m_chk_int > 0) {
            write_checkpoint_file();
        }

        open_forces_file(true);
        compute_eb_forces();
    } else {
        // restart from a checkpoint
        read_checkpoint_file();

        open_forces_file(false);
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
        pp.query("chk_file", m_chk_file);
        pp.query("chk_int", m_chk_int);
        pp.query("restart", m_restart_chkfile);
        pp.query("file_name_digits", m_file_name_digits);
    }

    {
        amrex::ParmParse pp("lbm");
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

        pp.query("save_streaming", m_save_streaming);
        pp.query("save_derived", m_save_derived);

        pp.query("compute_forces", m_compute_forces);
        pp.query("forces_file", m_forces_file);

        pp.query("initial_temperature", m_initialTemperature);

        pp.query("adiabatic_exponent", m_adiabaticExponent);
        pp.query("mean_molecular_mass", m_m_bar);

        m_speedOfSound_Ref = std::sqrt(
            m_adiabaticExponent * (m_R_u / m_m_bar) * m_initialTemperature);

        m_mesh_speed = m_dx_outer / m_dt_outer;
        m_cs = m_mesh_speed / constants::ROOT3;

        m_cs_2 = m_cs * m_cs;
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

        amrex::Print() << "\n==============================================="
                          "==============================="
                       << std::endl;
        amrex::Print() << "Step: " << step << " dt : " << m_dts[0]
                       << " time: " << cur_time << " to " << cur_time + m_dts[0]
                       << std::endl;

        m_fillpatch_op->fillpatch(0, cur_time, m_f[0]);

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

        m_fillpatch_g_op->fillpatch(lev + 1, m_ts_new[lev + 1], m_g[lev + 1]);

        for (int i = 1; i <= m_nsubsteps[lev + 1]; ++i) {
            m_fillpatch_op->physbc(lev + 1, m_ts_new[lev + 1], m_f[lev + 1]);

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

    stream(lev, m_f);

    stream(lev, m_g);

    if (lev < finest_level) {
        average_down_to(lev, amrex::IntVect(1));
    }

    collide(lev);
}

void LBM::post_time_step()
{
    BL_PROFILE("LBM::post_time_step()");

    for (int lev = 0; lev <= finest_level; ++lev) {
        compute_derived(lev);
    }

    compute_eb_forces();
}

// Stream the information to the neighbor particles
void LBM::stream(const int lev, amrex::Vector<amrex::MultiFab>& fs)
{
    BL_PROFILE("LBM::stream()");

    amrex::MultiFab f_star(
        boxArray(lev), DistributionMap(lev), constants::N_MICRO_STATES,
        fs[lev].nGrow(), amrex::MFInfo(), *(m_factory[lev]));
    f_star.setVal(-1.0);

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
            if (is_fluid_arrs[nbx](iv, 0) == 1) {
                const auto f_arr = f_arrs[nbx];
                const auto fs_arr = fs_arrs[nbx];
                const auto& lb = amrex::lbound(f_arr);
                const auto& ub = amrex::ubound(f_arr);
                const amrex::Box fbox(
                    amrex::IntVect(AMREX_D_DECL(lb.x, lb.y, lb.z)),
                    amrex::IntVect(AMREX_D_DECL(ub.x, ub.y, ub.z)));
                if (fbox.contains(ivn)) {
                    if (is_fluid_arrs[nbx](ivn, 0) != 0) {
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
            if (is_fluid_arrs[nbx](iv, 0) == 1) {

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

    amrex::ParallelFor(
        m_f[lev], m_eq[lev].nGrowVect(), constants::N_MICRO_STATES,
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k),
            int q) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, 0) == 1) {
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

                f_arr(iv, q) += omega * (eq_arr(iv, q) - f_arr(iv, q));

                g_arr(iv, q) += omega * (eq_arr_g(iv, q) - g_arr(iv, q));
            }
        });
    amrex::Gpu::synchronize();
    m_f[lev].FillBoundary(Geom(lev).periodicity());
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

    const stencil::Stencil stencil;
    const auto& evs = stencil.evs;
    amrex::ParallelFor(
        m_macrodata[lev], m_macrodata[lev].nGrowVect(),
        [=] AMREX_GPU_DEVICE(
            int nbx, int i, int j, int AMREX_D_PICK(, /*k*/, k)) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
            if (is_fluid_arrs[nbx](iv, 0) == 1) {

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

            if (if_arr(iv, 0) == 1) {
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

            if (if_arr(iv, 0) == 1) {
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
                if ((is_fluid_arrs[nbx](iv, 1) == 1) &&
                    (mask_arrs[nbx](iv) == 0)) {
                    for (int q = 0; q < constants::N_MICRO_STATES; q++) {
                        const auto& ev = evs[q];
                        const amrex::IntVect ivr(iv + evs[bounce_dirs[q]]);

                        for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
                            fs[idir] += 2.0 * ev[idir] * f_arrs[nbx](ivr, q) *
                                        is_fluid_arrs[nbx](ivr, 0);
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
    m_g[lev].define(
        ba, dm, m_g[lev - 1].nComp(), m_g[lev - 1].nGrow(), amrex::MFInfo(),
        *(m_factory[lev]));
    m_is_fluid[lev].define(
        ba, dm, m_is_fluid[lev - 1].nComp(), m_is_fluid[lev - 1].nGrow());
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

    m_ts_new[lev] = time;
    m_ts_old[lev] = constants::LOW_NUM;

    initialize_is_fluid(lev);
    initialize_mask(lev);
    m_fillpatch_op->fillpatch_from_coarse(lev, time, m_f[lev]);

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
    m_g[lev].define(
        ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
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

    m_ts_new[lev] = time;
    m_ts_old[lev] = constants::LOW_NUM;

    // Initialize the data
    initialize_is_fluid(lev);
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

    fill_f_inside_eb(lev);

    m_f[lev].FillBoundary(Geom(lev).periodicity());
    m_g[lev].FillBoundary(Geom(lev).periodicity());
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
            is_fluid_arrs[nbx](i, j, k, 0) =
                !(flag_arrs[nbx](i, j, k).isRegular() ||
                  flag_arrs[nbx](i, j, k).isSingleValued())
                    ? 0
                    : 1;
        });

    initialize_from_stl(Geom(lev), m_is_fluid[lev]);

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
                    all_covered &= (if_arr(iv - n * dimvec, 0) == 0) &&
                                   (if_arr(iv + n * dimvec, 0) == 0);
                }
            }

            if ((all_covered) || (if_arr(iv, 0) == 1)) {
                if_arr(iv, 1) = 0;
            } else {
                if_arr(iv, 1) = 1;
            }
        });

    m_is_fluid[lev].FillBoundary(Geom(lev).periodicity());
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
            if (is_fluid_arrs[nbx](i, j, k, 0) == 0) {

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
    amrex::MultiFab new_g(
        ba, dm, constants::N_MICRO_STATES, m_f_nghost, amrex::MFInfo(),
        *(m_factory[lev]));

    m_fillpatch_op->fillpatch(lev, time, new_f);
    m_fillpatch_g_op->fillpatch(lev, time, new_g);

    std::swap(new_f, m_f[lev]);
    std::swap(new_g, m_g[lev]);

    m_macrodata[lev].define(
        ba, dm, constants::N_MACRO_STATES, m_macrodata_nghost, amrex::MFInfo(),
        *(m_factory[lev]));
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
    m_plt_mf[lev].clear();
    m_mask[lev].clear();
}

// Set the user defined BC functions
void LBM::set_bcs()
{

    BL_PROFILE("LBM::set_bcs()");
    const bool is_an_energy_lattice(true);

    if (m_velocity_bc_type == "noop") {

        using VelBCOp = bc::BCOpCreator<bc::NoOp>;

        m_fillpatch_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect()), m_f);

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
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect()), m_f);

        m_fillpatch_g_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(
                m_mesh_speed, m_bc_type, m_g[0].nGrowVect(),
                is_an_energy_lattice),
            m_g);

    } else if (m_velocity_bc_type == "channel") {

        using VelBCOp = bc::BCOpCreator<bc::Channel>;

        m_fillpatch_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect()), m_f);

        m_fillpatch_g_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(
                m_mesh_speed, m_bc_type, m_g[0].nGrowVect(),
                is_an_energy_lattice),
            m_g);

    } else if (m_velocity_bc_type == "parabolic") {

        using VelBCOp = bc::BCOpCreator<bc::Parabolic>;

        m_fillpatch_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(m_mesh_speed, m_bc_type, m_f[0].nGrowVect()), m_f);

        m_fillpatch_g_op = std::make_unique<FillPatchOps<VelBCOp>>(
            geom, refRatio(), m_bcs,
            VelBCOp(
                m_mesh_speed, m_bc_type, m_g[0].nGrowVect(),
                is_an_energy_lattice),
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
}

// Check if a field exists
bool LBM::check_field_existence(const std::string& name)
{
    BL_PROFILE("LBM::check_field_existence()");

    {
        const auto vnames = {
            m_macrodata_varnames, m_microdata_varnames, m_microdata_g_varnames,
            m_deriveddata_varnames, m_idata_varnames};

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
        amrex::MultiFab::Copy(*mf, m_derived[lev], srccomp_mid, 0, nc, ngrow);
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
}

void LBM::write_checkpoint_file() const
{
    BL_PROFILE("LBM::write_checkpoint_file()");
    const auto& varnames = m_microdata_varnames;
    const auto& varnames_g = m_microdata_g_varnames;

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
}

void LBM::read_checkpoint_file()
{
    BL_PROFILE("LBM::read_checkpoint_file()");
    const auto& varnames = m_microdata_varnames;
    const auto& varnames_g = m_microdata_g_varnames;

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
} // namespace lbm
