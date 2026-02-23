Input Files and Controls
------------------------

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.

This file needs to specified along with the executable as an `argv` option, for example::

  mpirun -np 64 ./marbles inputs


Entries can be overwritten on the command line: ``./marbles inputs amr.plot_int=10``.

Here is an example input file::

  # maximum number of time steps at base AMR level
  max_step = 1000

  # coordinates of domain's lower corner
  geometry.prob_lo = 0.0 0.0 0.0

  # coordinates of domain's upper corner
  geometry.prob_hi = 96.0 16.0 16.0

  # flag for periodicity (here y direction is periodic)
  geometry.is_periodic = 0 1 0

  # number of cells along each direction at base level
  amr.n_cell  = 96 16 16

  # maximum level number allowed
  amr.max_level = 0

  # maximum number of cells per box along x,y,z
  amr.max_grid_size = 16

  # number of timesteps between plot files
  amr.plot_int = 10

  # number of timesteps between checkpoint files
  amr.chk_int = 10

  # LBM parameteris
  lbm.bc_lo = 2 0 1
  lbm.bc_hi = 3 0 1
  lbm.dx_outer = 1.0
  lbm.dt_outer = 1.0
  lbm.nu = 0.01733333333333333
  lbm.save_streaming = 0

  lbm.velocity_bc_type = "channel"
  velocity_bc_channel.initial_density = 1.0
  velocity_bc_channel.Mach_ref = 0.01
  velocity_bc_channel.initial_temperature = 0.03

  lbm.ic_type = "constant"
  ic_constant.density = 1.0
  ic_constant.velocity = 0.0 0.0 0.0

  # embedded boundary
  eb2.geom_type = "cylinder"
  eb2.cylinder_radius = 2.0
  eb2.cylinder_center = 32.0 8.0 8.0
  eb2.cylinder_has_fluid_inside = 0
  eb2.cylinder_height = 256.0
  eb2.cylinder_direction = 1

  # amrex options for trapping FPEs
  amrex.fpe_trap_invalid = 1
  amrex.fpe_trap_zero = 1
  amrex.fpe_trap_overflow = 1
  
Here is an examples of setting up tagging criteria::

  # Tag inside a box and a velocity magnitude value
  tagging.refinement_indicators = yLow vel_mag
  tagging.yLow.in_box_lo = -0.1  -0.52  -0.85
  tagging.yLow.in_box_hi =  3.1 -0.45    0.85
  
  tagging.vel_mag.max_level     = 2
  tagging.vel_mag.value_greater = 1.2e4
  tagging.vel_mag.field_name    = vel_mag

The following keys are implemented: `value_greater`, `value_less`, `adjacent_difference_greater`, `in_box_lo` and `in_box_hi` (to specify a refinement region), `max_level`, `start_time`, and `end_time`. The `field_name` key can be any available variable.
