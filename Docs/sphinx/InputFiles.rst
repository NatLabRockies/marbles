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

Comprehensive List of Runtime Parameters
----------------------------------------

Below is a categorized list of the available runtime parameters that can be specified in the input file or overriden on the command line.

AMReX Settings
~~~~~~~~~~~~~~

* ``amr.n_cell`` (Array of Ints): Number of cells along each direction at the base AMR level.
* ``amr.max_level`` (Int): Maximum AMR level allowed (0 means no refinement).
* ``amr.max_grid_size`` (Int): Maximum number of cells per box/grid along the x, y, and z directions.
* ``geometry.prob_lo`` (Array of Reals): Coordinates of the domain's lower corner.
* ``geometry.prob_hi`` (Array of Reals): Coordinates of the domain's upper corner.
* ``geometry.is_periodic`` (Array of Ints): Flags denoting periodicity in the x, y, and z directions (1 for true, 0 for false).
* ``amr.max_step`` (Int): Maximum number of time steps (at the base AMR level) to run.
* ``amr.stop_time`` (Real): Maximum physical simulation time to reach.
* ``amr.regrid_int`` (Int): Number of base-level timesteps between AMR regridding operations.
* ``amr.plot_file`` (String): Base name for plot files.
* ``amr.plot_int`` (Int): Number of base-level timesteps between writing plot files.
* ``amr.chk_file`` (String): Base name for checkpoint files (for restarting).
* ``amr.chk_int`` (Int): Number of base-level timesteps between checkpoint files.
* ``amr.restart`` (String): File path to a checkpoint directory to restart the simulation from.
* ``amr.file_name_digits`` (Int): Minimum number of digits used in checkpoint and plot file names.

Lattice Boltzmann (LBM) Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``lbm.dx_outer`` (Real): Lattice spacing (dx) at the coarsest (outer) level.
* ``lbm.dt_outer`` (Real): Time step (dt) at the coarsest (outer) level.
* ``lbm.nu`` (Real): Kinematic viscosity of the fluid.
* ``lbm.alpha`` (Real): Thermal diffusivity of the fluid.
* ``lbm.bc_lo``, ``lbm.bc_hi`` (Array of Ints): Boundary condition types for the lower and upper boundaries in each coordinate direction. Typical values include 0 (periodic), 1 (slip wall), 2 (no-slip wall), 3/4 (velocity inlet/outlet).
* ``lbm.ic_type`` (String): Type of initial condition to apply (e.g., "constant", "ic_taylorgreen").
* ``lbm.velocity_bc_type`` (String): Type of velocity boundary condition to use (e.g., "constant", "channel", "parabolic").
* ``lbm.save_streaming``, ``lbm.save_derived`` (Int): Flags to save post-streaming distributions and derived quantities into plot files (1 = true, 0 = false).
* ``lbm.compute_forces`` (Int): Flag to enable computation of forces on embedded boundaries (1 = true).
* ``lbm.forces_file`` (String): Name of the file to save computed forces.
* ``lbm.initial_temperature`` (Real): Initial global reference temperature.
* ``lbm.body_is_isothermal`` (Int): Flag indicating if the embedded body maintains a constant temperature (1 = true).
* ``lbm.body_temperature`` (Real): Temperature of the embedded isothermal body.
* ``lbm.adiabatic_exponent`` (Real): Adiabatic exponent (gamma) of the fluid.
* ``lbm.mean_molecular_mass`` (Real): Mean molecular mass of the fluid.

Embedded Boundary (EB) Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prefix: ``eb2.*``

* ``eb2.geom_type`` (String): Type of embedded geometry. Supported types include "cylinder", "box", and "stl".
* ``eb2.cylinder_radius`` (Real): Radius of the cylinder.
* ``eb2.cylinder_direction`` (Int): Coordinate axis alignment of the cylinder (0=X, 1=Y, 2=Z).
* ``eb2.cylinder_center`` (Array of Reals): Coordinates of the cylinder's center.
* ``eb2.cylinder_has_fluid_inside`` (Int): Flag to indicate if fluid exists inside the cylinder.
* ``eb2.cylinder_rotation``, ``eb2.cylinder_rotation_axe`` (Real, Int): Parameters to rotate the cylinder.
* ``eb2.box_lo``, ``eb2.box_hi`` (Array of Reals): Lower and upper bounds of a box geometry.
* ``eb2.box_has_fluid_inside`` (Int): Flag to indicate if fluid exists inside the box.
* ``eb2.box_rotation``, ``eb2.box_rotation_axe``, ``eb2.box_rotation_center``: Rotation parameters for a box.
* ``eb2.stl_file`` (String): Path to an STL file for complex boundary definitions.
* ``eb2.stl_scale`` (Real): Scaling factor applied to the STL geometry.
* ``eb2.stl_reverse_normal`` (Int): Reverses the surface normals of the STL object.
* ``eb2.stl_center`` (Array of Reals): Coordinates to center the STL object.

Prefix: ``voxel_cracks.*``

* ``voxel_cracks.use_voxel_cracks`` (Int): Flag to enable voxelized crack geometries.
* ``voxel_cracks.crack_file`` (String): TIFF/Image file path indicating the 3D crack mask.

Velocity Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These parameters depend on the ``lbm.velocity_bc_type`` chosen. Prefixes include ``velocity_bc_constant.*``, ``velocity_bc_channel.*``, ``velocity_bc_parabolic.*``.

* ``dir``, ``normal_dir``, ``tangential_dir`` (Int): Specification of the normal and tangential axes for the boundary velocity profile.
* ``Mach_ref`` (Real): Reference Mach number at the boundary.
* ``initial_density`` (Real): Initial density at the inlet/boundary.
* ``initial_temperature`` (Real): Initial temperature at the inlet/boundary.
* ``adiabatic_exponent`` (Real): Gamma parameter used for the local speed of sound calculation.
* ``mean_molecular_mass`` (Real): Mean molecular mass used for local ideal gas evaluations.

Initial Conditions (IC)
~~~~~~~~~~~~~~~~~~~~~~~

These parameters override global fluid initializations based on the ``lbm.ic_type`` string (e.g., using "constant" translates to the identifier ``ic_constant``). Common prefixes follow ``ic_constant.*``, ``ic_taylorgreen.*``, etc.

* ``density``, ``initial_temperature`` (Real): Overridden initial density and temperature for the domain block.
* ``velocity``, ``mach_components`` (Array of Reals): Component-wise specification of velocity vector or Mach numbers.
* ``adiabatic_exponent``, ``mean_molecular_mass`` (Real): Local thermodynamic properties for the initial condition region.
* ``rho0``, ``v0``, ``omega``, ``wave_length``: Specific scalar parameters for analytical formulations (e.g., for Taylor-Green vortices or tests).

Tagging and AMR Refinement
~~~~~~~~~~~~~~~~~~~~~~~~~~

Prefix: ``tagging.*`` or ``tagging.[indicator_name].*``

* ``tagging.refinement_indicators`` (Array of Strings): List of user-defined names for refinement criteria (e.g., ``yLow``, ``vel_mag``).
* ``max_level`` (Int): The maximum AMR level this specific tagging rule can trigger.
* ``field_name`` (String): The LBM state or derived variable to be checked (e.g., ``density``, ``vel_mag``).
* ``value_greater``, ``value_less`` (Real): Trigger refinement if the variable evaluated is above/below this threshold.
* ``adjacent_difference_greater`` (Real): Trigger refinement if the magnitude difference between neighbor cells exceeds this threshold (often used for discontinuity/shock detection).
* ``in_box_lo``, ``in_box_hi`` (Array of Reals): Bounding box spatial limits within which the refinement tagging rule is physically active.
* ``start_time``, ``end_time`` (Real): Activate the tagging rule only within a certain physical simulation time window.

