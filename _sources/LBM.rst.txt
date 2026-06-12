 .. role:: cpp(code)
    :language: c++

.. _LBM:



Lattice Boltzmann Method
========================

We briefly describe the Lattice Boltzmann (LBM) method. For more details, the reader is referred to the book "The Lattice Boltzmann Method" by Timm Krüger et al (Springer, 2016).

We introduce the distribution function in position-velocity space :math:`f(x, \xi, t)`, where :math:`\xi = \frac{\mathrm{d}x}{\mathrm{d}t}`, the ensemble average of velocities of particles at a location :math:`x`. In this notation, :math:`f` represents the density of particles with average velocity :math:`xi`, at location :math:`x`, at time :math:`t`. This enables a continuum description on a kinetic level, i.e., the macroscopic fluid properties are moments of :math:`f(x, \xi, t)`.

The Boltzmann equation describes the evolution of the distribution function:

.. math::

   \left(\frac{\partial}{\partial t} + \frac{\mathrm{d}x}{\mathrm{d}t} \frac{\partial}{\partial x} + \frac{\mathrm{d}\xi}{\mathrm{d}t} \frac{\partial}{\partial \xi} \right) f(x, \xi, t) = \Omega(f)

where :math:`\Omega(f)` is the collision operator. These collisions conserve mass, momentum and energy. If one uses the Bhatnagar-Gross-Krook (BGK) model with single relaxation time, it can be written:

.. math::

   \Omega(f) = -\frac{1}{\tau} \left( f - f^{\mathrm{eq}} \right)

where :math:`f^{\mathrm{eq}}` is the equilibrium distribution. The next step is to discretize in the velocity space. At each lattice intersection, particles are allowed to move in some specific lattice directions, thereby discretizing the velocity vector space. Because particles do not have equal probabilities of moving in certain directions, it is necessary to normalize the contributions by weights based on the direction of movement. The final step in the formulation is to replace the ensemble average of particles at a given lattice site by a finite number of particles at that site, i.e., :math:`f(x, \xi, t) = f(x, t)`. This last simplification leads to the lattice Boltzmann equation:

.. math::

   f_i(x+e_i\Delta t, t+\Delta t) - f_i(x,t) = -\frac{\Delta t}{\tau} \left( f_i(x,t) - f_i^{\mathrm{eq}}(x,t)\right)

where :math:`i` denotes the components of :math:`f` at a given lattice site. In two dimensions, usually denoted D2Q9, there are nine lattice vectors (i.e., nine components of :math:`f`) connected to each lattice points (the possible directions in which particles can move). The left hand side is referred to as the streaming step and it is a linear operation. It simply states that the distribution at :math:`x` and time :math:`t` and moving at velocity :math:`e_i` will be at :math:`x+e_i` at time :math:`t + \Delta t`. The right hand side is the collision step and is entirely local. These properties of linearity in the streaming and locality in the collisions lead to the computational strengths of LBM. It should be noted that, by using Chapman-Enskog analysis, the lattice-Boltzmann equation reduces to Navier-Stokes equations. The viscosity of the fluid is denoted:

.. math::

   \nu = c_s^2 \left(\tau - \frac{\Delta t}{2} \right)

In this notation, the equilibrium Maxwell distribution along direction :math:`i` is written:

.. math::

   f_i^{\mathrm{eq}} = w_i \rho \left(1+\frac{e_i \cdot u}{c_s^2} + \frac{(e_i \cdot u)^2}{2 c_s^4} - \frac{u \cdot u}{2 c_s^2} \right)

where :math:`u` is the velocity vector, :math:`\rho` is the density, :math:`w_i` are the normalization weights, and :math:`c_s = \frac{\Delta x}{\Delta t} \frac{1}{\sqrt{3}}` is the speed of sound. The macroscopic variables are derived from the distribution function as:

.. math::

   \rho(x, t) &= \sum_{i} f_i(x, t)

   u(x, t) &= \frac{1}{\rho} \sum_{i} e_i f_i(x, t)

MARBLES implements both the D3Q27 LBM approach for two dimensional and three dimensional flow, respectively.

For the thermal version of MARBLES, please refer to additional theory related to the additional energy lattice and cross coupling, 
in the PhD thesis  
"Sawant, Nilesh. 'Kinetic Modeling of Reactive Flows.'. ETH Zurich, 2023. https://doi.org/10.3929/ETHZ-B-000607045."

