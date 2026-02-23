.. highlight:: rst


Introduction
============

`MARBLES` is an open-source solver that implements the Lattice Boltzmann method (LBM). It is a high-performance and performance-portable C++ solver that leverages the AMReX library. It can run on GPUs for acceleration and is compatible with GPUs from major manufacturers, such as AMD, NVIDIA, and Intel. The solver can simulate 2D and 3D problems with arbitrarily complex geometry through the use of a level-set function.
Through adaptive mesh refinement (AMR), the solver can dynamically refine and coarsen the mesh in response to the evolving flow features within the simulation. This capability allows for improved resolution in regions of interest while conserving computational resources in less critical areas.

A variety of examples are included to provide a model setup for the various options. These are discussed further in the :ref:`Getting Started<GettingStarted>` section.


Dependencies
------------

`MARBLES` is built on AMReX (available at `https://github.com/AMReX-Codes/amrex <https://github.com/AMReX-Codes/amrex>`_), an adaptive mesh refinement software framework, which provides the underlying software infrastructure for block structured AMR operations. The full AMReX documentation can be found `here <https://amrex-codes.github.io/amrex>`_.


Acknowledgment
--------------

This work was authored by the National Renewable Energy Laboratory (NREL), operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by the Laboratory Directed Research and Development (LDRD) Program at NREL. The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.
