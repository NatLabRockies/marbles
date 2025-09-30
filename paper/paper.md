---
title: 'MARBLES: Multi-scale Adaptively Refined Boltzmann LatticE Solver'
tags:
  - C++
  - adaptive mesh refinement
  - computational fluid dynamics
  - lattice Boltzmann method
authors:
  - firstname: Marc T.
    surname: Henry de Frahan
    orcid: 0000-0001-7742-1565
    affiliation: 1
  - firstname: Nilesh
    surname: Sawant
    orcid: 0000-0001-8403-8943
    affiliation: 1
  - firstname: Ethan
    surname: Young
    orcid: 0000-0003-1106-7406
    affiliation: 1
  - firstname: Hariswaran
    surname: Sitaraman
    orcid: 0000-0001-5304-1664
    affiliation: 1
  - firstname: Ross E.
    surname: Larsen
    orcid: 0000-0002-2928-9835
    affiliation: "1, 2"
  - firstname: Jon
    surname: Rood
    orcid: 0000-0002-7513-3225
    affiliation: 1
affiliations:
  - name: Scalable Algorithms, Modeling and Simulation Group, National Renewable Energy Laboratory, USA
    index: 1
  - name: Renewable and Sustainable Energy Institute, University of Colorado Boulder, Boulder, Colorado 80309
    index: 2
date: 5 October 2023
bibliography: paper.bib
---

# Summary

MARBLES is a computational fluid dynamics solver that leverages the lattice Boltzmann method (LBM) and block-structured adaptive mesh refinement (AMR) to simulate flows in complex media. The solver leverages the AMReX [@AMReX] library, a library that provides underlying data structures, parallel paradigms, grids, domain decomposition, and portability programming models to run on massively parallel computing architectures. All major GPU architectures (e.g., NVIDIA, AMD, Intel) are supported through the use of performance portability functionalities implemented in AMReX. The MARBLES software is released in NREL Software Record SWR-23-37 “MARBLES (Multi-scale Adaptively Refined Boltzmann LatticE Solver)” [@marbles_swr].

MARBLES implements LBM, a mesoscopic approach to computational fluid dynamics in which the computation of the Boltzmann transport equation by advection and collision of fictitious particles results into a solution of the Navier Stokes equation at the macroscale. The evolution equation of LBM proceeds in two parts: a "streaming" step akin to space discretization in which particles representing probability distribution functions at grid cells transfer information about their current state to their nearest neighbors and a "collision" step akin to time discretization in which these advected quantities give rise to new state variables by subsequently relaxing toward their equilibrium distributions over some characteristic time. Macroscopic state variables such as density, momentum and total energy are obtained by locally computing moments of the probability distribution functions. In MARBLES, these steps are carried out on parallelized AMReX blocks in which the phase space is discretized using 27 grid neighbors in 3D (D3Q27) which provides high accuracy and slightly increased memory consumption compared to other LBM stencils. An equivalent stencil is implemented for 2D simulations. The stencil can be sub-divided and applied locally in a region with half the grid spacing using an implementation of the explode and coalesce algorithm [@chen2006grid], where refined cells process two time steps per every step of the parent grid assuming a factor of two grid refinement ratio. This cycling between refinement levels and interpolation between grids is managed by AMReX functionality and minimizes communication between parallel blocks to enable excellent scaling performance.

MARBLES uses the Embedded Boundary (EB) formulation in AMReX to represent complex geometries. In this approach, an arbitrary surface is defined by the user, using either compositions of simple shapes or an STL file. This is used to intersect the Cartesian mesh and define cells that are covered by the geometry (i.e., inside the body). Bounce back conditions are imposed on these surfaces as part of the streaming step to capture the effect of the geometry on the fluid flow. This leads to a very robust handling of very complex geometry, including flow through porous media.

MARBLES is written in C++ and uses an MPI+X approach: the AMR grid patches are distributed on different CPU ranks using MPI. Each grid can be either (i) logically tiled and solved on different threads on a multi-core CPU machine using OpenMP, or (ii) solved on GPU threads on GPU nodes using CUDA, HIP, and SYCL.

# Statement of Need

Several software tools that implement the Lattice Boltzmann method can found online, including a PyTorch based solver, Lettuce [@bedrunka2021lettuce], the commonly used OpenLB [@krause2021openlb], Palabos [@latt2021palabos], waLBerla [@bauer2021walberla; @godenschwager2013framework], and others [@schmieschek2017lb3d; @mora2020concise; @pastewka2019hpc].

In contrast with these solvers, MARBLES is the only code to the authors' knowledge to be built on AMReX data structures which allows us to tackle exascale problems and take advantage of both CPU and GPU hardware with little change to the underlying code. Because of this, MARBLES is easily extensible to other capabilities and has the potential to perform well on emerging architectures. MARBLES, through its use of vendor-agnostic programming models implemented in AMReX, achieves high performance from a small desktop station to the world's largest supercomputer. Additionally, AMReX provides tools to naturally track regions of the fluid in need of refinement (e.g., around an evolving fluid-vapor interface or near the surface of a deforming solid obstacle) and communicate variables between grid refinement levels in a highly performant manner. 

MARBLE's unique features consist in implementing the two-population compressible thermal LBM [@sawant_consistent_2022] in a highly scalable and performance portable framework. The LBM formulation implemented in MARBLES uses one lattice for mass and momentum equations and another lattice for the energy equation, keeping the whole solver within the lattice Boltzmann framework without the need for a hybrid discretization approach. A fully lattice Boltzmann framework combines exact conservation as in finite volume solvers with the flexibility of a structured grid like in immersed boundary finite difference solvers, without their respective drawbacks of grid generation and conservation errors, respectively. Using a nearest-neighbor based correction to equilibrium pressure, the model maintains Galilean invariance up to third order in velocity space, sufficient to faithfully simulate the correct viscosity and thermal diffusivity over a wide range of temperatures, while still using the standard D3Q27 lattice. 

MARBLES is intended for students, researchers and engineers interested in simulating mesoscale and macroscale flows with the Lattice Boltzmann method on modern high performance computing hardware. With performance portability and being vendor agnostic foremost considerations in the design of the code, MARBLES can leverage the computational resources available on the latest heterogeneous exascale platforms, including GPUs. In this context, MARBLES is a valuable tool for studying fluid dynamics and heat transfer scientific and engineering problems involving complex geometries. Its applications include flows in porous media such as electrodes and geological microstructures, aerodynamic calculations such as forces on turbines blades and structures subject to external flows, and non-equilibrium flows such as those involved in semiconductor manufacturing.

# Acknowledgments

This work was authored by the National Renewable Energy Laboratory (NREL), operated under Contract No. DE-AC36-08GO28308. This work was supported by the Laboratory Directed Research and Development (LDRD) Program at NREL. The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.

# References
