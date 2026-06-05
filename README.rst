MARBLES: Multi-scale Adaptively Refined Boltzmann LatticE Solver
----------------------------------------------------------------

|CI Badge| |Documentation Badge| |License Badge| |AMReX Badge| |C++ Badge|

.. |CI Badge| image:: https://github.com/NREL/marbles/workflows/MARBLES-CI/badge.svg
   :target: https://github.com/NREL/marbles/actions

.. |Documentation Badge| image:: https://github.com/NREL/marbles/workflows/MARBLES-Docs/badge.svg
   :target: https://natlabrockies.github.io/marbles

.. |License Badge| image:: https://img.shields.io/badge/License-Apache%20v2.0-blue.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0

.. |AMReX Badge| image:: https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22
   :target: https://amrex-codes.github.io/amrex/

.. |C++ Badge| image:: https://img.shields.io/badge/language-C%2B%2B17-blue
   :target: https://isocpp.org/

MARBLES (Multi-scale Adaptively Refined Boltzmann LatticE Solver) is a
computational fluid dynamics solver that leverages the lattice
Boltzmann method (LBM) and block-structured adaptive mesh refinement
(AMR) to simulate flows in complex media. MARBLES supports various
parallel decomposition strategies, including the use of the Message
Passing Interface (MPI) and OpenMP threading. All major GPU
architectures (e.g., Intel, AMD, NVIDIA) are supported through the use
of performance portability functionalities implemented in AMReX. The
MARBLES software is released in NREL Software Record `SWR-23-37
<https://doi.org/10.11578/dc.20231009.2>`_ “MARBLES (Multi-scale
Adaptively Refined Boltzmann LatticE Solver)”.

Getting Started
~~~~~~~~~~~~~~~

To compile and run `MARBLES`, one needs a C++ compiler that supports the C++17 standard, and then execute ::

    $ git clone --recursive git@github.com:NREL/marbles.git
    $ cd Build
    $ make realclean && make -j
    $ ./marbles3d.xxx.yyy.ex example.inp

To compile and run `MARBLES` with GPU support, e.g. with CUDA, clone the repository as above, and then execute ::

    $ cd Build
    $ make USE_CUDA=TRUE
    $ ./marbles3d.xxx.yyy.CUDA.ex example.inp    

For additional flags to enable different compilers, architectures, and options, see `AMReX build options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_ . 

.. note::
   A. In the exec line above, xxx.yyy is a tag identifying your compiler and various build options, and will vary across pltaform.  (Note that GNU compilers must be at least version 7, and MPI should be at least of standard version 3).
   B. The example file can be any file from the `Tests/test_files` directories.
   C. In addition to informative output to the terminal, periodic plotfiles are written in the run folder.  These may be viewed with AMReX's `Amrvis <https://amrex-codes.github.io/amrex/docs_html/Visualization.html>`_ or `VisIt <https://visit-dav.github.io/visit-website/>`_:

      1. In VisIt, direct the File->Open dialogue to select the file named "Header" that is inside each plotfile folder.
      2. With Amrvis, `$ amrvis3d plt00030`, for example.

Alternatively, one can use cmake to build the code ::

    $ cd Build
    $ ./cmake.sh
    $ ./marbles example.inp

Dependencies
~~~~~~~~~~~~

`MARBLES` is built on the `AMReX <https://github.com/AMReX-Codes/amrex>`_ library.


Documentation
~~~~~~~~~~~~~

The full documentation for `MARBLES` exists in the Docs directory; at present this is maintained inline using `Sphinx <http://www.sphinx-doc.org>`_. To build the documentation ::

    $ cd Docs && mkdir build && cd build && sphinx-build -M html ../sphinx .

Or, using cmake ::

    $ cd Build && cmake -B build-docs ../Docs && cmake --build build-docs

Funding
~~~~~~~

This work was authored by the National Renewable Energy Laboratory (NREL) for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by the Laboratory Directed Research and Development (LDRD) Program at NREL. The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes. A portion of the research was performed using computational resources sponsored by the Department of Energy’s Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory.
