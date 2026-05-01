#!/bin/bash
# Script to compile with necessary modules loaded

# Load required modules
module load PrgEnv-gnu/8.5.0
module load cuda/12.3
module load craype-x86-milan

# Show loaded modules
echo "Loaded modules:"
module list

# AMReX selection:
#   - Installed 25.11 (Nov 2025, newer): /projects/hpcapps/nsawant/marblesLBM/amrex
#   - Pinned submodule (Aug 2025): ../Submodules/AMReX  (run: git submodule update --init)
# Unset AMREX_HOME to use the Makefile default (Submodules/AMReX if populated)
# or set it explicitly to override.
export AMREX_HOME=/projects/hpcapps/nsawant/marblesLBM/amrex

# Compile
echo ""
echo "Starting compilation with AMREX_HOME=${AMREX_HOME} ..."
make -j8 USE_CUDA=TRUE

