#!/bin/bash
# Script to compile with necessary modules loaded

# Load required modules
module load PrgEnv-gnu/8.5.0
module load cuda/12.3
module load craype-x86-milan

# Show loaded modules
echo "Loaded modules:"
module list

# AMReX: uses AMREX_HOME from your environment (set in ~/.bashrc or ~/.bash_profile).
# To use the pinned submodule instead, run:
#   git submodule update --init Submodules/AMReX
# and temporarily set: export AMREX_HOME=$(realpath ../Submodules/AMReX)

# Compile
echo ""
echo "Starting compilation with AMREX_HOME=${AMREX_HOME} ..."
make -j8 USE_CUDA=TRUE

