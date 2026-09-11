#!/bin/bash
# Script to compile with necessary modules loaded

# Load required modules
module load PrgEnv-gnu/8.5.0
module load cuda/12.3
module load craype-x86-milan

# Show loaded modules
echo "Loaded modules:"
module list

# AMReX: use the pinned submodule (Submodules/AMReX).
# The GNUmakefile sets: AMREX_HOME ?= $(MARBLES_HOME)/Submodules/AMReX
# Unsetting any env AMREX_HOME here ensures local builds use the submodule,
# matching CI behavior exactly (GitHub Actions never has AMREX_HOME set).
unset AMREX_HOME

# Compile
echo ""
echo "Starting compilation with submodule AMReX ($(git -C ../Submodules/AMReX describe --tags 2>/dev/null)) ..."
make -j8 USE_CUDA=TRUE

