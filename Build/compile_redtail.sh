#!/bin/bash
# ============================================================================
# Marbles-thermal build script for the Utah CHPC Redtail AI supercomputer.
#
# System: HPE Cray XD670 nodes, 8x NVIDIA H200 SXM5 (141 GB HBM3e each),
#         NVLink/NVSwitch 900 GB/s intra-node, InfiniBand NDR 400 Gb/s
#         inter-node (rail-optimized fat-tree), Cray MPICH.
#
# Same GPU-aware MPI story as Kestrel: Cray MPICH provides libmpi_gtl_cuda,
# and our Build/GNUmakefile already links it when PE_MPICH_GTL_DIR_amd_gfx90a
# is set (see the GTL block in GNUmakefile).  H200 is compute capability 9.0,
# same as H100 -- no CUDA_ARCH change needed (Make.nrel already picks sm_90).
#
# TODO before first use on Redtail:
#   1. Verify the module names below against `module avail | grep -i cuda`
#      and `module avail | grep -i PrgEnv`.  Redtail may use different
#      versions than Kestrel (currently PrgEnv-gnu/8.5.0 + cuda/12.3).
#   2. Check `which_computer` detection in Submodules/AMReX/Tools/GNUMake/
#      Make.nrel -- it may reject Redtail's hostname.  If so, either:
#        a) add a case for the redtail hostname pattern to Make.nrel, or
#        b) create Submodules/AMReX/Tools/GNUMake/sites/Make.chpc and set
#           `WHICH_COMPUTER=redtail` in the environment, or
#        c) copy Make.nrel to Make.chpc and adjust hostname detection.
#   3. Check `echo $PE_MPICH_GTL_DIR_amd_gfx90a` after loading modules --
#      that env var must be set for the GNUmakefile GTL hook to fire.
#      If empty, the CPU-only fallback runs; multi-GPU perf will suffer.
# ============================================================================

set -e

# ---- modules (adjust names/versions to what Redtail actually provides) ----
module purge
module load PrgEnv-gnu        # provides gcc + Cray MPICH wrappers
module load cuda              # NVIDIA CUDA toolkit
# module load craype-x86-milan  # Kestrel-specific; Redtail is Xeon 8570 (Sapphire Rapids)
                                # -> either drop this line or load craype-x86-emeraldrapids
                                #    (or whatever Redtail's craype cpu module is called).

echo "Loaded modules:"
module list

echo
echo "Sanity checks:"
echo "  PE_MPICH_GTL_DIR_amd_gfx90a = ${PE_MPICH_GTL_DIR_amd_gfx90a:-<UNSET -- GTL will NOT link>}"
echo "  MPICH_DIR                    = ${MPICH_DIR:-<UNSET>}"
echo "  CUDA path                    = ${CUDA_ROOT:-${CUDA_HOME:-${CUDA_PATH:-<UNSET>}}}"
echo

# AMReX submodule
unset AMREX_HOME

# Force AMReX to think we're on Kestrel so Make.nrel picks up (fast path).
# Remove this line after adding a proper Redtail case to Make.nrel or a new
# Make.chpc site file.
export WHICH_COMPUTER=kestrel

echo
echo "Starting compilation with submodule AMReX ($(git -C ../Submodules/AMReX describe --tags 2>/dev/null)) ..."
make -j16 USE_CUDA=TRUE
