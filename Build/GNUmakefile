MARBLES_HOME = ../
AMREX_HOME ?= $(MARBLES_HOME)/Submodules/AMReX

EBASE = marbles

DEBUG	= FALSE

DIM	= 3

COMP    = gnu

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_EB    = TRUE

TINY_PROFILE = TRUE
USE_PARTICLES = TRUE

BL_NO_FORT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(MARBLES_HOME)/Source/Make.package

INCLUDE_LOCATIONS += $(MARBLES_HOME)/Source
VPATH_LOCATIONS   += $(MARBLES_HOME)/Source

# ------------------------------------------------------------------
# Cray-MPICH GPU-aware MPI (GTL) — required for multi-GPU scaling on
# Kestrel-GPU.  Without this the AMReX FillBoundary path stages ghost
# data through host pinned buffers (D2H → MPI → H2D) which dominates
# runtime at ~75% on 2 ranks.  Linking libmpi_gtl_cuda lets AMReX
# auto-detect GPU-aware MPI (see AMReX_ParallelDescriptor.cpp:1551
# MPIX_GPU_query_support(MPIX_GPU_SUPPORT_CUDA)) and skip staging.
#
# Note: env var is named PE_MPICH_GTL_DIR_amd_gfx90a even for NVIDIA
# nodes — the CUDA .so lives in the same gtl/lib directory as the
# HIP one and is picked up by -lmpi_gtl_cuda.
#
# Runtime also needs: export MPICH_GPU_SUPPORT_ENABLED=1
# ------------------------------------------------------------------
ifeq ($(USE_CUDA),TRUE)
  ifeq ($(USE_MPI),TRUE)
    ifdef PE_MPICH_GTL_DIR_amd_gfx90a
      LIBRARIES += $(PE_MPICH_GTL_DIR_amd_gfx90a) -lmpi_gtl_cuda
    endif
  endif
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


