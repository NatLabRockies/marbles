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
# Cray-MPICH pkg-config exports PE_MPICH_GTL_DIR_<accelerator> for
# every accelerator it knows about (nvidia70/80/90, amd_gfx90a/942,
# ponteVecchio, ...), but every DIR points to the SAME directory —
# /opt/cray/pe/mpich/*/gtl/lib — which contains libmpi_gtl_cuda.so
# alongside libmpi_gtl_hsa.so and libmpi_gtl_ze.so.  So we probe any
# NVIDIA-suffixed export (with an AMD fallback so this block still
# fires on Cray systems that ship only the ROCm hooks) and hardcode
# -lmpi_gtl_cuda since we're always building for NVIDIA here.
#
# Runtime also needs: export MPICH_GPU_SUPPORT_ENABLED=1
# ------------------------------------------------------------------
ifeq ($(USE_CUDA),TRUE)
  ifeq ($(USE_MPI),TRUE)
    CRAY_GTL_DIR := $(firstword                \
        $(PE_MPICH_GTL_DIR_nvidia90)           \
        $(PE_MPICH_GTL_DIR_nvidia80)           \
        $(PE_MPICH_GTL_DIR_nvidia70)           \
        $(PE_MPICH_GTL_DIR_amd_gfx90a))
    ifneq ($(CRAY_GTL_DIR),)
      LIBRARIES += $(CRAY_GTL_DIR) -lmpi_gtl_cuda
    endif
  endif
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


