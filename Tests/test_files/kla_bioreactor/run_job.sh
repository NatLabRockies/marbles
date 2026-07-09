#!/bin/bash
#SBATCH -A hdcomb
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=160G
#SBATCH --gres=gpu:2
##SBATCH -p debug
#SBATCH -J wubbaLubbaDubDub
#SBATCH -o marbles_%j.out
#SBATCH -e marbles_%j.err

# Load modules (same as compile.sh)
module load PrgEnv-gnu/8.5.0
module load cuda/12.3
module load craype-x86-milan

# Enable Cray MPICH GPU-aware paths (requires libmpi_gtl_cuda linked into
# the executable; see Build/GNUmakefile).  Without this, ghost data would
# stage through host pinned buffers and 2-GPU runs would be ~3x SLOWER
# than 1-GPU due to D2H/MPI/H2D staging in FillBoundary.
export MPICH_GPU_SUPPORT_ENABLED=1

cd /scratch/nsawant/movingBody/marblesThermal/Tests/test_files/kla_bioreactor

# 2-GPU X-cut launch (default in the inp file: max_grid_size_x=90,
# _y=_z=180 -> 2 grids of 90x180x180, one per rank).  Do NOT pass
# --gpus-per-task=1: it forces per-task cgroup GPU isolation and breaks
# CUDA IPC for the GTL library ("cuIpcOpenMemHandle: invalid argument").
# Instead, let both ranks see both GPUs; AMReX round-robins device
# assignment (harmless "Multiple GPUs are visible" warning at startup).
srun --ntasks=2 --cpu-bind=cores \
    /scratch/nsawant/movingBody/marblesThermal/Build/marbles3d.gnu.TPROF.MPI.CUDA.ex \
    kla_bioreactor.inp
