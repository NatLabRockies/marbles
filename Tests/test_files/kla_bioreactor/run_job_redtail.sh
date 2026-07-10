#!/bin/bash
# ============================================================================
# Marbles-thermal kLa bioreactor run script for Utah CHPC Redtail.
#
# Redtail node: 8x NVIDIA H200 SXM5 (141 GB HBM3e each), NVSwitch 900 GB/s
# intra-node, InfiniBand NDR 400 Gb/s inter-node.  Same GPU-aware MPI story
# as Kestrel (Cray MPICH + libmpi_gtl_cuda).
#
# TODO before first submission:
#   1. Replace #SBATCH -A with your Redtail allocation account.
#   2. Verify the partition name (Kestrel uses `debug`, `gpu-h100`; Redtail
#      may use something like `redtail-gpu` -- check `sinfo`).
#   3. Verify --gres syntax: some Cray clusters need --gres=gpu:h200:8
#      instead of --gres=gpu:8.
#   4. Check the module names in the load block match compile_redtail.sh.
#   5. If Redtail's MPI wasn't linked with GTL, comment out MPICH_GPU_SUPPORT
#      and expect ~3x slower FillBoundary (see LBM.cpp:refill_and_spill
#      comment header).
# ============================================================================

#SBATCH -A hdcomb                       # TODO Redtail allocation account
#SBATCH -t 24:00:00
#SBATCH --nodes=1                       # bump to 2/4/8 for scaling runs
#SBATCH --ntasks-per-node=8             # 1 MPI rank per H200
#SBATCH --gres=gpu:8                    # TODO may need gpu:h200:8 on Redtail
#SBATCH --mem=1500G                     # 2 TB per node; leave headroom
#SBATCH --cpus-per-task=14              # 56 cores / 8 ranks = 7 per rank
                                        # x2 hyperthread pair = 14 (adjust if
                                        # hyperthreading is off on Redtail)
#SBATCH -J marbles_kla_8gpu
#SBATCH -o marbles_%j.out
#SBATCH -e marbles_%j.err

# ---- modules ----
module purge
module load PrgEnv-gnu
module load cuda
# module load craype-x86-emeraldrapids  # TODO check Redtail CPU module name

# ---- GPU-aware MPI ----
export MPICH_GPU_SUPPORT_ENABLED=1

# ---- inputs ----
cd /path/to/marblesThermal/Tests/test_files/kla_bioreactor    # TODO Redtail path

# ============================================================================
# STRONG-SCALING RUN on 8 GPUs.  Default inp uses X-cut mgs=90 which was
# tuned for 2 ranks; override for 8 ranks with an isotropic 2x2x2 cube
# split (mgs=90 with no per-axis overrides).
#
#   180^3 domain / mgs=90 -> 8 grids of 90^3, one per rank.
#   Each rank has 730k cells (well-tuned for H200 kernel occupancy).
#   Each rank has 3 face + 12 edge + 8 corner ghost exchanges -- all
#   over NVSwitch (~900 GB/s) since everything is on one node.
#
# For WEAK-SCALING demo (constant 5.8M cells/rank), scale n_cell with the
# rank count:
#   1 GPU: amr.n_cell = 180 180 180  amr.max_grid_size = 180
#   2 GPU: amr.n_cell = 226 226 226  amr.max_grid_size_x = 113 ...
#   4 GPU: amr.n_cell = 285 285 285  ...
#   8 GPU: amr.n_cell = 360 360 360  amr.max_grid_size = 180 (2x2x2 split)
#
# DO NOT pass --gpus-per-task=1: it forces per-task cgroup GPU isolation
# which breaks CUDA IPC for the GTL library (cuIpcOpenMemHandle error).
# ============================================================================
srun --ntasks=8 --cpu-bind=cores \
    /path/to/marblesThermal/Build/marbles3d.gnu.TPROF.MPI.CUDA.ex \
    kla_bioreactor.inp \
    amr.max_grid_size=90 \
    amr.max_grid_size_x=90 \
    amr.max_grid_size_y=90 \
    amr.max_grid_size_z=90
