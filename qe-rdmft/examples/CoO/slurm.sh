#!/usr/bin/env bash
#SBATCH --job-name=CoO-Power
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --time=10-01:00:00
#SBATCH --partition=cpu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

ulimit -s unlimited
export OMP_NUM_THREADS=2

module load tbb/2023.0 tcm/1.5 umf/1.1.0 compiler-rt/2026.0.0 compiler/2026.0.0 mkl/2026.0 mpi/2021.18

export PATH="/home/kluo/work/apps/qe-rdmft/bin:${PATH}"
export LD_LIBRARY_PATH="/home/kluo/work/apps/qe-rdmft/lib:${LD_LIBRARY_PATH:-}"

mkdir -p tmp
# K-point pool parallelism (-nk / -npool > 1) is supported: with
# nspin = 2 the pools split by spin (here 2 pools = spin up / spin
# down).  Note -nk is just a synonym for -npool, so pass only ONE of
# them (the original "-npool 1 -nk 2" was contradictory -- the trailing
# -nk 2 won and selected 2 pools).
# srun --mpi=pmi2 pw.x -nk 2 -in rdmft.in 2>&1 | tee rdmft.log
mpirun -np 2 pw.x -nk 2 -in rdmft.in 2>&1 | tee rdmft.log
