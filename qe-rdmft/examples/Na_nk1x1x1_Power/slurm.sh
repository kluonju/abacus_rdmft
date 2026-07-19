#!/usr/bin/env bash
#SBATCH --job-name=Na-Power
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=10-01:00:00
#SBATCH --partition=cpu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

ulimit -s unlimited
export OMP_NUM_THREADS=1

module load tbb/2023.0 tcm/1.5 umf/1.1.0 compiler-rt/2026.0.0 compiler/2026.0.0 mkl/2026.0 mpi/2021.18

export PATH="/home/kluo/work/apps/qe-rdmft/bin:${PATH}"
export LD_LIBRARY_PATH="/home/kluo/work/apps/qe-rdmft/lib:${LD_LIBRARY_PATH:-}"

mkdir -p tmp
mpirun -np 16 pw.x -nk 1 -in rdmft.in 2>&1 | tee rdmft.log
