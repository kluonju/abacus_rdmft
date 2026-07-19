#!/usr/bin/env bash
#SBATCH --job-name=NiO-Power
#SBATCH --ntasks=24
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

# mkdir -p tmp
# srun --mpi=pmi2 pw.x -in rdmft.in 2>&1 | tee rdmft.log
mpirun -np 24 rdmft_dos.x -nk 2 -in dos.in 2>&1 | tee dos.log
