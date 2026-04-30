#!/usr/bin/env bash
#SBATCH --job-name=h2-d_1_5
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=10-01:00:00
#SBATCH --partition=cpu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

ulimit -s unlimited
export OMP_NUM_THREADS=1
module load compiler/2022.1.0 mpi/2021.14 mkl/2022.1.0

export PATH=/home/kluo/work/apps/repos/abacus-develop/build:$PATH
ABACUS="${ABACUS:-abacus_std_para}"

# /opt/intel/oneapi/mpi/2021.14/bin/mpirun -np 4 "${ABACUS}" 2>&1 | tee log

srun --mpi=pmi2 "${ABACUS}"  2>&1  | tee log
