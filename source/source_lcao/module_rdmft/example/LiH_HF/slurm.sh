#!/bin/bash 
#SBATCH -J abacus
#SBATCH -o %j.out
#SBATCH -e %j.err 
#SBATCH -t 10-24:00:00

#SBATCH --partition amd_a8_384
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8

## SBATCH -x m4ci1604

ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load oneAPI/2022.1  mpi/oneAPI/2022.1 intel/2021.2 

# export PATH=/home/kluo/work/apps/repos/abacus-develop/build:$PATH
export PATH=~/kluo/app/kluo-abacus/build:$PATH
ABACUS="${ABACUS:-abacus_std_para}"

# /opt/intel/oneapi/mpi/2021.14/bin/mpirun -np 4 "${ABACUS}" 2>&1 | tee log
# srun --mpi=pmi2 "${ABACUS}"  2>&1  | tee log
# srun "${ABACUS}"  2>&1  | tee log
mpirun  "${ABACUS}"  2>&1  | tee log
