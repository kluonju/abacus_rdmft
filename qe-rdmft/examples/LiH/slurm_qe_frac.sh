#!/bin/bash
#SBATCH -J pw_lih_nspin2_q-2
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -t 10-24:00:00

#SBATCH --partition amd_a8_384
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

ulimit -s unlimited
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

module load oneAPI/2022.1 mpi/oneAPI/2022.1 intel/2021.2

export PATH="${HOME}/kluo/app/qe-rdmft/bin:${PATH}"
export LD_LIBRARY_PATH="${HOME}/kluo/app/qe-rdmft/lib:${LD_LIBRARY_PATH:-}"

mkdir -p tmp

NP="${SLURM_NTASKS:-4}"
mpirun -np "${NP}" pw.x -in lih.scf.in 2>&1 | tee lih.scf.log
