#!/bin/bash
# Run the ZnO RDMFT example.
#
# Assumes abacus has been built with RDMFT support (ENABLE_RDMFT=ON)
# and is either on PATH or invoked via an absolute path.

# Uncomment to source an Intel oneAPI environment:
source /opt/intel/oneapi/setvars.sh

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# OpenMPI in containers/CI usually refuses to run as root without these
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_MCA_btl_vader_single_copy_mechanism=none

# Pick up the binary from PATH, or set ABACUS_PATH to your build tree:
#   e.g. ABACUS_PATH=$HOME/abacus-develop/build
#
export PATH=/home/kluo/Documents/repo/abacus-develop/build:$PATH
ABACUS=${ABACUS:-abacus_3p}

mpirun -np ${NP:-8} ${ABACUS} | tee log
