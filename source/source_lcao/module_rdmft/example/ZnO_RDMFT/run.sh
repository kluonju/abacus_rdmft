#!/bin/bash

#source /home/kluo/Documents/repos/kluo/abacus-develop/toolchain/abacus_env.sh
#source ~/.bashrc
#activate_oneapi
source /opt/intel/oneapi/setvars.sh 
ABACUS_PATH=/home/kluo/Documents/repo/abacus-develop/build
#ABACUS_PATH=/home/kluo/Documents/repos/kluo/abacus-develop/build-lts
export OMP_NUM_THREADS=1
#mpirun -np 4 $ABACUS_PATH/abacus 
#mpirun -np 2 $ABACUS_PATH/abacus  1> log 2>err
#mpirun -np 1 $ABACUS_PATH/abacus_3p  1> log 
mpirun -np 1 $ABACUS_PATH/abacus_3p  |tee  log 
#$ABACUS_PATH/abacus 
#mpirun -np 4 abacus 
