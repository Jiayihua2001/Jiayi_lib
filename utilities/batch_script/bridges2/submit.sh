#!/bin/sh
#SBATCH -J test 
#SBATCH -n 128
#SBATCH -p RM 
#SBATCH -N 1 
#SBATCH -A mat210008p
#SBATCH -t 12:00:00
 
ulimit -s unlimited
ulimit -v unlimited
# Directory for aims binary and the env 
AIMS_DIR="/ocean/projects/mat210008p/shared/bin/"
AIMS_BIN="$AIMS_DIR/aims.221103_1.scalapack.mpi.x"
AIMS_ENV="$AIMS_DIR/aims_env.sh"

source $AIMS_ENV
export OMP_NUM_THREADS=1

mpirun -np 128 $AIMS_BIN > aims.out
