#!/bin/bash

#SBATCH -N 1
#SBATCH --time=00:10:00
#SBATCH -q RM
#SBATCH -A mat210008p
#SBATCH -J dimer_generation

# lets the stack grow unlimited to prevent program from crashing
ulimit -s unlimited
ulimit -v unlimited

# Limit number of threads in numpy
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# activate conda environment
module load anaconda3
module load openmpi/4.0.2-intel20.4
conda activate gnrs_env

chmod +x ui.conf
mpirun -np 128 gnrs ui.conf

