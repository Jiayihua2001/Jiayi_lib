#!/bin/bash

#SBATCH -N 1
#SBATCH --time=1:00:00
#SBATCH -q RM
#SBATCH -A mat210008p
#SBATCH -J ase_aims_test

# lets the stack grow unlimited to prevent program from crashing
ulimit -s unlimited
ulimit -v unlimited

# Limit number of threads in numpy
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

# activate conda environment

module load intelmpi/2021.3.0-intel2021.3.0
source /opt/intel/oneapi/intelpython/latest/bin/activate
conda activate aims_run

chmod +x ase_aims.py
python ase_aims.py

