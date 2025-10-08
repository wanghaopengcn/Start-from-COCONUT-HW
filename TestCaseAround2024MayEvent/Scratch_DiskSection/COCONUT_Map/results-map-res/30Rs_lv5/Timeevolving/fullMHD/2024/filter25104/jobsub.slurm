#!/bin/bash -l
##SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_vswmc"
#SBATCH --cluster=genius
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
#SBATCH --nodes="1"
#SBATCH --ntasks-per-node="1"
#SBATCH --ntasks="1"
#SBATCH --time="48:00:00"

module load CMake/3.20.1-GCCcore-10.3.0
module load Boost/1.76.0-GCC-10.3.0
module load ParMETIS/4.0.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a

export OMP_NUM_THREADS=1
mpirun ./sperical.YY
