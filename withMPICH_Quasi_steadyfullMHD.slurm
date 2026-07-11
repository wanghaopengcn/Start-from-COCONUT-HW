#!/bin/bash -l
#SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_vswmc"
#SBATCH --cluster=wice
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
#SBATCH --nodes="5"
#SBATCH --ntasks-per-node="96"
#SBATCH --ntasks="480"
#SBATCH --time="10:00:00"

module purge
module load CMake/3.29.3-GCCcore-13.3.0
module load PETSc/3.24.4-gfbf-2024a-MPICH


export OMP_NUM_THREADS=1


mpiexec ./coolfluid-solver --scase /scratch/leuven/357/vsc35794/Coronal_modelling/TimeDependent/COCONUT_Master/map_steady_fullMHD.CFcase -malloc_debug



