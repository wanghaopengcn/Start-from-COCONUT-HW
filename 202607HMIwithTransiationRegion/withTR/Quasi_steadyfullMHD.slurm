#!/bin/bash -l
##SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_vswmc"
#SBATCH --cluster=mindwell
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
#SBATCH --nodes="5"
#SBATCH --ntasks-per-node="144"
#SBATCH --ntasks="720"
#SBATCH --time="10:00:00"

module purge
module load cluster/mindwell/batch
module load CMake/3.29.3-GCCcore-13.3.0
module load PETSc/3.24.4-gfbf-2024a-MPICH


export OMP_NUM_THREADS=1


mpiexec ./coolfluid-solver --scase /lustre1/project/stg_00007/HaopengWang/COCONUT_TR/map_steady_fullMHD.CFcase -malloc_debug



