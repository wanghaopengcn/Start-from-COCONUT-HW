#!/bin/bash -l
#SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_vswmc"
#SBATCH --cluster=wice
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
#SBATCH --nodes="30"
#SBATCH --ntasks-per-node="18"
#SBATCH --ntasks="540"
#SBATCH --time="72:00:00"
##SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=12000M
##SBATCH --exclude=l33c25n2


module purge
##module load cluster/mindwell/batch
module load CMake/3.29.3-GCCcore-13.3.0
module load PETSc/3.24.4-gfbf-2024a-MPICH


export OMP_NUM_THREADS=1


mpiexec ./coolfluid-solver --scase /lustre1/project/stg_00007/HaopengWang/COCONUT_TR/map_unsteadyeclipse_fullMHD.CFcase -malloc_debug


