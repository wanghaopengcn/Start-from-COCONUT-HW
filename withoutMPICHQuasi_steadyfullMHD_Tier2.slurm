#!/bin/bash -l
#SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_vswmc"
#SBATCH --cluster=wice
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
#SBATCH --nodes="20"
#SBATCH --ntasks-per-node="96"
#SBATCH --ntasks="1920"
#SBATCH --time="4:00:00"

PMIX_MCA_gds=hash
module load CMake/3.26.3-GCCcore-12.3.0
module load Boost/1.82.0-GCC-12.3.0
module load ParMETIS/4.0.3-gompi-2023a
module load PETSc/3.20.3-foss-2023a
##module load PETSc/3.15.1-foss-2021a

mpirun ./coolfluid-solver --scase /scratch/leuven/357/vsc35794/Coronal_modelling/TimeDependent/COCONUT_Map/CME_Relax_fullMHD.CFcase -malloc_debug



