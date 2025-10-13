#!/bin/bash -l
##SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_sidc"
#SBATCH --cluster=genius
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
#SBATCH --nodes="10"
#SBATCH --ntasks-per-node="36"
#SBATCH --ntasks="360"
#SBATCH --time="2:00:00"

PMIX_MCA_gds=hash
module load CMake/3.26.3-GCCcore-12.3.0
module load Boost/1.82.0-GCC-12.3.0
module load ParMETIS/4.0.3-gompi-2023a
module load PETSc/3.20.3-foss-2023a
##module load PETSc/3.15.1-foss-2021a

mpirun ./coolfluid-solver --scase /scratch/leuven/357/vsc35794/Coronal_modelling/TimeDependent/COCONUT_Map/map_steady_fullMHD.CFcase -malloc_debug



