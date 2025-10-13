#!/bin/bash -l
##SBATCH --partition=batch_sapphirerapids
#SBATCH --account="lp_sidc"
#SBATCH --cluster=wice
#SBATCH --job-name="def"
#SBATCH --mail-type="END,FAIL"
##SBATCH --nodes="30"
##SBATCH --ntasks-per-node="10"
#SBATCH --nodes="18"
#SBATCH --ntasks-per-node="20"
#SBATCH --ntasks="360"
#SBATCH --time="48:59:00"
##SBATCH --partition=bigmem
#SBATCH --mem-per-cpu=10000M


PMIX_MCA_gds=hash
module load CMake/3.20.1-GCCcore-10.3.0
module load Boost/1.76.0-GCC-10.3.0
module load ParMETIS/4.0.3-gompi-2021a
module load PETSc/3.20.3-foss-2023a

mpirun ./coolfluid-solver --scase /scratch/leuven/357/vsc35794/Coronal_modelling/TimeDependent/COCONUT_Map/map_unsteadyeclipse_fullMHD.CFcase -malloc_debug


