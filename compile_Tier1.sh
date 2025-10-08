module load CMake/3.26.3-GCCcore-12.3.0
module load PETSc/3.20.3-foss-2023a
export TOP_DIR="/readonly${VSC_SCRATCH_PROJECTS_BASE}/2025_006/HaopengWang/"
export COOLFLUID_TOP_DIR="${TOP_DIR}/COOLFluiD_HortenseMaster"

export BUILD_MODE=optim
export CONF_FILE="COOLFluid_Hortense_nocuda.conf"


export COOLFLUID_BASEBUILD_DIR="${COOLFLUID_TOP_DIR}/OPENMPI"
export COOLFLUID_CONF_FILE="${COOLFLUID_TOP_DIR}/${CONF_FILE}"
export COOLFLUID_INSTALL_DIR="${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}/INSTALL"
export ALL_ACTIVE=1

cd $COOLFLUID_TOP_DIR
#./prepare.pl --config-file=${COOLFLUID_CONF_FILE} --build=${BUILD_MODE}

cd ${COOLFLUID_BASEBUILD_DIR}/${BUILD_MODE}
make -j 4
#make install
