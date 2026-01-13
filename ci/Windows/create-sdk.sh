#! /bin/bash
set -e
TARGET_DIR="${TARGET_DIR:-febio4-sdk}"
FEBIO_REPO="${FEBIO_REPO:-.}"
mkdir -p ${TARGET_DIR}/{include,lib,bin}
mkdir ${TARGET_DIR}/lib/{Release,Debug}
mkdir ${TARGET_DIR}/bin/Debug

# Copy in FEBioConfig.cmake
mkdir -p ${TARGET_DIR}/lib/cmake/FEBio
cp $FEBIO_REPO/FEBioConfig.cmake ${TARGET_DIR}/lib/cmake/FEBio

sdkDirs=(
    FECore
    FEBioMech
    FEBioMix
    FEBioFluid
    FEBioRVE
    FEBioPlot
    FEBioXML
    FEBioLib
    FEAMR
    FEBioOpt
    FEImgLib
)

for item in ${sdkDirs[@]}; do
    mkdir ${TARGET_DIR}/include/$item
    cp $FEBIO_REPO/$item/*.h ${TARGET_DIR}/include/$item
    
    # don't exit if there aren't .hpp files
    if ls "$FEBIO_REPO/$item"/*.hpp >/dev/null 2>&1; then
        cp "$FEBIO_REPO/$item"/*.hpp "${TARGET_DIR}/include/$item"
    fi 

    cp $FEBIO_REPO/cmbuild/lib/Release/$item.lib ${TARGET_DIR}/lib/Release
    cp $FEBIO_REPO/cmbuild/lib/Debug/$item.lib ${TARGET_DIR}/lib/Debug
done

cp $FEBIO_REPO/cmbuild/bin/Debug/febio4.exe ${TARGET_DIR}/bin/Debug
cp $FEBIO_REPO/cmbuild/bin/Debug/*.dll ${TARGET_DIR}/bin/Debug

