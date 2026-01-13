#! /bin/bash
set -e
TARGET_DIR="${TARGET_DIR:-febio4-sdk}"
FEBIO_REPO="${FEBIO_REPO:-.}"
mkdir -p ${TARGET_DIR}/{include,lib}

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

sdkLibs=(
    libfecore.so
    libfebiomech.so
    libfebiomix.so
    libfebiofluid.so
    libfebiorve.so
    libfebioplot.a
    libfebioxml.a
    libfebiolib.so
    libfeamr.so
    libfebioopt.so
    libfeimglib.so
)

for item in ${sdkDirs[@]}; do
    mkdir ${TARGET_DIR}/include/$item
    cp $FEBIO_REPO/$item/*.h ${TARGET_DIR}/include/$item
    
    # don't exit if there aren't .hpp files
    if ls "$FEBIO_REPO/$item"/*.hpp >/dev/null 2>&1; then
        cp "$FEBIO_REPO/$item"/*.hpp "${TARGET_DIR}/include/$item"
    fi 
done

for item in ${sdkLibs[@]}; do
    cp $FEBIO_REPO/cmbuild/lib/$item ${TARGET_DIR}/lib
done

