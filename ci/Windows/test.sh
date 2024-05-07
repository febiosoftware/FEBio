#! /usr/bin/bash
# Uncomment next line if not global on target machine
set -e
set -x
FEBIO_XML="$(dirname $0)/febio.xml"
FEBIO_DIR=$(realpath ./cmbuild/bin/Release)
FEBIO_LIB=$(realpath ./cmbuild/bin/Release)
FEBIO_BIN="${FEBIO_DIR}/febio4.exe"

FEBIOHEAT=$(realpath ./febioheat/Release/FEBioHeat.dll)
FEBIOCHEM=$(realpath ./febiochem/Release/FEBioChem.dll)

ONEAPI=$(cygpath -u "$ONEAPI_ROOT") #Convert $ONEAPI_ROOT to posix path
PYTHON="${ONEAPI}intelpython/latest/python"

IOMP_LIB="${ONEAPI}compiler/latest/windows/redist/intel64_win/compiler/libiomp5md.dll"

cp -a "$IOMP_LIB" "$FEBIO_LIB"

TESTSUITE=./TestSuite

if [[ ! -d $TESTSUITE ]]; then
	echo "Error: TestSuite was not located" >&2; exit 1
fi

if [[ -f "$FEBIOHEAT" ]]; then
    cp $FEBIOHEAT $FEBIO_LIB
fi

if [[ -f "$FEBIOCHEM" ]]; then
    cp $FEBIOCHEM $FEBIO_LIB
fi

# Copy configuration in
cp $FEBIO_XML $FEBIO_DIR
"$PYTHON" ./TestSuite/code/tools.py -r $FEBIO_BIN
