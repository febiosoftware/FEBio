#! /bin/bash
set -e

source "/opt/intel/oneapi/setvars.sh" --force
. $(dirname $0)/cmake.sh
pushd cmbuild
make -j $(sysctl -n hw.ncpu)
popd
