#! /bin/bash
# Uncomment next line if not global on target machine
# git config --global --add safe.directory /FEBio
. $(dirname $0)/cmake.sh

source "/opt/intel/oneapi/setvars.sh" --force
. ./cmake.sh
pushd cmbuild
make -j $(sysctl -n hw.ncpu)
popd
