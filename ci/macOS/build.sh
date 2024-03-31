#! /bin/bash
# Uncomment next line if not global on target machine
# git config --global --add safe.directory /FEBio

source "/opt/intel/oneapi/setvars.sh" --force
. $(dirname $0)/cmake.sh
pushd cmbuild
make -j $(sysctl -n hw.ncpu)
popd
