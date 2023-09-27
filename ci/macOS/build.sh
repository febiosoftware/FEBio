#! /bin/bash
# Uncomment next line if not global on target machine
git config --global --add safe.directory /FEBio
source "/opt/intel/oneapi/setvars.sh" --force
cmake . -B cmbuild -LA -DCMAKE_OSX_ARCHITECTURES="x86_64"
pushd cmbuild
make -j $(sysctl -n hw.ncpu)
popd
