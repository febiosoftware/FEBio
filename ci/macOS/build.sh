#! /bin/bash
# Uncomment next line if not global on target machine
git config --global --add safe.directory /FEBio
source "/opt/intel/oneapi/setvars.sh" --force
cmake . -B cmbuild -L -DCMAKE_OSX_ARCHITECTURES="x86_64" -DCMAKE_OSX_DEPLOYMENT_TARGET=10.13
pushd cmbuild
make -j8 #$(sysctl -n hw.ncpu)
#popd
