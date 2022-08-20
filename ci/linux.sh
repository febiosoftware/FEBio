#! /bin/bash
# Uncomment next line if not global on target machine
source "/opt/intel/oneapi/setvars.sh" --force
cmake . -B cmbuild -LA
pushd cmbuild
make -j $(nproc)
popd
