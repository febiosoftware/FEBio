#! /bin/bash
# Uncomment next line if not global on target machine
git config --global --add safe.directory /FEBio
source "/opt/intel/oneapi/setvars.sh" --force
cmake . -B febio -LA
pushd febio
make -j $(nproc)
popd
