#! /bin/bash
# Uncomment next line if not global on target machine
set -e
#git config --global --add safe.directory /FEBio
source "/opt/intel/oneapi/setvars.sh" --force
cmake . -B cmbuild -LA \
	-DSET_DEVCOMMIT=ON \
	-DUSE_FFTW=ON \
	-DUSE_HYPRE=ON \
	-DUSE_LEVMAR=ON \
	-DUSE_MKL=ON \
	-DUSE_MMG=ON \
	-DUSE_STATIC_STDLIBS=ON \
	-DUSE_ZLIB=ON
pushd cmbuild
make -j $(nproc)
popd
