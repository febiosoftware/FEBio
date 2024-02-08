#!/bin/zsh
pushd $SOURCE_PATH

git clone https://github.com/hypre-space/hypre.git
pushd hypre/src
cmake .  -L -B cmbuild \
	-DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" \
	-DHYPRE_HAVE_MPI=Off \
	-DHYPRE_WITH_MPI=Off \
	-DCMAKE_OSX_ARCHITECTURES="x86_64" \
	-DCMAKE_OSX_DEPLOYMENT_TARGET=10.15

pushd cmbuild
make install -j $(sysctl -n hw.ncpu)
popd
popd
rm -rf hypre
popd
