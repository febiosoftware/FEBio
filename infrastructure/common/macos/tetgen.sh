#!/bin/zsh
pushd $SOURCE_PATH

git clone "https://github.com/ufz/tetgen.git"
pushd tetgen
cmake .  -L -B cmbuild \
	-DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" \
	-DCMAKE_OSX_ARCHITECTURES="x86_64" \
	-DCMAKE_OSX_DEPLOYMENT_TARGET=10.15

pushd cmbuild
make -j $(sysctl -n hw.ncpu)
make install

popd
popd
rm -rf tetgen
popd
