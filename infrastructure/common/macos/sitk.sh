#!/bin/zsh
pushd $SOURCE_PATH

git clone --depth 1 --branch "v2.1.1.2" "https://github.com/SimpleITK/SimpleITK.git"
pushd SimpleITK
cmake .  -L -B cmbuild \
	-DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" \
  -DWRAP_DEFAULT:BOOL=OFF \
  -DBUILD_EXAMPLES:BOOL=OFF \
  -DBUILD_TESTING:BOOL=OFF \
  -DBUILD_SHARED_LIBS:BOOL=OFF \
  -DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_OSX_ARCHITECTURES="x86_64" \
	-DCMAKE_OSX_DEPLOYMENT_TARGET=10.15

pushd cmbuild
make install -j $(sysctl -n hw.ncpu)
popd
popd
rm -rf SimpleITK
popd
