#!/bin/zsh
pushd $SOURCE_PATH

git clone --depth 1 --branch "v5.2.1" "https://github.com/InsightSoftwareConsortium/ITK.git"
pushd ITK
cmake .  -L -B cmbuild \
	-DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" \
  -DBUILD_EXAMPLES:BOOL=OFF \
  -DBUILD_SHARED_LIBS:BOOL=OFF \
  -DBUILD_TESTING:BOOL=OFF \
  -DITK_USE_SYSTEM_EIGEN:BOOL=ON \
  -DITK_WRAP_PYTHON:BOOL=OFF \
  -DITK_DOXYGEN_HTML:BOOL=OFF \
  -DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_OSX_ARCHITECTURES="x86_64" \
	-DCMAKE_OSX_DEPLOYMENT_TARGET=10.15

pushd cmbuild
make install -j $(sysctl -n hw.ncpu)
popd
popd
rm -rf ITK
popd
