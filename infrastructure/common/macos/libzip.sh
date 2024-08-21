#!/bin/zsh
pushd $SOURCE_PATH

repo="https://github.com/nih-at/libzip.git"
branch="v1.10.1"

git clone --depth 1 --branch "$branch" "$repo" "$branch"
pushd $branch
cmake .  -LA -B cmbuild \
  -DCMAKE_INSTALL_PREFIX="$INSTALLATION_PATH" \
  -DBUILD_DOC=OFF \
  -DBUILD_EXAMPLES=OFF \
  -DBUILD_OSSFUZZ=OFF \
  -DBUILD_REGRESS=OFF \
  -DBUILD_TOOLS=OFF \
  -DBUILD_SHARED_LIBS=ON \
  -DENABLE_BZIP2=OFF \
  -DENABLE_COMMONCRYPTO=OFF \
  -DENABLE_FDOPEN=OFF \
  -DENABLE_GNUTILS=OFF \
  -DENABLE_LZMA=OFF \
  -DENABLE_MBEDTLS=OFF \
  -DENABLE_OPENSSL=OFF \
  -DENABLE_WINDOWS_CRYPTO=OFF \
  -DLIBZIP_DO_INSTALL=ON \
  -DSHARED_LIB_VERSIONNING=ON \
  -DCMAKE_OSX_ARCHITECTURES="x86_64" \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.15

pushd cmbuild
make install -j $(sysctl -n hw.ncpu)
popd
popd
rm -rf $branch
popd
