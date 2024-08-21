#!/bin/zsh
set -e

REPO="https://github.com/FFmpeg/FFmpeg.git"
BRANCH="n6.1"

pushd $SOURCE_PATH
rm -rf "${SOURCE_PATH}/${BRANCH}"

export MACOSX_DEPLOYMENT_TARGET=10.15
export MACOSX_DEPLOYMENT_ARCHITECTURES=x86_64

git clone --depth 1 --branch $BRANCH $REPO $BRANCH

pushd $BRANCH

arch -x86_64 ./configure \
  --disable-everything \
  --disable-programs \
  --disable-doc \
  --disable-static \
  --disable-debug \
  --enable-shared \
  --enable-encoder=mpeg1video \
  --enable-muxer=mpeg1video \
  --prefix="$INSTALLATION_PATH"

arch -x86_64 make -j $(sysctl -n hw.ncpu)
arch -x86_64 make install

popd
popd
rm -rf $BRANCH
