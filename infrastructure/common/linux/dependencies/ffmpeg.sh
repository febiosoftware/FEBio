#!/bin/bash
set -e

REPO="https://github.com/FFmpeg/FFmpeg.git"
BRANCH="n6.1"

build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"

  pushd $BRANCH

  ./configure \
    --disable-everything \
    --disable-programs \
    --disable-doc \
    --disable-static \
    --disable-debug \
    --enable-shared \
    --enable-encoder=mpeg1video \
    --enable-muxer=mpeg1video \
    --prefix="/usr/local"

  sudo make -j $(nproc)
  sudo make install
  popd
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$REPO" "$BRANCH"
	popd || exit 1
}

main
