#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

MMG_SOURCE=https://github.com/MmgTools/mmg/archive/master.zip
MMG_ARCHIVE=$(basename $MMG_SOURCE)
MMG_PATH="mmg-${MMG_ARCHIVE%.*}"

build_and_install() {
	local source=$1
	pushd "$source" || exit 1
	cmake . -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DCMAKE_POSITION_INDEPENDENT_CODE=ON
	pushd cmbuild || exit 1
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	download_source "$MMG_SOURCE"
	extract_source "$MMG_ARCHIVE"
	build_and_install "$MMG_PATH"
	popd || exit 1
}

main
