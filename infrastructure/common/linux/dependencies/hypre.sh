#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

HYPRE_SOURCE="https://github.com/hypre-space/hypre/archive/refs/tags/v2.23.0.zip"
HYPRE_ARCHIVE=$(basename $HYPRE_SOURCE)
HYPRE_PATH="hypre-2.23.0"

build_and_install() {
	local source=$1
	pushd "$source" || exit 1
	pushd src || exit 1
	cmake .  -LA -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DHYPRE_HAVE_MPI=Off -DHYPRE_WITH_MPI=Off \
		-DCMAKE_POSITION_INDEPENDENT_CODE=On
	pushd cmbuild
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	download_source "$HYPRE_SOURCE"
	extract_source "$HYPRE_ARCHIVE"
	build_and_install "$HYPRE_PATH"
	popd || exit 1
}

main
