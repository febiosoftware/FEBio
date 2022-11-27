#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

NETGEN="https://github.com/NGSolve/netgen.git"
BRANCH="v6.2.2204"
build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
	cmake .  -LA -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DUSE_PYTHON=OFF \
		-DUSE_GUI=OFF \
		-DUSE_NATIVE_ARCH=OFF \
		-DUSE_OCC=OFF
	pushd cmbuild
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$NETGEN" "$BRANCH"
	popd || exit 1
}

main
