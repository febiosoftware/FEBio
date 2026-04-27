#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

MMG="https://github.com/MmgTools/mmg.git"
BRANCH="v5.7.3"

build_and_install() {
    local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
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
	build_and_install "$MMG" "$BRANCH"
	popd || exit 1
}

main
