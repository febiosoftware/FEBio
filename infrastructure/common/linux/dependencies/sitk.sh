#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

SITK="https://github.com/SimpleITK/SimpleITK.git"
BRANCH="v2.1.1.2"
build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
	cmake .  -LA -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DWRAP_DEFAULT:BOOL=OFF \
		-DBUILD_EXAMPLES:BOOL=OFF \
		-DBUILD_TESTING:BOOL=OFF \
		-DBUILD_SHARED_LIBS:BOOL=OFF \
		-DCMAKE_BUILD_TYPE=Release

	pushd cmbuild
	make -j "$(nproc --ignore 2)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$SITK" "$BRANCH"
	popd || exit 1
}

main
