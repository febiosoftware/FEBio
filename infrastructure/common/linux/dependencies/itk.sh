#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

ITK="https://github.com/InsightSoftwareConsortium/ITK.git"
BRANCH="v5.2.1"
build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
	cmake . -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DBUILD_EXAMPLES:BOOL=OFF \
		-DBUILD_SHARED_LIBS:BOOL=OFF \
		-DBUILD_TESTING:BOOL=OFF \
		-DITK_WRAP_PYTHON:BOOL=OFF \
		-DITK_DOXYGEN_HTML:BOOL=OFF \
		-DCMAKE_BUILD_TYPE=Release

	pushd cmbuild
	make -j "$(nproc --ignore 2)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$ITK" "$BRANCH"
	popd || exit 1
}

main
