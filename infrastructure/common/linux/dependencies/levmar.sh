#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

LEVMAR_SOURCE=http://users.ics.forth.gr/~lourakis/levmar/levmar-2.6.tgz
LEVMAR_ARCHIVE=$(basename $LEVMAR_SOURCE)
LEVMAR_PATH="${LEVMAR_ARCHIVE%.*}"

build_and_install() {
	local source=$1
	local patchfile=$2
	pushd "$source" || exit 1
	dos2unix CMakeLists.txt CMakeLists.txt
	patch --ignore-whitespace -p0 < "$patchfile"
	cmake  . -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DCMAKE_POSITION_INDEPENDENT_CODE=ON \
		-DBUILD_DEMO:BOOLEAN=false
	pushd cmbuild || exit 1
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	# shellcheck disable=2155
	local patchfile="$(pwd)/levmar-install.patch"
	pushd "$BUILD_PATH" || exit 1
	download_source "$LEVMAR_SOURCE"
	tar xzf "$LEVMAR_ARCHIVE"
	build_and_install "$LEVMAR_PATH" "$patchfile"
	popd || exit 1
}

main
