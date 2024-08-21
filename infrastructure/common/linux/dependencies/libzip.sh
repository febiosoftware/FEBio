#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

SOURCE="https://github.com/nih-at/libzip.git"
BRANCH="v1.10.1"

build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
	cmake .  -LA -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
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
		-DSHARED_LIB_VERSIONNING=ON
	pushd cmbuild
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$SOURCE" "$BRANCH"
	popd || exit 1
}

main
