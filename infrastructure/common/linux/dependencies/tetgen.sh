#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

TETGEN_GIT="https://github.com/ufz/tetgen.git"
TETGEN_FILENAME=$(basename $TETGEN_GIT)
TETGEN_SOURCE_DIR="${TETGEN_FILENAME%.*}"

build_and_install() {
	local git_remote=$1
	local src_dir=$2

	git clone $git_remote $src_dir
	pushd "${src_dir}" || exit 1
	cmake .  -LA -B cmbuild -DCMAKE_INSTALL_PREFIX="/usr/local"
	pushd cmbuild
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$TETGEN_GIT" "$TETGEN_SOURCE_DIR"
	popd || exit 1
}

main
