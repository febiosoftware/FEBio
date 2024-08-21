#! /bin/bash
set -o errexit
set -o verbose

REPO=https://github.com/jturney/levmar.git
DIR=levmar

build_and_install() {
	local repo=$1
	local dir=$2
	git clone $repo $dir
	pushd $dir
	cmake  . -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DCMAKE_POSITION_INDEPENDENT_CODE=ON \
		-DBUILD_DEMO:BOOLEAN=false
	pushd cmbuild
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd $BUILD_PATH
	build_and_install $REPO $DIR
	popd || exit 1
}

main
