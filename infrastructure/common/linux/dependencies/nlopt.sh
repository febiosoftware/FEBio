#! /bin/bash
set -o errexit
set -o verbose

NLOPT="https://github.com/stevengj/nlopt.git"
BRANCH="v2.10.1"

build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
	cmake .  -LA -B cmbuild \
		-DBUILD_SHARED_LIBS=OFF \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DNLOPT_GUILE=OFF \
		-DNLOPT_JAVA=OFF \
		-DNLOPT_MATLAB=OFF \
		-DNLOPT_OCTAVE=OFF \
		-DNLOPT_PYTHON=OFF \
		-DNLOPT_PYTHON_SABI=OFF \
		-DNLOPT_SWIG=OFF \
		-DNLOPT_TESTS=OFF 

	pushd cmbuild
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$NLOPT" "$BRANCH"
	popd || exit 1
}

main
