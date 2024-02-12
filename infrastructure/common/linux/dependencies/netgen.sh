#! /bin/bash
set -o errexit
set -o verbose

NETGEN="https://github.com/NGSolve/netgen.git"
BRANCH="v6.2.2307"

build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
	git submodule update --init --recursive
	cmake .  -LA -B cmbuild \
		-DCMAKE_BUILD_TYPE=Release \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DUSE_CCACHE:BOOL=OFF \
		-DUSE_CGNS:BOOL=OFF \
		-DUSE_CSG:BOOL=ON \
		-DUSE_GEOM2D:BOOL=ON \
		-DUSE_GUI:BOOL=OFF \
		-DUSE_INTERFACE:BOOL=ON \
		-DUSE_INTERNAL_TCL:BOOL=OFF \
		-DUSE_JPEG:BOOL=OFF \
		-DUSE_MPEG:BOOL=OFF \
		-DUSE_MPI:BOOL=OFF \
		-DUSE_MPI4PY:BOOL=OFF \
		-DUSE_NATIVE_ARCH:BOOL=OFF \
		-DUSE_NUMA:BOOL=OFF \
		-DUSE_OCC:BOOL=ON \
		-DUSE_PYTHON:BOOL=OFF \
		-DUSE_STLGEOM:BOOL=ON \
		-DUSE_SUPERBUILD:BOOL=OFF \
		-DENABLE_CPP_CORE_GUIDELINES_CHECK:BOOL=OFF \
		-DENABLE_UNIT_TESTS:BOOL=OFF
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
