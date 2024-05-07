#! /bin/bash
set -o errexit
set -o verbose

# shellcheck disable=1091
. ./common-functions.sh

OCCT="https://github.com/Open-Cascade-SAS/OCCT.git"
BRANCH="V7_7_2"
build_and_install() {
	local source=$1
	local branch=$2

	git clone --depth 1 --branch "$branch" "$source" "$branch"
	pushd $branch || exit 1
	cmake . -B cmbuild \
		-DCMAKE_INSTALL_PREFIX="/usr/local" \
		-DBUILD_DOC_Overview:BOOL=OFF \
		-DBUILD_ENABLE_FPE_SIGNAL_HANDLER:BOOL=OFF \
		-DBUILD_Inspector:BOOL=OFF \
		-DBUILD_LIBRARY_TYPE:STRING=Shared \
		-DBUILD_MODULE_ApplicationFramework:BOOL=OFF \
		-DBUILD_MODULE_DETools:BOOL=OFF \
		-DBUILD_MODULE_DataExchange:BOOL=ON \
		-DBUILD_MODULE_Draw:BOOL=OFF \
		-DBUILD_MODULE_FoundationClasses:BOOL=TRUE \
		-DBUILD_MODULE_ModelingAlgorithms:BOOL=TRUE \
		-DBUILD_MODULE_ModelingData:BOOL=OFF \
		-DBUILD_MODULE_Visualization:BOOL=OFF \
		-DBUILD_OPT_PROFILE:STRING=Default \
		-DBUILD_RELEASE_DISABLE_EXCEPTIONS:BOOL=OFF \
		-DBUILD_RESOURCES:BOOL=OFF \
		-DBUILD_SAMPLES_QT:BOOL=OFF \
		-DBUILD_USE_PCH:BOOL=OFF \
		-DUSE_DRACO:BOOL=OFF \
		-DUSE_FFMPEG:BOOL=OFF \
		-DUSE_FREEIMAGE:BOOL=OFF \
		-DUSE_FREETYPE:BOOL=OFF \
		-DUSE_GLES2:BOOL=OFF \
		-DUSE_OPENGL:BOOL=OFF \
		-DUSE_OPENVR:BOOL=OFF \
		-DUSE_RAPIDJSON:BOOL=OFF \
		-DUSE_TBB:BOOL=OFF \
		-DUSE_TK:BOOL=OFF \
		-DUSE_VTK:BOOL=OFF \
		-DUSE_XLIB:BOOL=OFF

	pushd cmbuild
	make -j "$(nproc)"
	sudo make install
	popd || exit 1
	popd || exit 1
}

main() {
	pushd "$BUILD_PATH" || exit 1
	build_and_install "$OCCT" "$BRANCH"
	popd || exit 1
}

main
