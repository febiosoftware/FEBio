#! /bin/bash
set -o errexit
set -o verbose

export BUILD_PATH=/tmp/src
mkdir -p $BUILD_PATH

SETVARS="${SETVARS:-/opt/intel/oneapi/setvars.sh}"
. $SETVARS

main() {
	local dir=$1
	pushd $dir

	local installers=(hypre levmar mmg tetgen itk sitk occt netgen libzip ffmpeg)
	for installer in ${installers[@]}; do
		./${installer}.sh
	done

	popd
}

main $DEPENDENCIES_PATH
