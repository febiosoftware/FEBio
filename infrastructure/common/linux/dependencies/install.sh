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

	local installers=(quazip hypre levmar mmg tetgen itk sitk occt netgen )
	for installer in ${installers[@]}; do
		./${installer}.sh
	done

	popd
}

main $DEPENDENCIES_PATH
