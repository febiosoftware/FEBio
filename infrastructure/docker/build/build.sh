#! /bin/bash
set -o errexit
set -o verbose

source /opt/intel/oneapi/setvars.sh

main() {
	mkdir -p cmbuild-linux
	pushd cmbuild-linux
	cmake .. -B .
	make -j "$(nproc)"
	popd || exit 1
}

main
