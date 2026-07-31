#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive
SUDO=""
if command -v sudo &> /dev/null
then
	SUDO=$(which sudo)
fi

INSTALLER="intel-oneapi-base-toolkit-2025.0.1.46_offline.sh"
aws s3 cp "s3://febiosoftware/linux/oneapi/${INSTALLER}" .
chmod +x ./$INSTALLER
$SUDO ./intel-oneapi-base-toolkit-2025.0.1.46_offline.sh -a --cli --silent --eula accept --action install --components intel.oneapi.lin.mkl.devel
rm ./$INSTALLER