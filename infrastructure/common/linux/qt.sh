#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive
SUDO=""
if command -v sudo &> /dev/null
then
	SUDO=$(which sudo)
fi

$SUDO pip install aqtinstall -v
$SUDO aqt install-qt --outputdir /opt/Qt linux desktop 6.9.3 -m all
