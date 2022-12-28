#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive
SUDO=""
if command -v sudo &> /dev/null
then
	SUDO=$(which sudo)
fi
$SUDO apt-get update
$SUDO apt-get install -y \
	qt6-base-dev \
	libqt6widgets6 \
	libqt6core5compat6-dev \
	libquazip5-dev
