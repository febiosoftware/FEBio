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
	libavformat-dev \
	libavcodec-dev \
	libavfilter-dev \
	libavutil-dev \
	libavdevice-dev \
	libswresample-dev