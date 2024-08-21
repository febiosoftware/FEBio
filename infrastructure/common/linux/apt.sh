#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive
SUDO=""
if command -v sudo &> /dev/null
then
	SUDO=$(which sudo)
fi

$SUDO apt-get update
$SUDO apt-get install linux-headers-generic -y
$SUDO apt-get install software-properties-common wget gpg sudo -y
$SUDO apt-get update --fix-missing


wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
$SUDO apt-get update

DEBIAN_FRONTEND=noninteractive $SUDO apt-get install lua5.3 -y

xargs $SUDO apt-get install </tmp/linux/packages.txt -y
