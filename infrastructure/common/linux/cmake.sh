#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive

wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'

sudo apt update
sudo apt install cmake cmake-curses-gui -y
cmake --version
