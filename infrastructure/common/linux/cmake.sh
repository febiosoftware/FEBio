#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive

wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null
sudo apt-add-repository -y 'deb https://apt.kitware.com/ubuntu/ jammy main'

sudo apt update
sudo apt install libssl3 cmake cmake-curses-gui -y
cmake --version
