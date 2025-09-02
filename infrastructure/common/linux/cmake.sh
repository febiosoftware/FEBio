#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive

sudo apt update
sudo apt install libssl3 cmake cmake-curses-gui -y
cmake --version
