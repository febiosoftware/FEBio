#!/bin/bash
set -e
export DEBIAN_FRONTEND=noninteractive

if command -v git &> /dev/null
then
	git --version
fi
sudo add-apt-repository ppa:git-core/ppa -y
sudo apt-get update
sudo apt-get install git -y
git --version
