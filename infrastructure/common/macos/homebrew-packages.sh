#!/bin/zsh
set -ex

packages=(jq awscli eigen glew libssh libssh2)

arch -x86_64 $HOMEBREW_BIN  install "${packages[@]}"
