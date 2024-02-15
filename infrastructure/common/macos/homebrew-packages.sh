#!/bin/zsh
set -ex

packages=(jq awscli eigen glew libssh libssh2 yasm zstd pcre2 harfbuzz freetype pkg-config jpeg-turbo)

arch -x86_64 $HOMEBREW_BIN  install "${packages[@]}"
