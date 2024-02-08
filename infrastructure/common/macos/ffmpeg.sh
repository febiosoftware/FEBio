#!/bin/zsh
set -e

arch -x86_64 $HOMEBREW_BIN tap homebrew-ffmpeg/ffmpeg
arch -x86_64 $HOMEBREW_BIN install homebrew-ffmpeg/ffmpeg/ffmpeg --with-openssl
