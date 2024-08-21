#!/bin/bash
set -ex

if [ -f "${HOMEBREW_BIN}" ]; then
	echo "Homebrew already installed at: ${HOMEBREW_PREFIX}"
	exit
fi

git clone https://github.com/Homebrew/brew.git $HOMEBREW_PREFIX

cat << EOF >> /Users/$SSH_USER/.zprofile
export HOMEBREW_PREFIX="${HOMEBREW_PREFIX}"
export PATH="\$PATH:${HOMEBREW_PREFIX}:${HOMEBREW_PREFIX}/bin:${HOMEBREW_PREFIX}/Cellar:${HOMEBREW_PREFIX}/opt"
EOF
