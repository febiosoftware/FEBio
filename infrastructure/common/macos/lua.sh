#!/bin/zsh
set -ex

arch -x86_64 $HOMEBREW_BIN install lua@5.3

cat << EOF >> /Users/$SSH_USER/.zprofile
export PATH="${HOMEBREW_PREFIX}/opt/lua@5.3/bin:\$PATH"
EOF
