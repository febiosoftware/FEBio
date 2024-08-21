#!/bin/zsh
set -ex

if [ -d "${INSTALLATION_PATH}" ]; then
  echo "Paths already configured"
  exit
fi

mkdir -p "${SOURCE_PATH}"
mkdir -p "${INSTALLATION_PATH}"
chown -R "${SSH_USER}" "${INSTALLATION_PATH}"

cat << EOF >> /Users/$SSH_USER/.zprofile
export PATH="\$PATH:${INSTALLATION_PATH}"
EOF
