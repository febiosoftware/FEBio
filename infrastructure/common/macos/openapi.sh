#!/bin/bash
set -e

if [ -f /opt/intel/oneapi/setvars.sh ]; then
	echo "Oneapi was previously installed"
	exit 0
fi

TARGET_DIR="/tmp/intel"
ONEAPI_FILE="m_BaseKit_p_2023.2.0.49398_offline"
ONEAPI_DMG="$ONEAPI_FILE.dmg"
ONEAPI_URI="https://registrationcenter-download.intel.com/akdlm/IRC_NAS/cd013e6c-49c4-488b-8b86-25df6693a9b7/$ONEAPI_DMG"
curl --create-dirs -O --output-dir "$TARGET_DIR" "$ONEAPI_URI"
hdiutil attach "$TARGET_DIR/$ONEAPI_DMG"
pushd /Volumes/$ONEAPI_FILE/bootstrapper.app/Contents/MacOS/
chmod +x ./install.sh
arch -x86_64 ./install.sh --silent --eula accept
popd
hdiutil detach "/Volumes/$ONEAPI_FILE"
