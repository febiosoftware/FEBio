#!/bin/bash
set -ex
export DEBIAN_FRONTEND=noninteractive
pushd /tmp
INSTALLER="installbuilder-enterprise-23.11.0-linux-x64-installer.run"
aws s3 cp "s3://febiosoftware/linux/installbuilder/${INSTALLER}" .
chmod +x ./$INSTALLER

sudo ./$INSTALLER \
	--mode unattended \
	--installer-language en

sudo ln -s /opt/installbuilder-23.11.0/bin/builder /usr/local/bin/
builder --version
