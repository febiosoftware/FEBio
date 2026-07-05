#! /bin/bash

REMOTE_PATH="/serverRoot/update2/FEBioStudio2Dev/Linux/stage"
if [ $# == 1 ] && [ "$1" != "develop" ]; then
    REMOTE_PATH="/serverRoot/update2/FEBioStudio2Dev/branches/$1/Linux/stage"
fi


set -e
scp cmbuild/bin/* febio-web:$REMOTE_PATH/bin
scp cmbuild/lib/* febio-web:$REMOTE_PATH/lib
ssh febio-web "chmod +x $REMOTE_PATH/bin/febio4"

# package and upload sdk
pushd sdk
zip -r sdk.zip include
zip -r sdk.zip lib
scp sdk.zip febio-web:$REMOTE_PATH/
popd