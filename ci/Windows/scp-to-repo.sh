#! /bin/bash

REMOTE_PATH="/serverRoot/update2/FEBioStudio2Dev/Windows/stage"
if [ $# == 1 ] && [ "$1" != "develop" ]; then
    REMOTE_PATH="/serverRoot/update2/FEBioStudio2Dev/branches/$1/Windows/stage"
fi

scp cmbuild/bin/Release/* repo:$REMOTE_PATH/bin

# package and upload sdk
pushd sdk
/c/WINDOWS/system32/tar -acf sdk.zip include lib bin
scp sdk.zip repo:$REMOTE_PATH/
popd

if [ -f testLogs/Logs/* ]; then
    scp testLogs/Logs/* repo:/serverRoot/TestSuite/Logs/windows.txt
fi