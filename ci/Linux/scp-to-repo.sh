#! /bin/bash

REMOTE_PATH="update2/FEBioStudio2Dev/Linux/stage"
if [ $# == 1 ] && [ "$1" != "develop" ]; then
    REMOTE_PATH="update2/FEBioStudio2Dev/branches/$1/Linux/stage"
fi


set -e
scp cmbuild/bin/* repo:~/$REMOTE_PATH/bin
scp cmbuild/lib/* repo:~/$REMOTE_PATH/lib
ssh repo "chmod +x $REMOTE_PATH/bin/febio4"

# package and upload sdk
zip -r sdk.zip sdk/include
zip -r sdk.zip sdk/lib
scp sdk.zip repo:~/$REMOTE_PATH/

if [ -f testLogs/Logs/* ]; then
    scp testLogs/Logs/* repo:~/TestSuite/Logs/linux.txt
fi