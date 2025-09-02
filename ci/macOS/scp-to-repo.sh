#! /bin/bash

REMOTE_PATH="update2/FEBioStudio2Dev/macOS/stage"
if [ $# == 1 ] && [ "$1" != "develop" ]; then
    REMOTE_PATH="update2/FEBioStudio2Dev/branches/$1/macOS/stage"
fi

scp cmbuild/bin/* repo:~/$REMOTE_PATH/FEBioStudio.app/Contents/MacOS
scp cmbuild/lib/* repo:~/$REMOTE_PATH/FEBioStudio.app/Contents/Frameworks
ssh repo "chmod +x $REMOTE_PATH/FEBioStudio.app/Contents/MacOS/febio4"

if [ -f testLogs/Logs/* ]; then
    scp testLogs/Logs/* repo:~/TestSuite/Logs/macOS.txt
fi