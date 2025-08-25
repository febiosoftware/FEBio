#! /bin/bash

REMOTE_PATH="update2/FEBioStudio2Dev/Windows/stage"
if [ $# == 1 ] && [ "$1" != "develop" ]; then
    REMOTE_PATH="update2/FEBioStudio2Dev/branches/$1/Windows/stage"
fi

scp cmbuild/bin/Release/* repo:~/$REMOTE_PATH/bin

if [ -f testLogs/Logs/* ]; then
    scp testLogs/Logs/* repo:~/TestSuite/Logs/windows.txt
fi