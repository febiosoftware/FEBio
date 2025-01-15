#! /bin/bash

REMOTE_PATH="update2/FEBioStudio2Dev/Linux/stage"
if [ $# == 1 ] && [ "$1" != "develop" ]; then
    REMOTE_PATH="update2/FEBioStudio2Dev/branches/$1/Linux/stage"
fi


set -e
scp cmbuild/bin/* repo:~/$REMOTE_PATH/bin
scp cmbuild/lib/* repo:~/$REMOTE_PATH/lib
ssh repo "chmod +x $REMOTE_PATH/bin/febio4"


if [ -f testLogs/Logs/* ]; then
    scp testLogs/Logs/* repo:~/TestSuite/Logs/linux.txt
fi

if [ -f ChemArtifacts/lib/* ]; then
    scp ChemArtifacts/lib/* repo:~/$REMOTE_PATH/lib
fi

if [ -f HeatArtifacts/lib/* ]; then
    scp HeatArtifacts/lib/* repo:~/$REMOTE_PATH/lib
fi