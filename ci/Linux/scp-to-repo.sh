#! /bin/bash
set -e
scp cmbuild/bin/* repo:~/update2/FEBioStudio2Dev/Linux/stage/bin
scp cmbuild/lib/* repo:~/update2/FEBioStudio2Dev/Linux/stage/lib
ssh repo "chmod +x update2/FEBioStudio2Dev/Linux/stage/bin/febio4"


if [ -f testLogs/Logs/* ]; then
    scp testLogs/Logs/* repo:~/TestSuite/Logs/linux.txt
fi

if [ -f ChemArtifacts/lib/* ]; then
    scp ChemArtifacts/lib/* repo:~/update2/FEBioStudio2Dev/Linux/stage/lib
fi

if [ -f HeatArtifacts/lib/* ]; then
    scp HeatArtifacts/lib/* repo:~/update2/FEBioStudio2Dev/Linux/stage/lib
fi