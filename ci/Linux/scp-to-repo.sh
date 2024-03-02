#! /bin/bash
set -e
scp cmbuild/bin/* repo:~/update2/FEBioStudio2Dev/Linux/stage/bin
scp cmbuild/lib/* repo:~/update2/FEBioStudio2Dev/Linux/stage/lib
ssh repo "chmod +x update2/FEBioStudio2Dev/Linux/stage/bin/febio4"
