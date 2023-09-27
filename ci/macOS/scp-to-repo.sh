#! /bin/bash
scp cmbuild/bin/* repo:~/update2/FEBioStudio2Dev/macOS/stage/bin
scp cmbuild/lib/* repo:~/update2/FEBioStudio2Dev/macOS/stage/lib
ssh repo "chmod +x update2/FEBioStudio2Dev/macOS/stage/bin/febio4"
