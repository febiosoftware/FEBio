#! /bin/bash
scp cmbuild/bin/* repo:~/update2/FEBioStudio2Dev/macOS/stage/FEBioStudio.app/Contents/MacOS
scp cmbuild/lib/* repo:~/update2/FEBioStudio2Dev/macOS/stage/FEBioStudio.app/Contents/Frameworks
ssh repo "chmod +x update2/FEBioStudio2Dev/macOS/stage/FEBioStudio.app/Contents/MacOS/febio4"
