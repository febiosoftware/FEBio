#! /bin/bash
scp cmbuild/bin/Release/* repo:~/update2/FEBioStudio2Dev/Windows/stage/bin

if [ -f testLogs/Logs/* ]; then
    scp testLogs/Logs/* repo:~/TestSuite/Logs/windows.txt
fi