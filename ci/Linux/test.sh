#! /bin/bash
# Uncomment next line if not global on target machine
set -e
source "/opt/intel/oneapi/setvars.sh" --force
ls -lr cmbuild
chmod +x cmbuild/bin/febio4
ci/test-suite/nightly.py
