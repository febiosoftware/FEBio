#! /bin/bash
# Uncomment next line if not global on target machine
source "/opt/intel/oneapi/setvars.sh" --force
ls -lr febio
chmod +x cmbuild/bin/febio4
ci/test-suite/nightly.py
