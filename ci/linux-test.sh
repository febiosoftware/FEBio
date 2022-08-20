#! /bin/bash
# Uncomment next line if not global on target machine
source "/opt/intel/oneapi/setvars.sh" --force
ls -lr cmbuild/bin/*
chmod +x cmbuild/bin/febio3
ci/test-suite/nightly.py
