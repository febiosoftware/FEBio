#! /bin/bash
mkdir -p febio-sdk/{include,lib,bin}
find . -name "*.h*" -type f -exec cp --parents -a {} ./febio-sdk/include/ \;
find cmbuild/lib -name "*" -type f -exec cp -a {} ./febio-sdk/lib/ \;
find cmbuild/bin -name "*" -type f -exec cp -a {} ./febio-sdk/bin/ \;
