#! /bin/bash
mkdir -p febio-sdk/{include,lib,bin}
# Gotta use rsync macos cp doesn't support --parent
find . -name "*.h*" -type f -not -path "./infrastructure/*" -exec rsync -R -dir {} ./febio-sdk/include/ \;
find cmbuild/lib -name "*" -type f -exec cp -a {} ./febio-sdk/lib/ \;
find cmbuild/bin -name "*" -type f -exec cp -a {} ./febio-sdk/bin/ \;
