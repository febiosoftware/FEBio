#! /bin/bash
set -e
TARGET_DIR="${TARGET_DIR:-febio4-sdk}"
mkdir -p ${TARGET_DIR}/{include,lib,bin}
# Gotta use rsync macos cp doesn't support --parent
find . -name "*.h*" -type f -not -path "./infrastructure/*" -not -path "./${TARGET_DIR}/*" -exec rsync -R -dir {} "./${TARGET_DIR}/include/" \;
find cmbuild/lib -name "*" -type f -exec cp -a {} "./${TARGET_DIR}/lib/" \;
find cmbuild/bin -name "*" -type f -exec cp -a {} "./${TARGET_DIR}/bin/" \;
