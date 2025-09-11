#! /bin/bash
TARGET_DIR="${TARGET_DIR:-febio4-sdk}"
mkdir -p ${TARGET_DIR}/{include,lib,bin}

find . \( -name "*.h" -o -name "*.hpp" \) -type f -not -path "./${TARGET_DIR}/*" -exec cp --parents -a -t "./${TARGET_DIR}/include/" {} +
find cmbuild/lib -type f -exec cp -a -t "./${TARGET_DIR}/lib/" {} +
find cmbuild/bin -type f -exec cp -a -t "./${TARGET_DIR}/bin/" {} +
