#! /bin/bash
TARGET_DIR="${TARGET_DIR:-febio4-sdk}"
mkdir -p ${TARGET_DIR}/{include,lib,bin}
find . -name "*.h*" -type f -not -path "./infrastructure/*" -not -path "./${TARGET_DIR}/*" -exec cp --parents -a {} "./${TARGET_DIR}/include/" \;
find cmbuild/lib -name "*" -type f -exec cp -a {} "./${TARGET_DIR}/lib/" \;
find cmbuild/bin -name "*" -type f -exec cp -a {} "./${TARGET_DIR}/bin/" \;
