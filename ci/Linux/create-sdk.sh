#! /bin/bash
set -e
TARGET_DIR="${TARGET_DIR:-febio4-sdk}"
mkdir -p ${TARGET_DIR}/{include,lib,bin}

rsync -a --prune-empty-dirs --include="*.h" --include="*.hpp" --exclude="*" ./ "./${TARGET_DIR}/include/"
find cmbuild/lib -type f -exec cp -at "./${TARGET_DIR}/lib/" {} +
find cmbuild/bin -type f -exec cp -at "./${TARGET_DIR}/bin/" {} +

