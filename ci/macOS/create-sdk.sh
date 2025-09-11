#! /bin/bash
TARGET_DIR="${TARGET_DIR:-febio4-sdk}"
mkdir -p ${TARGET_DIR}/{include,lib,bin}

rsync -a --prune-empty-dirs --include="*/" --include="*.h" --include="*.hpp" --exclude="*" ./ "./${TARGET_DIR}/include/"
find cmbuild/lib -name "*" -type f -exec cp -a {} "./${TARGET_DIR}/lib/" \;
find cmbuild/bin -name "*" -type f -exec cp -a {} "./${TARGET_DIR}/bin/" \;
