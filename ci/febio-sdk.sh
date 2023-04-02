#! /bin/bash
mkdir -p febio-sdk/{include,lib,bin}
find . -name "*.h*" -type f -exec cp --parents -a {} ./febio-sdk/include/ \;
find febio/lib -name "*" -type f -exec cp -a {} ./febio-sdk/lib/ \;
find febio/bin -name "*" -type f -exec cp -a {} ./febio-sdk/bin/ \;
