#! /bin/bash

chmod +x artifacts/febio/bin/febio4
ci/common/linux/publish-to-s3.sh artifacts/febio
ci/common/linux/publish-to-s3.sh artifacts/febio-sdk
