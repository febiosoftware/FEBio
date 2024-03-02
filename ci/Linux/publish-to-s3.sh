#! /bin/bash
set -e
chmod +x artifacts/febio/bin/febio4
ci/common/linux/publish-to-s3.sh artifacts/febio4
ci/common/linux/publish-to-s3.sh artifacts/febio4-sdk
