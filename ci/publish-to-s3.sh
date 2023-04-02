#! /bin/bash
set -x
SOURCE=$1
GIT_TAG=${GIT_TAG}
PACKAGE="${SOURCE}-${GIT_TAG}.tgz"
BUCKET_NAME=${BUCKET_NAME:-febio-packages}
REF_NAME=${REF_NAME:-$(git rev-parse --abbrev-ref HEAD)}
PACKAGE_NAME=${PACKAGE_NAME:-febio}
OS=${OS:-Linux}
tar -czf $PACKAGE $SOURCE
TARGET="s3://${BUCKET_NAME}/${PACKAGE_NAME}/${REF_NAME}/${OS}/${PACKAGE}"
aws s3 cp $PACKAGE $TARGET
