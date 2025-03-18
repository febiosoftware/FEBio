#! /bin/bash
OS="${OS}"
BUCKET="${BUCKET:-febio-packages}"
PACKAGE_PATH="${PACKAGE_PATH:-febio4}"
PACKAGE="${PACKAGE:-febio4-sdk}"
PACKAGE_BRANCH="${PACKAGE_BRANCH:-develop}"
VERSION="${VERSION:-v}"

PACKAGE_PREFIX="${PACKAGE_PATH}/${PACKAGE_BRANCH}/${OS}/"
PACKAGE_SEARCH="${PACKAGE}-${VERSION}"

PACKAGE_URI=$(aws --output json s3api list-objects \
	--bucket "$BUCKET" \
	--prefix "$PACKAGE_PREFIX" \
	--query "reverse(sort_by(Contents,&LastModified)) && Contents[?contains(Key, '$PACKAGE_SEARCH')]" \
	| jq -r -e '. | max_by(.LastModified) | .Key')

STATUS=$?

if [[ ! -z "${PACKAGE_URI}" ]] && [[ $PACKAGE_URI != null ]]; then
	echo "SDK found at ${PACKAGE_URI}"
	ARCHIVE="${PACKAGE_URI##*/}"
	aws s3 cp "s3://$BUCKET/$PACKAGE_URI" .
	tar xvzf "$ARCHIVE"
else
	echo "SDK not found at ${PACKAGE_PREFIX}${PACKAGE_SEARCH}*"; exit $STATUS
fi
