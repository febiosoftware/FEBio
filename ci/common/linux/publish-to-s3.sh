#! /bin/bash
SOURCE=$1
DIR=$(dirname "$SOURCE")
SOURCE=$(basename "$SOURCE")
TARGET_NAME=$SOURCE
pushd "$DIR" || exit 1
OS=${OS:-Linux}
PACKAGE="${TARGET_NAME}-${GIT_TAG}-${OS}.tar.gz"
BUCKET_NAME=${BUCKET_NAME:-febio-packages}
REF_NAME=${REF_NAME:-$(git rev-parse --abbrev-ref HEAD)}
PACKAGE_NAME=${PACKAGE_NAME:-febio}
tar -czf "$PACKAGE" "$SOURCE"

FEBIO_VERSION="${FEBIO_VERSION:-v}"
PREFIX="${PACKAGE_NAME}/${REF_NAME}/${OS}"

# Fetch previous version so we can mark it for expiration
echo "Fetching previous for expiration"
PREVIOUS_PACKAGE=$(aws --output json s3api list-objects \
	--bucket "$BUCKET_NAME" \
	--prefix "${PREFIX}/" \
	--query "reverse(sort_by(Contents,&LastModified)) && Contents[?contains(Key, '$TARGET_NAME-${FEBIO_VERSION}')]" \
	| jq -r -e '. | max_by(.LastModified) | .Key')

# Upload new package
TARGET="s3://${BUCKET_NAME}/${PREFIX}/${PACKAGE}"
aws s3 cp "$PACKAGE" "$TARGET"

# aws query fails if prefix doesn't exist and value is set to empty string ''
# jq sets empty value to null, hence the test for null
if [[ ! -z "${PREVIOUS_PACKAGE}" ]] && [[ $PREVIOUS_PACKAGE != null ]]; then
	echo "Tagging previous for expiration"
	aws s3api put-object-tagging \
	--bucket $BUCKET_NAME \
	--key $PREVIOUS_PACKAGE \
	--tagging '{"TagSet": [{ "Key": "ttl", "Value": "expire" }]}'
fi


popd || exit 1
