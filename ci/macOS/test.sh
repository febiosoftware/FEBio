#! /bin/bash
# Uncomment next line if not global on target machine
set -e

FEBIO_XML=$(realpath ./ci/febio.xml)
PLUGIN_DIR=$(realpath ./plugins)
FEBIO_DIR=$(realpath ./cmbuild/bin)
FEBIO_LIB=$(realpath ./cmbuild/lib)
FEBIO_BIN="${FEBIO_DIR}/febio4"
chmod +x $FEBIO_BIN

# Copy the plugins from their subdirectories directly to
# the root of the plugin dir
cp $PLUGIN_DIR/*/*.dylib $PLUGIN_DIR

# Set the plugin dir in the FEBio XML
sed -i '' "s@PLUGINS_FOLDER@${PLUGIN_DIR}@g" "$FEBIO_XML"

# Copy febio xml into febio dir
cp $FEBIO_XML $FEBIO_DIR

# Run the test suite
./TestSuite/code/tools.py -r $FEBIO_BIN -n

