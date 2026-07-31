#! /bin/bash
# Uncomment next line if not global on target machine
set -e

FEBIO_XML=$(realpath ./ci/febio.xml)
FEBIO_DIR=$(realpath ./cmbuild/bin)
FEBIO_BIN="${FEBIO_DIR}/febio4"
chmod +x $FEBIO_BIN

# If PLUGIN_DIRS is not empty, the plugins were built and we need to copy them 
# and set the plugin folder in the febio xml. 
if [[ -n "$PLUGIN_DIRS" ]]; then
    PLUGIN_DIR=$(realpath ./plugins)

    # Copy the plugins from their subdirectories directly to
    # the root of the plugin dir
    cp $PLUGIN_DIR/*/*.dylib $PLUGIN_DIR 

    # Set the plugin dir in the FEBio XML
    sed -i '' 's@PLUGINS_FOLDER@'$PLUGIN_DIR'@' $FEBIO_XML

# Otherwise, remove the import_folder tag from the febio xml or FEBio will fail to start
else
    sed -i '' 's@<import_folder>PLUGINS_FOLDER</import_folder>@@' $FEBIO_XML
fi

# Copy febio xml into febio dir
cp $FEBIO_XML $FEBIO_DIR

# Run the test suite
./TestSuite/code/tools.py -r $FEBIO_BIN -n -p $PLUGIN_DIRS

# remove the febio.xml file to prevent it from getting uploaded to the repo server.
rm $FEBIO_DIR/febio.xml 