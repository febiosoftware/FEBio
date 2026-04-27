#! /usr/bin/bash
set -e

FEBIO_XML=$(realpath ./ci/febio.xml)
FEBIO_DIR=$(realpath ./cmbuild/bin/Release)
FEBIO_LIB=$(realpath ./cmbuild/bin/Release)
FEBIO_BIN="${FEBIO_DIR}/febio4.exe"

# If PLUGIN_DIRS is not empty, the plugins were built and we need to copy them 
# and set the plugin folder in the febio xml. 
if [[ -n "$PLUGIN_DIRS" ]]; then
    PLUGIN_DIR=$(realpath ./plugins)

    # Convert plugin dir to Windows style
    PLUGIN_DIR_WIN=$(cygpath -w "$PLUGIN_DIR" | sed 's|\\|/|g')

    # Copy the plugins from their subdirectories directly to
    # the root of the plugin dir
    cp $PLUGIN_DIR/*/*.dll $PLUGIN_DIR

    # Set the plugin dir in the FEBio XML
    sed -i "s@PLUGINS_FOLDER@${PLUGIN_DIR_WIN}@g" "$FEBIO_XML"

# Otherwise, remove the import_folder tag from the febio xml or FEBio will fail to start
else
    sed -i 's@<import_folder>PLUGINS_FOLDER</import_folder>@@' $FEBIO_XML
fi

# Copy febio xml into febio dir
cp $FEBIO_XML $FEBIO_DIR

# Copy iomp lib into febio dir
ONEAPI=$(cygpath -u "$ONEAPI_ROOT") #Convert $ONEAPI_ROOT to posix path
IOMP_LIB="${ONEAPI}compiler/latest/windows/redist/intel64_win/compiler/libiomp5md.dll"
cp -a "$IOMP_LIB" "$FEBIO_LIB"

# Copy fftw lib into febio dir
FFTW_LIB="/c/usr/local/febio/vcpkg_installed/x64-windows/bin/fftw3.dll"
cp -a "$FFTW_LIB" "$FEBIO_LIB"

ZLIB="/c/usr/local/febio/vcpkg_installed/x64-windows/bin/zlib1.dll"
cp -a "$ZLIB" "$FEBIO_LIB"

# Run the test suite
python ./TestSuite/code/tools.py -r $FEBIO_BIN -n -p $PLUGIN_DIRS
