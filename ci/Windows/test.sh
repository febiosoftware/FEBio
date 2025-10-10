#! /usr/bin/bash
# Uncomment next line if not global on target machine
set -e
set -x
FEBIO_XML=$(realpath ./ci/febio.xml)
PLUGIN_DIR=$(realpath ./plugins)
FEBIO_DIR=$(realpath ./cmbuild/bin/Release)
FEBIO_LIB=$(realpath ./cmbuild/bin/Release)
FEBIO_BIN="${FEBIO_DIR}/febio4.exe"

# Convert plugin dir to Windows style
PLUGIN_DIR_WIN=$(cygpath -w "$PLUGIN_DIR" | sed 's|\\|/|g')

# Copy the plugins from their subdirectories directly to
# the root of the plugin dir
cp $PLUGIN_DIR/*/*.dll $PLUGIN_DIR

# Set the plugin dir in the FEBio XML
sed -i "s@PLUGINS_FOLDER@${PLUGIN_DIR_WIN}@g" "$FEBIO_XML"

# Copy febio xml into febio dir
cp $FEBIO_XML $FEBIO_DIR

# Copy iomp lib into febio dir
ONEAPI=$(cygpath -u "$ONEAPI_ROOT") #Convert $ONEAPI_ROOT to posix path
IOMP_LIB="${ONEAPI}compiler/latest/windows/redist/intel64_win/compiler/libiomp5md.dll"
cp -a "$IOMP_LIB" "$FEBIO_LIB"

# Copy fftw lib into febio dir
FFTW_LIB="/c/vcpkg/packages/fftw3_x64-windows/bin/fftw3.dll"
cp -a "$FFTW_LIB" "$FEBIO_LIB"

# Run the test suite
PYTHON="${ONEAPI}intelpython/latest/python"
"$PYTHON" ./TestSuite/code/tools.py -r $FEBIO_BIN -n
