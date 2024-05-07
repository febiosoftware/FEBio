call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
git clone --depth 1 --branch v1.10.1 https://github.com/nih-at/libzip.git
pushd libzip
cmake . -LA -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX=%INSTALLATION_PATH% ^
  -DBUILD_DOC=OFF ^
  -DBUILD_EXAMPLES=OFF ^
  -DBUILD_OSSFUZZ=OFF ^
  -DBUILD_REGRESS=OFF ^
  -DBUILD_TOOLS=OFF ^
  -DBUILD_SHARED_LIBS=ON ^
  -DENABLE_BZIP2=OFF ^
  -DENABLE_COMMONCRYPTO=OFF ^
  -DENABLE_FDOPEN=OFF ^
  -DENABLE_GNUTILS=OFF ^
  -DENABLE_LZMA=OFF ^
  -DENABLE_MBEDTLS=OFF ^
  -DENABLE_OPENSSL=OFF ^
  -DENABLE_WINDOWS_CRYPTO=OFF ^
  -DLIBZIP_DO_INSTALL=ON ^
  -DSHARED_LIB_VERSIONNING=ON
pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
RD /S /Q libzip
popd
