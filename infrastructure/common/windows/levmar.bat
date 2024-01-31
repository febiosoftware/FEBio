call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
git clone https://github.com/jturney/levmar.git
pushd levmar
cmake . -LA -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX=%INSTALLATION_PATH% ^
  -DBUILD_DEMO:BOOLEAN=false
pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
RD /S /Q levmar
popd
