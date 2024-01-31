call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
git clone https://github.com/MmgTools/mmg.git
pushd mmg
cmake . -LA -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX=%INSTALLATION_PATH%
pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
RD /S /Q mmg
popd
