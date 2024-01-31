call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
set SOURCE="https://github.com/InsightSoftwareConsortium/ITK.git"
set BRANCH="v5.2.1"
git clone --depth 1 --branch "%BRANCH%" "%SOURCE%" "%BRANCH%"
pushd %BRANCH%
cmake .  -LA -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX=%INSTALLATION_PATH% ^
  -DBUILD_EXAMPLES:BOOL=OFF ^
  -DBUILD_SHARED_LIBS:BOOL=OFF ^
  -DBUILD_TESTING:BOOL=OFF ^
  -DITK_USE_SYSTEM_EIGEN:BOOL=ON ^
  -DITK_WRAP_PYTHON:BOOL=OFF ^
  -DITK_DOXYGEN_HTML:BOOL=OFF ^
  -DCMAKE_BUILD_TYPE=Release
pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
RD /S /Q %BRANCH%
popd
