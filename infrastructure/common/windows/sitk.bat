call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
set SOURCE="https://github.com/SimpleITK/SimpleITK.git"
set BRANCH="v2.1.1.2"
git clone --depth 1 --branch "%BRANCH%" "%SOURCE%" "%BRANCH%"
pushd %BRANCH%
cmake .  -LA -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX=%INSTALLATION_PATH% ^
  -DWRAP_DEFAULT:BOOL=OFF ^
  -DBUILD_EXAMPLES:BOOL=OFF ^
  -DBUILD_TESTING:BOOL=OFF ^
  -DBUILD_SHARED_LIBS:BOOL=OFF ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DUSE_CCACHE=OFF ^
  -DWRAP_CSHARP=OFF ^
  -DWRAP_JAVA=OFF ^
  -DWRAP_LUA=OFF ^
  -DWRAP_PYTHON=OFF ^
  -DWRAP_R=OFF ^
  -DWRAP_RUBY=OFF ^
  -DWRAP_TCL=OFF
pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
REM RD /S /Q %BRANCH%
popd
