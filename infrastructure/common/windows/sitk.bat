call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
set SOURCE="https://github.com/SimpleITK/SimpleITK.git"
set BRANCH="v2.1.1.2"
git clone --depth 1 --branch "%BRANCH%" "%SOURCE%" "%BRANCH%"
cd %BRANCH%
cmake .  -LA -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX="c:\local" ^
  -DWRAP_DEFAULT:BOOL=OFF ^
  -DBUILD_EXAMPLES:BOOL=OFF ^
  -DBUILD_TESTING:BOOL=OFF ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DSimpleITK_USE_ELASTIX:BOOL=ON ^
  -DSimpleITK_USE_SYSTEM_ITK:BOOL=OFF
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
cd %TEMP%
REM RD /S /Q %BRANCH%
