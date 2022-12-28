call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
set SOURCE="https://github.com/InsightSoftwareConsortium/ITK.git"
set BRANCH="v5.2.1"
git clone --depth 1 --branch "%BRANCH%" "%SOURCE%" "%BRANCH%"
cd %BRANCH%
cmake .  -LA -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX="c:\local"
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
cd %TEMP%
REM RD /S /Q %BRANCH%
