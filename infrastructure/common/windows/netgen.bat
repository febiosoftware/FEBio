call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
set SOURCE="https://github.com/NGSolve/netgen.git"
set BRANCH="v6.2.2204"
REM git clone --depth 1 --branch "%BRANCH%" "%SOURCE%" "%BRANCH%"
cd %BRANCH%
git submodule update --init --recursive
cmake --build .  ^
  --config Release ^
  --target install ^
  -DUSE_SUPERBUILD=ON ^
  -DCMAKE_INSTALL_PREFIX="c:\local\netgen" ^
  -DUSE_PYTHON=OFF ^
  -DUSE_GUI=OFF ^
  -DUSE_NATIVE_ARCH=OFF ^
  -DUSE_MPEG=ON ^
  -DUSE_OCC=ON 
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
cd %TEMP%
REM RD /S /Q %BRANCH%
