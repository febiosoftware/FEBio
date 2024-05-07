call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
set SOURCE="https://github.com/NGSolve/netgen.git"
set BRANCH="v6.2.2307"

git clone --depth 1 --branch "%BRANCH%" "%SOURCE%" "%BRANCH%"
pushd %BRANCH%
git submodule update --init --recursive
cmake .  -LA -B cmbuild ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_INSTALL_PREFIX="c:/usr/local" ^
  -DUSE_CCACHE:BOOL=OFF ^
  -DUSE_CGNS:BOOL=OFF ^
  -DUSE_CSG:BOOL=ON ^
  -DUSE_GEOM2D:BOOL=ON ^
  -DUSE_GUI:BOOL=OFF ^
  -DUSE_INTERFACE:BOOL=ON ^
  -DUSE_INTERNAL_TCL:BOOL=OFF ^
  -DUSE_JPEG:BOOL=OFF ^
  -DUSE_MPEG:BOOL=OFF ^
  -DUSE_MPI:BOOL=OFF ^
  -DUSE_MPI4PY:BOOL=OFF ^
  -DUSE_NATIVE_ARCH:BOOL=OFF ^
  -DUSE_NUMA:BOOL=OFF ^
  -DUSE_OCC:BOOL=ON ^
  -DUSE_PYTHON:BOOL=OFF ^
  -DUSE_STLGEOM:BOOL=ON ^
  -DUSE_SUPERBUILD:BOOL=OFF ^
  -DENABLE_CPP_CORE_GUIDELINES_CHECK:BOOL=OFF ^
  -DENABLE_UNIT_TESTS:BOOL=OFF

pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
RD /S /Q %BRANCH%
popd
