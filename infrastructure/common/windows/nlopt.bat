call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
set SOURCE="https://github.com/stevengj/nlopt.git"
set BRANCH="v2.10.1"

git clone --depth 1 --branch "%BRANCH%" "%SOURCE%" "%BRANCH%"
pushd %BRANCH%
cmake . -LA -B cmbuild ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_INSTALL_PREFIX="c:/usr/local" ^
  -DBUILD_SHARED_LIBS=OFF ^	
	-DNLOPT_GUILE=OFF ^
	-DNLOPT_JAVA=OFF ^
	-DNLOPT_MATLAB=OFF ^
	-DNLOPT_OCTAVE=OFF ^
	-DNLOPT_PYTHON=OFF ^
	-DNLOPT_PYTHON_SABI=OFF ^
	-DNLOPT_SWIG=OFF ^
	-DNLOPT_TESTS=OFF 

pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
RD /S /Q %BRANCH%
popd
