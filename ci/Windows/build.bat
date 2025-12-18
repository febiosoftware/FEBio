call "%ONEAPI_ROOT%\setvars.bat"

set CONFIG_TYPE="Release"

IF "%~1"=="-d" set CONFIG_TYPE="Debug"

cmake . -LA -B cmbuild ^
  -DSET_DEVCOMMIT=ON ^
  -DUSE_FFTW=ON ^
  -DUSE_HYPRE=ON ^
  -DUSE_LEVMAR=ON ^
  -DUSE_MKL=ON ^
  -DUSE_MMG=ON ^
  -DUSE_ZLIB=ON ^
  -DFFTW_LIB="C:\usr\local\febio\vcpkg_installed\x64-windows\lib" ^
  -DHYPRE_LIB="C:\usr\local\lib" ^
  -DMMG_LIB="C:\usr\local\lib" ^
  -DLEVMAR_INC="C:\usr\local\include\levmar" ^
  -DLEVMAR_LIB="C:\usr\local\lib"
  

cd "cmbuild"
msbuild /p:configuration=Release /maxCpuCount:%NUMBER_OF_PROCESSORS% ALL_BUILD.vcxproj
msbuild /p:configuration=Debug /maxCpuCount:%NUMBER_OF_PROCESSORS% ALL_BUILD.vcxproj
cd ..

exit /b %errorlevel%
