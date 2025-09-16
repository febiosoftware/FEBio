call "%ONEAPI_ROOT%\setvars.bat"

set CONFIG_TYPE="Release"

IF "%~1"=="-d" set CONFIG_TYPE="Debug"

cmake . -LA -B cmbuild ^
  -DUSE_MMG=ON ^
  -DUSE_LEVMAR=ON ^
  -DUSE_HYPER=ON

cd "cmbuild"
msbuild /p:configuration=Release /maxCpuCount:%NUMBER_OF_PROCESSORS% ALL_BUILD.vcxproj
msbuild /p:configuration=Debug /maxCpuCount:%NUMBER_OF_PROCESSORS% ALL_BUILD.vcxproj
cd ..

exit /b %errorlevel%
