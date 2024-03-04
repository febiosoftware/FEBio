call "%ONEAPI_ROOT%\setvars.bat"

set CONFIG_TYPE="Release"

IF "%~1"=="-d" set CONFIG_TYPE="Debug"

cmake . -LA -B "cmbuild
cd "cmbuild
msbuild /p:configuration=%CONFIG_TYPE% /maxCpuCount:%NUMBER_OF_PROCESSORS% ALL_BUILD.vcxproj
copy bin\febio.xml bin\%CONFIG_TYPE%\
cd ..
