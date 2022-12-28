call "%ONEAPI_ROOT%\setvars.bat"
cmake . -LA -B "cmbuild
cd "cmbuild
msbuild /p:configuration=Release /maxCpuCount:%NUMBER_OF_PROCESSORS% ALL_BUILD.vcxproj
copy bin\febio.xml bin\Release\
cd ..
