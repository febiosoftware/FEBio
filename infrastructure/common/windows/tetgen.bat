call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
SET TETGEN="https://github.com/libigl/tetgen.git"
SET TETGEN_SRC="tetgen"
git clone %TETGEN%
cd "%TETGEN_SRC%"
cmake .  -LA -B cmbuild -DCMAKE_INSTALL_PREFIX="c:\local"
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% ALL_BUILD.vcxproj
cd ..
xcopy cmbuild\Release\tetgen.lib "c:\local\lib" /E /C /H /R /K /O /Y
xcopy tetgen.h "c:\local\include" /E /C /H /R /K /O /Y
cd %TEMP%
RD /S /Q %TETGEN_SRC%
