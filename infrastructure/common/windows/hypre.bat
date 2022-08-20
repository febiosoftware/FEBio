call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
git clone https://github.com/hypre-space/hypre.git
cd hypre/src
cmake .  -LA -B cmbuild -DCMAKE_INSTALL_PREFIX="c:\local" -DHYPRE_HAVE_MPI=Off -DHYPRE_WITH_MPI=Off
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
cd %TEMP%
RD /S /Q hypre
