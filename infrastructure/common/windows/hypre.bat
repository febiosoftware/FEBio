call "%ONEAPI_ROOT%setvars.bat" --force
pushd %SOURCE_PATH%
git clone https://github.com/hypre-space/hypre.git
pushd hypre\src
cmake .  -L -B cmbuild ^
  -DCMAKE_INSTALL_PREFIX=%INSTALLATION_PATH% ^
  -DHYPRE_HAVE_MPI=Off ^
  -DHYPRE_WITH_MPI=Off

pushd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
popd
popd
RD /S /Q hypre
popd
