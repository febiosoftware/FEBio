call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
git clone https://github.com/MmgTools/mmg.git
cd mmg
cmake . -LA -B cmbuild -DCMAKE_INSTALL_PREFIX="c:\local"
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
cd %TEMP%
RD /S /Q mmg
