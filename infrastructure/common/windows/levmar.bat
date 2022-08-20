call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
git clone https://github.com/jturney/levmar.git
cd levmar
cmake . -LA -B cmbuild -DCMAKE_INSTALL_PREFIX="c:\local" -DBUILD_DEMO:BOOLEAN=false
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
cd %TEMP%
RD /S /Q levmar
