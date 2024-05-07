call "%ONEAPI_ROOT%setvars.bat" --force
cd %TEMP%
set ZLIB=zlib-1.2.12.tar.gz
set ZLIB_HOST=http://zlib.net
set ZLIB_SRC=zlib-1.2.12
set ZLIB_URI=%ZLIB_HOST%/%ZLIB%
curl -O -L %ZLIB_URI%
bash --login -c "tar xzf %ZLIB%"
cd %ZLIB_SRC%
cmake .  -LA -B cmbuild -DCMAKE_INSTALL_PREFIX="c:\local"
cd cmbuild
msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
cd %TEMP%
DEL %ZLIB%
RD /S /Q %ZLIB_SRC%
