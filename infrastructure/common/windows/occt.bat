call "%ONEAPI_ROOT%setvars.bat"
REM cd %TEMP%
set OCCT=opencascade-7.7.0.tgz
set OCCT_HOST=https://dev.opencascade.org/system/files/occt/OCC_7.7.0_release
set OCCT_SRC=opencascade-7.7.0.tgz
set OCCT_URI=%OCCT_HOST%/%OCCT%
echo "curl -O -L %OCCT_URI%"
REM bash --login -c "tar xzf %OCCT%"
REM cd %OCCT_SRC%
REM cmake .  -LA -B cmbuild -DCMAKE_INSTALL_PREFIX="c:\local"
REM cd cmbuild
REM msbuild /P:Configuration=Release /m:%NUMBER_OF_PROCESSORS% INSTALL.vcxproj
REM cd %TEMP%
REM DEL %OCCT%
REM RD /S /Q %OCCT_SRC%
