call "%ONEAPI_ROOT%setvars.bat" --force

set AQT_DISABLE_UPDATES=1

call pip install aqtinstall -v
call aqt list-qt windows desktop
call aqt install-qt --outputdir "%QT_INSTALL_DIR%" windows desktop 6.9.3 win64_msvc2022_64 -m all
