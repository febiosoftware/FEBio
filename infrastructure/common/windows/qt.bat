call "%ONEAPI_ROOT%setvars.bat" --force
call pip install aqtinstall -v
call aqt list-qt windows desktop
call aqt install-qt --outputdir "%QT_INSTALL_DIR%" windows desktop 6.7.3 win64_msvc2019_64 -m all
