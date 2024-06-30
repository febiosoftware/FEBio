call "%ONEAPI_ROOT%setvars.bat" --force
call pip install aqtinstall -v
call aqt list-qt windows desktop
call aqt install-qt --outputdir "%QT_INSTALL_DIR%" windowsdesktop 6.7.2 win64_msvc2019_64 -m all
