set PYTHON="%ONEAPI_ROOT%intelpython\latest\python"
copy "%ONEAPI_ROOT%\compiler\latest\windows\redist\intel64_win\compiler\libiomp5md.dll" cmbuild\bin\Release
%PYTHON% ./TestSuite/code/tools.py -r ./cmbuild/bin/Release/febio4.exe
