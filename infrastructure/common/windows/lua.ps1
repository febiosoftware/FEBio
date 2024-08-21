choco install lua53 -y --version 5.3.5
$luaPath = "C:\ProgramData\chocolatey\lib\lua53\tools"
New-Item -Path "${luaPath}\lua.exe" -ItemType SymbolicLink -Value "${luaPath}\lua53.exe"
