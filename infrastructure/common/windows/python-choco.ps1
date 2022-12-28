Set-PSDebug -Trace 1
choco feature enable -n allowGlobalConfirmation
choco install vcredist140 --x86 -y --force --no-progress --reduce-package-size
choco install python3 -y --force
$env:Path = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
python -V
