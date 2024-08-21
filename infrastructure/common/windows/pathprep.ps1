$installPath = $Env:INSTALLATION_PATH
$srcPath = $Env:SOURCE_PATH
$vcpkgPackagePath = $Env:VCPKG_PACKAGE_PATH
md $installPath
md $srcPath
md $vcpkgPackagePath
$systempath=[Environment]::GetEnvironmentVariable("Path", "Machine")
$systempath = $systempath -split ';'
$systempath= $systempath + $installPath
$systempath = $systempath  -join ';'
[Environment]::SetEnvironmentVariable("Path", $systempath, "MACHINE")

[Environment]::GetEnvironmentVariable("Path", "Machine") -Split ';'
