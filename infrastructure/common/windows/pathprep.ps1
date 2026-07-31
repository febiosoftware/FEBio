$installPath = $Env:INSTALLATION_PATH
$srcPath = $Env:SOURCE_PATH
$vcpkgPackagePath = $Env:VCPKG_PACKAGE_PATH
md $installPath
md $srcPath
md $vcpkgPackagePath

$gitPath = "C:\Program Files\Git\bin"
$gnu32Path = "C:\Program Files\Git\usr\bin"

$systempath=[Environment]::GetEnvironmentVariable("Path", "Machine")
$systempath = $systempath -split ';'
$systempath= $systempath + $installPath + $gitPath + $gnu32Path
$systempath = $systempath  -join ';'
[Environment]::SetEnvironmentVariable("Path", $systempath, "MACHINE")

[Environment]::GetEnvironmentVariable("Path", "Machine") -Split ';'
