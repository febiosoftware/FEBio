$vcpkgPackagePath = $Env:VCPKG_PACKAGE_PATH
cd $vcpkgPackagePath
vcpkg install --triplet=x64-windows
$machinepath=[Environment]::GetEnvironmentVariable("Path", "MACHINE")
$paths = $machinepath -split ';'
$paths = $paths + "${vcpkgPackagePath}\vcpkg_installed\x64-windows"
$value = $paths -join ';'
[Environment]::SetEnvironmentVariable("Path", $value, "MACHINE")
