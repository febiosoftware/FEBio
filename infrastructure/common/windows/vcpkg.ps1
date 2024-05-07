$vcpkgInstallPath="C:\usr\local\vcpkg"
md $vcpkgInstallPath
git clone https://github.com/MicroSoft/vcpkg.git $vcpkgInstallPath
& "$vcpkgInstallPath\bootstrap-vcpkg.bat"
[Environment]::SetEnvironmentVariable("VCPKG_ROOT", $vcpkgInstallPath, "MACHINE")
$machinepath=[Environment]::GetEnvironmentVariable("Path", "MACHINE")

$paths = $machinepath -split ';'
$paths = $paths + $vcpkgInstallPath
$value = $paths -join ';'
[Environment]::SetEnvironmentVariable("Path", $value, "MACHINE")
