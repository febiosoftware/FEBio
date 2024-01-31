$vcpkgInstallPath="c:\vcpkg"
git clone https://github.com/MicroSoft/vcpkg.git $vcpkgInstallPath
& "$vcpkgInstallPath\bootstrap-vcpkg.bat"

# Add vcpkg_root to environment
[Environment]::SetEnvironmentVariable("VCPKG_ROOT", $vcpkgInstallPath, "MACHINE")
$machinepath=[Environment]::GetEnvironmentVariable("Path", "MACHINE")

# Append Release build only to x64-windows triplet set(VCPKG_BUILD_TYPE release)
$tripletPath = "${vcpkgInstallPath}\triplets\x64-windows.cmake"
Get-Content $tripletPath
Add-Content -Path $tripletPath "set(VCPKG_BUILD_TYPE release)"
Get-Content $tripletPath

# Add vcpkg to machine path
$paths = $machinepath -split ';'
$paths = $paths + $vcpkgInstallPath
$value = $paths -join ';'
[Environment]::SetEnvironmentVariable("Path", $value, "MACHINE")
