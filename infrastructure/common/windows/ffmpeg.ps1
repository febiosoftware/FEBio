choco install ffmpeg-shared -y --version 6.1.1
[Environment]::SetEnvironmentVariable("VCPKG_ROOT", $vcpkgInstallPath, "MACHINE")
# Set path for ffmpeg
$machinepath=[Environment]::GetEnvironmentVariable("Path", "MACHINE")
$ffmpegPath = "C:\Program Files\FFMpeg"
$paths = $machinepath -split ';'
$paths = $paths + $ffmpegPath
$paths = $paths -join ';'
[Environment]::SetEnvironmentVariable("Path", $paths, "MACHINE")
