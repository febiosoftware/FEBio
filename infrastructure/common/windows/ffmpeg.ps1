choco install ffmpeg-shared -y --version 6.1.1

# Conditionally set path for ffmpeg
$machinepath=[Environment]::GetEnvironmentVariable("Path", "MACHINE")
$ffmpegPath = "C:\Program Files\FFMpeg"
$paths = $machinepath -split ';'

if ($ffmpegPath -notin $paths) {
  $paths = $paths + $ffmpegPath
  $paths = $paths -join ';'
  [Environment]::SetEnvironmentVariable("Path", $paths, "MACHINE")
}
