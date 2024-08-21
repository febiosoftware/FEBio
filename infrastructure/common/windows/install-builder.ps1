aws s3 cp s3://febiosoftware/windows/installbuilder/installbuilder-enterprise-23.11.0-windows-x64-installer.exe .
$installBuilderPath = "C:\Program Files\InstallBuilder Enterprise 23.11.0"
Start-Process .\installbuilder-enterprise-23.11.0-windows-x64-installer.exe `
  -ArgumentList '--mode unattended', '--installer-language en', "--prefix ""${installBuilderPath}"""`
  -NoNewWindow `
  -Wait

$builderPath = "${installBuilderPath}/bin"

# Conditionally add builderpath to environment
$machinepath=[Environment]::GetEnvironmentVariable("Path", "MACHINE")
$paths = $machinepath -split ';'

if ($builderPath -notin $paths) {
  $paths = $paths + $builderPath
  $paths = $paths -join ';'
  [Environment]::SetEnvironmentVariable("INSTALL_BUILDER_ROOT", $installBuilderPath, "MACHINE")
  [Environment]::SetEnvironmentVariable("Path", $paths, "MACHINE")
}
