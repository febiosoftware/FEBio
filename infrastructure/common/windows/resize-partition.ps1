Get-Partition
"rescan" | diskpart
$currentSize = $(Get-Partition -DriveLetter C).Size
$maxSize = $(Get-PartitionSupportedSize -DriveLetter C).SizeMax
if($currentSize -lt $maxSize) {
    "Resizing Partition"
    Resize-Partition -DriveLetter C -Size $maxSize
    Get-Partition
}

# TODO: Move this to choco.ps1
$machinePath = [Environment]::GetEnvironmentVariable("Path", "MACHINE")
$machinePath = $machinePath -Split ';'
$gnu32 = "C:\Program Files\Git\usr\bin"
if($machinePath -notcontains $gnu32) {
  $machinePath = $machinePath + $gnu32
  $machinePath = $machinePath -join ';'
  [Environment]::SetEnvironmentVariable("Path", $machinePath, "MACHINE")
}
