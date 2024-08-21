Set-ExecutionPolicy Bypass -Scope Process -Force; [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072; iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))

choco install git -y
# Set Path for git-bash and perl
$gitPath = "C:\Program Files\Git\bin"
$gnu32Path = "C:\Program Files\Git\usr\bin"
$machinePath = [Environment]::GetEnvironmentVariable("Path", "MACHINE")
$machinePath = $machinePath -Split ';'
$machinePath + $machinePath + $gitPath + $gnu32Path
$machinePath = $machinePath -Join ';'
[Environment]::SetEnvironmentVariable("Path", $machinePath, "MACHINE")
