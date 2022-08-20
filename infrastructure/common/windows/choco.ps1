Set-ExecutionPolicy Bypass -Scope Process -Force; [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072; iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))

choco install git -y
# Set Path for git-bash
$userpath=[Environment]::GetEnvironmentVariable("Path", "User")
setx PATH "$userpath;C:\Program Files\Git\bin"
