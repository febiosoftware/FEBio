Set-PSDebug -Trace 1
$uri="https://www.python.org/ftp/python/3.11.1/python-3.11.1-amd64.exe" 
$installer="python-3.11.1-amd64.exe"
$arguments="/quiet InstallAllUsers=0 PrependPath=1 Include_test=0 Include_doc=0"
Invoke-WebRequest -Uri $uri -OutFile $installer
$process=Start-Process -NoNewWindow -passthru .\$installer $arguments
Write-Verbose 'Sending ENTER keystrokes until the window closes...'
while (-not $process.HasExited) {
  # To be safe, activate the external program's window. If that fails, it must be closed already.
  try { [Microsoft.VisualBasic.Interaction]::AppActivate($process.Id) } catch { break }
  # Send the keystroke.
  [System.Windows.Forms.SendKeys]::SendWait('{Enter}')
  Start-Sleep -Milliseconds 200  # Sleep a little between attempts.
}
$env:Path = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
python -V
Exit
