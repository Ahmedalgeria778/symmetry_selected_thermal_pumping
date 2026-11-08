# Idempotent launcher for the pending DNS campaign D.3/D.4/D.5.
# Restarts are safe: the runner resumes from the checkpoint CSV and skips
# every case_id already present there.
$ErrorActionPreference = 'Stop'
$root = 'C:\Users\user\Documents\Documents - Copie\mixed_conv_verification'
$dir = Join-Path $root 'campaign_v2_symmetry_selection'
$py = 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe'
$script = Join-Path $dir 'run_dns_D3_D4_D5.py'
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$out = Join-Path $dir ('log_DNS_D345_' + $stamp + '.txt')
$err = Join-Path $dir ('log_DNS_D345_' + $stamp + '.err.txt')

$alive = Get-CimInstance Win32_Process -Filter "Name like '%python%'" |
         Where-Object { $_.CommandLine -match 'run_dns_D3_D4_D5' }
if ($alive) {
    Write-Output ('already running: pid ' + (($alive.ProcessId | Sort-Object) -join ','))
    exit 0
}

$p = Start-Process -FilePath $py `
     -ArgumentList @('-u', '-X', 'utf8', ('"' + $script + '"'), '--go') `
     -WorkingDirectory $root -WindowStyle Hidden `
     -RedirectStandardOutput $out -RedirectStandardError $err -PassThru
Write-Output ('launched pid ' + $p.Id + ' log ' + (Split-Path $out -Leaf))
