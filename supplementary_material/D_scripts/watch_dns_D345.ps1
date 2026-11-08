# Supervisor for the DNS campaign D.3/D.4/D.5.
# Polls every 15 minutes: if the runner process is gone while the checkpoint
# still has fewer than the expected rows, it relaunches (resume is automatic).
# It never kills a live process, because a single case legitimately integrates
# for hours without writing to its log.
$ErrorActionPreference = 'Continue'
$root = 'C:\Users\user\Documents\Documents - Copie\mixed_conv_verification'
$dir = Join-Path $root 'campaign_v2_symmetry_selection'
$ckpt = Join-Path $dir 'DNS_D345_results.csv'
$supLog = Join-Path $dir 'log_DNS_D345_supervisor.txt'
$launch = Join-Path $dir 'launch_dns_D345.ps1'
$target = 13
$maxRestarts = 10
$deadline = (Get-Date).AddHours(96)
$restarts = 0

while ($true) {
    $rows = 0
    if (Test-Path $ckpt) { $rows = @(Import-Csv $ckpt).Count }
    $alive = Get-CimInstance Win32_Process -Filter "Name like '%python%'" |
             Where-Object { $_.CommandLine -match 'run_dns_D3_D4_D5' }
    $ts = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'

    if ($rows -ge $target) {
        Add-Content -Path $supLog -Value ($ts + ' COMPLETE ' + $rows + '/' + $target)
        break
    }
    if ((Get-Date) -gt $deadline) {
        Add-Content -Path $supLog -Value ($ts + ' DEADLINE reached at ' + $rows + '/' + $target)
        break
    }
    if ($alive) {
        Add-Content -Path $supLog -Value ($ts + ' alive pid ' +
            (($alive.ProcessId | Sort-Object) -join ',') + ' rows ' + $rows + '/' + $target)
    } else {
        if ($restarts -ge $maxRestarts) {
            Add-Content -Path $supLog -Value ($ts + ' GIVING UP after ' +
                $restarts + ' relaunches, rows ' + $rows + '/' + $target)
            break
        }
        $restarts = $restarts + 1
        Add-Content -Path $supLog -Value ($ts + ' runner absent at ' + $rows +
            '/' + $target + ' -> relaunch #' + $restarts)
        & $launch | Add-Content -Path $supLog
    }
    Start-Sleep -Seconds 900
}
