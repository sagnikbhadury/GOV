param()

Write-Host "FDR ANALYSIS PIPELINE - PowerShell"
Write-Host ""

$basePath = "C:\Users\bhadury\University of Michigan Dropbox\Sagnik Bhadury\CODES\SBTurboCodes\H&E Lung TCGA Graph Dataset"
$codesDir = Join-Path $basePath "codes"
$resultsDir = Join-Path $basePath "CLAUDE_REVISIONS_JULY2026\analysis_results"
$tablesDir = Join-Path $resultsDir "supplementary_tables"

New-Item -ItemType Directory -Path $tablesDir -Force -EA SilentlyContinue | Out-Null

Write-Host "Loading data..."

$pathways = Import-Csv (Join-Path $codesDir "luad_Activity_lsGPRVS_pip_longformat.csv")
$kinases = Import-Csv (Join-Path $codesDir "luad_kinase_lsGPRVS_pip_longformat.csv")
$tfs = Import-Csv (Join-Path $codesDir "luad_tf_lsGPRVS_pip_longformat.csv")

if ($pathways -is [pscustomobject]) { $pathways = @($pathways) }
if ($kinases -is [pscustomobject]) { $kinases = @($kinases) }
if ($tfs -is [pscustomobject]) { $tfs = @($tfs) }

Write-Host "Pathways: $($pathways.Count)"
Write-Host "Kinases: $($kinases.Count)"
Write-Host "TFs: $($tfs.Count)"
Write-Host ""

function ApplyBH {
    param([double[]] $p)
    $n = $p.Length
    if ($n -eq 0) { return @() }

    $idx = @()
    for ($i = 0; $i -lt $n; $i++) {
        $idx += @{i = $i; p = $p[$i]}
    }

    $srt = $idx | Sort-Object -Property p
    $q = New-Object double[] $n
    $cm = 1.0

    for ($r = $n; $r -ge 1; $r--) {
        $it = $srt[$n - $r]
        $qv = [Math]::Min($it.p * $n / $r, 1.0)
        if ($qv -lt $cm) { $cm = $qv }
        else { $cm = $qv }
        $q[$it.i] = $cm
    }
    return $q
}

function Process {
    param($d, $nm)

    $p = @()
    foreach ($r in $d) {
        $p += (1 - [double]$r.PIP)
    }

    $q = ApplyBH $p
    $res = @()
    $p55 = 0
    $f05 = 0

    for ($i = 0; $i -lt $d.Count; $i++) {
        $r = $d[$i]
        $pi = [double]$r.PIP
        $qv = $q[$i]

        if ($pi -gt 0.5) { $p55++ }
        if ($qv -lt 0.05) { $f05++ }

        $ft = $r.Pathways
        if (!$ft) { $ft = $r.Kinase }
        if (!$ft) { $ft = $r.TF }

        $o = [PSCustomObject]@{
            HPC = $r.HPC
            Feature = $ft
            PIP = [Math]::Round($pi, 4)
            p_value = [Math]::Round($p[$i], 6)
            q_value = [Math]::Round($qv, 6)
            risk = $r.risk
        }

        $res += $o
    }

    $res = $res | Sort-Object -Property q_value
    return @{res = $res; p55 = $p55; f05 = $f05}
}

Write-Host "Processing..."
$pw = Process $pathways "Pathways"
$ki = Process $kinases "Kinases"
$tf = Process $tfs "TFs"

Write-Host ""
Write-Host "RESULTS"
Write-Host ""
Write-Host "Pathways: PIP gt 0.5 = $($pw.p55), FDR lt 0.05 = $($pw.f05)"
Write-Host "Kinases: PIP gt 0.5 = $($ki.p55), FDR lt 0.05 = $($ki.f05)"
Write-Host "TFs: PIP gt 0.5 = $($tf.p55), FDR lt 0.05 = $($tf.f05)"
Write-Host ""

Write-Host "Saving tables..."
$pw.res | Export-Csv -Path (Join-Path $tablesDir "Table_S2_Pathway_Associations_FDR.csv") -NoTypeInformation
Write-Host "Table_S2_Pathway_Associations_FDR.csv"

$ki.res | Export-Csv -Path (Join-Path $tablesDir "Table_S3_Kinase_Associations_FDR.csv") -NoTypeInformation
Write-Host "Table_S3_Kinase_Associations_FDR.csv"

$tf.res | Export-Csv -Path (Join-Path $tablesDir "Table_S4_TF_Associations_FDR.csv") -NoTypeInformation
Write-Host "Table_S4_TF_Associations_FDR.csv"

$comp = @(
    [PSCustomObject]@{Layer = "Pathways"; PIP_gt_05 = $pw.p55; FDR_lt_005 = $pw.f05; Reduction_Pct = if ($pw.p55 -gt 0) { [Math]::Round(($pw.p55 - $pw.f05)/$pw.p55*100,1) } else { 0 } }
    [PSCustomObject]@{Layer = "Kinases"; PIP_gt_05 = $ki.p55; FDR_lt_005 = $ki.f05; Reduction_Pct = if ($ki.p55 -gt 0) { [Math]::Round(($ki.p55 - $ki.f05)/$ki.p55*100,1) } else { 0 } }
    [PSCustomObject]@{Layer = "TFs"; PIP_gt_05 = $tf.p55; FDR_lt_005 = $tf.f05; Reduction_Pct = if ($tf.p55 -gt 0) { [Math]::Round(($tf.p55 - $tf.f05)/$tf.p55*100,1) } else { 0 } }
)

$comp | Export-Csv -Path (Join-Path $tablesDir "Threshold_Comparison.csv") -NoTypeInformation
Write-Host "Threshold_Comparison.csv"
Write-Host ""

Write-Host "THRESHOLD COMPARISON"
Write-Host ""
$comp | Format-Table -AutoSize
Write-Host ""

Write-Host "TOP HITS"
Write-Host ""
Write-Host "Pathways (FDR lt 0.05):"
$pw.res | Where-Object {[double]$_.q_value -lt 0.05} | Select-Object -First 10 | Format-Table -AutoSize

Write-Host "Kinases (FDR lt 0.05):"
$ki.res | Where-Object {[double]$_.q_value -lt 0.05} | Select-Object -First 10 | Format-Table -AutoSize

Write-Host "TFs (FDR lt 0.05):"
$tf.res | Where-Object {[double]$_.q_value -lt 0.05} | Select-Object -First 10 | Format-Table -AutoSize

Write-Host ""
Write-Host "COMPLETE"
Write-Host ""
Write-Host "Results saved to: $tablesDir"
