param(
    [string]$InputRoot = "data",
    [string]$OutputRoot = "handoff\rebuilt_txt_parser"
)

$ErrorActionPreference = "Stop"

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path

function Resolve-RepoPath {
    param([string]$PathText)

    if ([System.IO.Path]::IsPathRooted($PathText)) {
        return (Resolve-Path $PathText).Path
    }

    return (Resolve-Path (Join-Path $repoRoot $PathText)).Path
}

$inputPath = Resolve-RepoPath $InputRoot

if ([System.IO.Path]::IsPathRooted($OutputRoot)) {
    $outputPath = $OutputRoot
} else {
    $outputPath = Join-Path $repoRoot $OutputRoot
}

New-Item -ItemType Directory -Force -Path $outputPath | Out-Null

$requiredPrefix = @("t_sec", "gps_week", "tow_sec", "constellation", "prn")
$placeholderColumns = @("gt_pos_x", "gt_pos_y", "gt_pos_z")
$manifestRows = New-Object System.Collections.Generic.List[string]
$manifestRows.Add("source_file`toutput_file`tadded_columns")

$candidateFiles = Get-ChildItem -Path $inputPath -Recurse -File -Filter *.csv | Sort-Object FullName
$processedCount = 0

foreach ($file in $candidateFiles) {
    $lines = [System.IO.File]::ReadAllLines($file.FullName)
    if ($lines.Length -eq 0) {
        continue
    }

    $headerColumns = $lines[0] -split "`t", -1
    if ($headerColumns.Length -lt $requiredPrefix.Length) {
        continue
    }

    $matchesPrefix = $true
    for ($i = 0; $i -lt $requiredPrefix.Length; $i++) {
        if ($headerColumns[$i] -ne $requiredPrefix[$i]) {
            $matchesPrefix = $false
            break
        }
    }

    if (-not $matchesPrefix) {
        continue
    }

    $missingColumns = @($placeholderColumns | Where-Object { $headerColumns -notcontains $_ })

    $newLines = New-Object string[] $lines.Length
    $newHeader = $lines[0]
    if ($missingColumns.Count -gt 0) {
        $newHeader += "`t" + ($missingColumns -join "`t")
    }
    $newLines[0] = $newHeader

    $suffix = ""
    if ($missingColumns.Count -gt 0) {
        $suffix = "`t" + ((@("nan") * $missingColumns.Count) -join "`t")
    }

    for ($lineIndex = 1; $lineIndex -lt $lines.Length; $lineIndex++) {
        if ([string]::IsNullOrWhiteSpace($lines[$lineIndex])) {
            $newLines[$lineIndex] = $lines[$lineIndex]
            continue
        }

        $newLines[$lineIndex] = $lines[$lineIndex] + $suffix
    }

    $destinationFile = Join-Path $outputPath $file.Name
    [System.IO.File]::WriteAllLines($destinationFile, $newLines)

    $addedColumnsText = if ($missingColumns.Count -gt 0) { $missingColumns -join "," } else { "-" }
    $manifestRows.Add("$($file.Name)`t$([System.IO.Path]::GetFileName($destinationFile))`t$addedColumnsText")
    $processedCount++
}

$manifestPath = Join-Path $outputPath "manifest.tsv"
[System.IO.File]::WriteAllLines($manifestPath, $manifestRows)

Write-Host ("Prepared {0} handoff file(s) in {1}" -f $processedCount, $outputPath)
