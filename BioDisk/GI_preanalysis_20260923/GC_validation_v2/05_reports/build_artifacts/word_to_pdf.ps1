param([string]$DocPath, [string]$PdfPath)
$ErrorActionPreference = 'Stop'
$taskWord = $null
$taskDocument = $null
try {
    $taskWord = New-Object -ComObject Word.Application
    $taskWord.Visible = $false
    $taskWord.DisplayAlerts = 0
    $taskWord.AutomationSecurity = 3
    $taskDocument = $taskWord.Documents.Open($DocPath, $false, $true)
    $taskDocument.Repaginate()
    $taskDocument.ExportAsFixedFormat($PdfPath, 17)
    if (-not (Test-Path -LiteralPath $PdfPath)) { throw 'Word did not produce PDF' }
} finally {
    if ($null -ne $taskDocument) { $taskDocument.Close(0); [void][System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($taskDocument) }
    if ($null -ne $taskWord) { $taskWord.Quit(); [void][System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($taskWord) }
}
