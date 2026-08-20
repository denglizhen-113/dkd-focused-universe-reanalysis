param([Parameter(Mandatory=$true)][string]$DocxPath,[Parameter(Mandatory=$true)][string]$PdfPath)
$resolvedDocx=(Resolve-Path -LiteralPath $DocxPath).Path
$pdfParent=[System.IO.Path]::GetDirectoryName([System.IO.Path]::GetFullPath($PdfPath))
[System.IO.Directory]::CreateDirectory($pdfParent)|Out-Null
$word=$null;$document=$null
try {$word=New-Object -ComObject Word.Application;$word.Visible=$false;$word.DisplayAlerts=0;$document=$word.Documents.Open($resolvedDocx,$false,$true);$document.ExportAsFixedFormat($PdfPath,17)}
finally {if($null -ne $document){$document.Close($false);[System.Runtime.InteropServices.Marshal]::ReleaseComObject($document)|Out-Null};if($null -ne $word){$word.Quit();[System.Runtime.InteropServices.Marshal]::ReleaseComObject($word)|Out-Null};[GC]::Collect();[GC]::WaitForPendingFinalizers()}
if(-not(Test-Path -LiteralPath $PdfPath)){throw "Word did not create PDF: $PdfPath"}
