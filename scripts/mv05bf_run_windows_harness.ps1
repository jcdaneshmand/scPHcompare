param(
  [Parameter(Mandatory = $true)][string]$Library,
  [Parameter(Mandatory = $true)][string]$EvidenceDirectory
)

$ErrorActionPreference = "Stop"
$vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path -LiteralPath $vswhere)) { throw "vswhere.exe is unavailable" }
$installation = & $vswhere -latest -products * `
  -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
  -property installationPath
if (-not $installation) { throw "Visual Studio C tools are unavailable" }
$developerCommand = Join-Path $installation "Common7\Tools\VsDevCmd.bat"
$repository = (Resolve-Path ".").Path
$libraryPath = (Resolve-Path $Library).Path
$evidence = (Resolve-Path $EvidenceDirectory).Path
$source = Join-Path $repository "scripts\mv05bf_windows_ffi_harness.c"
$include = Join-Path $repository "rust\scph_landscape_kernel\include"
$executable = Join-Path $evidence "ffi-harness.exe"
$dependencies = Join-Path $evidence "native-dependencies-raw.txt"
$command = "call `"$developerCommand`" -arch=x64 -host_arch=x64" +
  " && cl /nologo /W4 /WX /I`"$include`" `"$source`" /Fe:`"$executable`"" +
  " && `"$executable`" `"$libraryPath`"" +
  " && dumpbin /dependents `"$libraryPath`" > `"$dependencies`""
cmd.exe /d /s /c $command
if ($LASTEXITCODE -ne 0) { throw "Windows FFI harness or dependency inventory failed" }
