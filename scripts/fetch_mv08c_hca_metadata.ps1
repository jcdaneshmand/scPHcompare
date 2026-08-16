param(
    [Parameter(Mandatory = $false)]
    [string]$OutputDirectory = "tmp/mv08c-hca-metadata"
)

$ErrorActionPreference = "Stop"
$ProjectId = "cc95ff89-2e68-4a08-a234-480eca21ce79"
$Catalog = "dcp60"
$ApiRoot = "https://service.azul.data.humancellatlas.org"
$Utf8NoBom = [System.Text.UTF8Encoding]::new($false)

function Resolve-RepositoryPath([string]$RelativePath) {
    $RepositoryRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
    return Join-Path $RepositoryRoot $RelativePath
}

function Save-JsonResponse(
    [string]$Name,
    [string]$Endpoint,
    [string]$FiltersJson
) {
    $Uri = $Endpoint + "?catalog=$Catalog&size=75&filters=" +
        [uri]::EscapeDataString($FiltersJson)
    $Response = Invoke-WebRequest -UseBasicParsing -Uri $Uri -Headers @{
        Accept = "application/json"
    }
    [IO.File]::WriteAllText(
        (Join-Path $ResolvedOutput $Name),
        $Response.Content,
        $Utf8NoBom
    )
}

$ContractPath = Resolve-RepositoryPath "docs/specifications/mv08c-hca-query-contract-v1.csv"
$Contract = Import-Csv -LiteralPath $ContractPath
if ($Contract.Count -ne 4 -or
    @($Contract.expression_download_authorized | Where-Object { $_ -ne "False" }).Count -ne 0 -or
    @($Contract.outcome_access_authorized | Where-Object { $_ -ne "False" }).Count -ne 0) {
    throw "MV8-C query contract is missing or authorizes forbidden payloads."
}

$ResolvedOutput = if ([IO.Path]::IsPathRooted($OutputDirectory)) {
    $OutputDirectory
} else {
    Resolve-RepositoryPath $OutputDirectory
}
New-Item -ItemType Directory -Force -Path $ResolvedOutput | Out-Null
if (Get-ChildItem -LiteralPath $ResolvedOutput -Force -ErrorAction SilentlyContinue) {
    throw "MV8-C metadata output directory must be empty."
}

$OpenApi = Invoke-WebRequest -UseBasicParsing -Uri "$ApiRoot/openapi.json" -Headers @{
    Accept = "application/json"
}
[IO.File]::WriteAllText(
    (Join-Path $ResolvedOutput "openapi.json"),
    $OpenApi.Content,
    $Utf8NoBom
)

$ProjectFilters = '{"projectId":{"is":["' + $ProjectId + '"]}}'
$WholeMarrowFilters = '{"projectId":{"is":["' + $ProjectId + '"]},' +
    '"organ":{"is":["bone marrow"]},' +
    '"selectedCellType":{"is":["mononuclear cell of bone marrow"]}}'
$SelectedMarrowFilters = '{"projectId":{"is":["' + $ProjectId + '"]},' +
    '"organPart":{"is":["bone marrow"]},' +
    '"selectedCellType":{"is":["bone marrow hematopoietic cell"]}}'
$PrimaryManifestFilters = '{"projectId":{"is":["' + $ProjectId + '"]},' +
    '"organ":{"is":["bone marrow"]},' +
    '"selectedCellType":{"is":["mononuclear cell of bone marrow"]},' +
    '"donorCount":{"is":[1]},"fileFormat":{"is":["h5"]}}'

Save-JsonResponse "project.json" "$ApiRoot/index/projects" $ProjectFilters
Save-JsonResponse "whole-marrow-samples.json" "$ApiRoot/index/samples" $WholeMarrowFilters
Save-JsonResponse "selected-marrow-samples.json" "$ApiRoot/index/samples" $SelectedMarrowFilters

$ManifestUri = "$ApiRoot/fetch/manifest/files?catalog=$Catalog&format=compact&filters=" +
    [uri]::EscapeDataString($PrimaryManifestFilters)
$State = Invoke-RestMethod -Method Put -Uri $ManifestUri -Headers @{
    Accept = "application/json"
}
for ($Attempt = 0; $Attempt -lt 40 -and [int]$State.Status -ne 302; $Attempt++) {
    if ([int]$State.Status -ne 301) {
        throw "Unexpected HCA manifest state $($State.Status)."
    }
    $Delay = [Math]::Max(1, [Math]::Min(10, [int]$State.'Retry-After'))
    Start-Sleep -Seconds $Delay
    $State = Invoke-RestMethod -Method Get -Uri $State.Location -Headers @{
        Accept = "application/json"
    }
}
if ([int]$State.Status -ne 302) {
    throw "HCA manifest job did not finish within the bounded polling window."
}

# The 302 location is an expiring manifest URL. It is used once and is never
# printed or persisted. The compact TSV contains metadata, not expression data.
$ManifestPath = Join-Path $ResolvedOutput "primary-h5-manifest.tsv"
Invoke-WebRequest -UseBasicParsing -Uri $State.Location -OutFile $ManifestPath
$Manifest = Import-Csv -Delimiter "`t" -LiteralPath $ManifestPath
if ($Manifest.Count -ne 8 -or
    @($Manifest.file_uuid | Sort-Object -Unique).Count -ne 8 -or
    @($Manifest.'donor_organism.provenance.document_id' | Sort-Object -Unique).Count -ne 8 -or
    @($Manifest | Where-Object { $_.file_format -ne "h5" }).Count -ne 0) {
    throw "The compact manifest does not contain the frozen eight-file cohort."
}

$ExpectedJsonCounts = @{
    "project.json" = 1
    "whole-marrow-samples.json" = 24
    "selected-marrow-samples.json" = 63
}
foreach ($Name in $ExpectedJsonCounts.Keys) {
    $Parsed = Get-Content -LiteralPath (Join-Path $ResolvedOutput $Name) -Raw |
        ConvertFrom-Json
    if ([int]$Parsed.pagination.count -ne $ExpectedJsonCounts[$Name]) {
        throw "Unexpected entity count in $Name."
    }
}

$Summary = Get-ChildItem -LiteralPath $ResolvedOutput -File | Sort-Object Name |
    ForEach-Object {
        [pscustomobject]@{
            name = $_.Name
            bytes = $_.Length
            sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $_.FullName).Hash.ToLowerInvariant()
        }
    }
$Summary | Format-Table -AutoSize
Write-Output "MV8-C metadata acquisition complete: 8 manifest rows; no expression payload requested."
