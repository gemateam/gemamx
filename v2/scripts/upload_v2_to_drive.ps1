param(
    [Parameter(Mandatory = $true)]
    [string]$DriveRoot
)

$ErrorActionPreference = "Stop"

function Ensure-Directory {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        New-Item -ItemType Directory -Path $Path | Out-Null
    }
}

$projectRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$sourceRoot = Join-Path $projectRoot "v2"
$targetRoot = Join-Path $DriveRoot "v2"

if (-not (Test-Path -LiteralPath $sourceRoot)) {
    throw "No se encontro la carpeta fuente: $sourceRoot"
}

Ensure-Directory -Path $targetRoot
Ensure-Directory -Path (Join-Path $targetRoot "data")
Ensure-Directory -Path (Join-Path $targetRoot "scripts")

$items = @(
    @{
        Type = "File"
        Source = Join-Path $sourceRoot "01_romero_structures_colab_gpu.ipynb"
        Target = Join-Path $targetRoot "01_romero_structures_colab_gpu.ipynb"
    },
    @{
        Type = "File"
        Source = Join-Path $sourceRoot "README.md"
        Target = Join-Path $targetRoot "README.md"
    },
    @{
        Type = "File"
        Source = Join-Path $sourceRoot "data\top_variants.json"
        Target = Join-Path $targetRoot "data\top_variants.json"
    },
    @{
        Type = "File"
        Source = Join-Path $sourceRoot "data\fasta_sequences.json"
        Target = Join-Path $targetRoot "data\fasta_sequences.json"
    },
    @{
        Type = "File"
        Source = Join-Path $sourceRoot "scripts\predict_structures_romero_boltz.py"
        Target = Join-Path $targetRoot "scripts\predict_structures_romero_boltz.py"
    },
    @{
        Type = "File"
        Source = Join-Path $sourceRoot "scripts\scan_boltz_v2.py"
        Target = Join-Path $targetRoot "scripts\scan_boltz_v2.py"
    },
    @{
        Type = "File"
        Source = Join-Path $sourceRoot "scripts\calc_structural_metrics_romero.py"
        Target = Join-Path $targetRoot "scripts\calc_structural_metrics_romero.py"
    },
    @{
        Type = "Directory"
        Source = Join-Path $sourceRoot "data\fasta"
        Target = Join-Path $targetRoot "data\fasta"
    }
)

Write-Output "Subiendo archivos de v2 a Google Drive..."
Write-Output "Fuente:  $sourceRoot"
Write-Output "Destino: $targetRoot"

foreach ($item in $items) {
    if (-not (Test-Path -LiteralPath $item.Source)) {
        Write-Warning "No existe: $($item.Source)"
        continue
    }

    $targetParent = Split-Path -Parent $item.Target
    Ensure-Directory -Path $targetParent

    if ($item.Type -eq "Directory") {
        if (Test-Path -LiteralPath $item.Target) {
            Remove-Item -LiteralPath $item.Target -Recurse -Force
        }
        Copy-Item -LiteralPath $item.Source -Destination $item.Target -Recurse -Force
        Write-Output "Directorio copiado: $($item.Target)"
    }
    else {
        Copy-Item -LiteralPath $item.Source -Destination $item.Target -Force
        Write-Output "Archivo copiado: $($item.Target)"
    }
}

Write-Output ""
Write-Output "Carga terminada."
Write-Output "Si en Colab vas a montar MyDrive/GEMAMX, la carpeta esperada queda en:"
Write-Output $targetRoot
