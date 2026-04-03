#!/usr/bin/env bash
# ============================================================
# DESCARGA DATOS GEO - Cáncer de mama en mujeres mexicanas
# Romero-Cordoba et al., Nature Communications 2021
# DOI: 10.1038/s41467-021-22478-5
# SuperSeries: GSE87049 (contiene GSE87048 + GSE87050)
# ============================================================
# Ejecutar desde tu máquina local (requiere ~2-5 GB de espacio)
# Dependencias: wget o curl, python3 + GEOparse/pandas/Biopython

set -euo pipefail

OUTDIR="$HOME/geo_cancer_mama_mx"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== Descargando datos GEO: GSE87049 (SuperSeries) ==="
echo "Directorio de salida: $OUTDIR"
echo ""

# --------------------------------------------------------
# OPCIÓN A: Descarga directa via wget (más rápida, archivos raw)
# --------------------------------------------------------
download_series_matrix() {
    local ACC="$1"
    local PREFIX="${ACC:0:6}nnn"   # GSE870 -> GSE87nnn

    echo ">> Descargando series matrix para $ACC..."

    BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACC}/matrix"

    # Listar archivos disponibles
    wget -q -r -l1 --no-parent -A "*.txt.gz" \
        --no-check-certificate \
        "$BASE/" \
        -P "$OUTDIR/${ACC}_matrix/" 2>/dev/null || true

    # Alternativa directa si conocemos el nombre del archivo
    wget -q -nc \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACC}/matrix/${ACC}_series_matrix.txt.gz" \
        -O "${ACC}_series_matrix.txt.gz" \
        --show-progress 2>&1 || echo "   (archivo único no encontrado, puede haber múltiples matrices)"
}

download_soft() {
    local ACC="$1"
    local PREFIX="${ACC:0:6}nnn"

    echo ">> Descargando SOFT (metadatos completos) para $ACC..."
    wget -q -nc \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACC}/soft/${ACC}_family.soft.gz" \
        -O "${ACC}_family.soft.gz" \
        --show-progress 2>&1
}

download_supplementary() {
    local ACC="$1"
    local PREFIX="${ACC:0:6}nnn"

    echo ">> Descargando archivos suplementarios para $ACC..."
    # Listar y descargar todo lo que haya en suppl/
    wget -q -r -l1 --no-parent -nd \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${ACC}/suppl/" \
        -P "${ACC}_suppl/" \
        --show-progress 2>&1 || echo "   (sin suplementarios o error de red)"
}

# --- SubSeries de expresión: GSE87050 (microarray, 109 tumores + normales) ---
echo ""
echo "=== [1/3] SubSeries GSE87050: Microarray de expresión (Affymetrix HuGene 1.0 ST) ==="
echo "     109 tumores de mama + tejido normal adyacente (mujeres mexicanas)"
download_series_matrix "GSE87050"
download_soft "GSE87050"
download_supplementary "GSE87050"

# --- SubSeries de copy-number: GSE87048 (SNP6 array, 78 pares tumor/sangre) ---
echo ""
echo "=== [2/3] SubSeries GSE87048: SNP6 arrays / copy-number (Affymetrix SNP 6.0) ==="
echo "     78 tumores + sangre pareada (variación genómica somática)"
download_series_matrix "GSE87048"
download_soft "GSE87048"
download_supplementary "GSE87048"

# --- SuperSeries completa: GSE87049 ---
echo ""
echo "=== [3/3] SuperSeries GSE87049: SOFT completo (todos los metadatos) ==="
download_soft "GSE87049"

echo ""
echo "=== Descarga completada ==="
ls -lh "$OUTDIR"/*.gz 2>/dev/null || true
echo ""
echo "Siguiente paso: ejecutar el script de análisis en R o Python (ver archivos adjuntos)"
