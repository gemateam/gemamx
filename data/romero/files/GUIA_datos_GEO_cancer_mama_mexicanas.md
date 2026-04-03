# Guía completa: Datos GEO - Cáncer de mama en mujeres mexicanas
## Romero-Cordoba et al., *Nature Communications* 2021
**DOI:** [10.1038/s41467-021-22478-5](https://doi.org/10.1038/s41467-021-22478-5)

---

## ✅ Accession Numbers confirmados (verificados en NCBI GEO)

| Accession | Tipo de dato | Plataforma | n muestras | Descripción |
|-----------|-------------|------------|------------|-------------|
| **GSE87049** | SuperSeries | — | 359 total | Contenedor de todas las series del estudio |
| **GSE87050** | Microarray expresión | Affymetrix HuGene 1.0 ST (GPL6244) | 109 tumores + normales adyacentes | Perfiles de expresión génica |
| **GSE87048** | Copy-number (SNP array) | Affymetrix SNP 6.0 (GPL6801) | 78 tumores + 100 sangres (buffy coat) | Alteraciones somáticas de número de copia |

> ⚠️ **Importante:** El estudio usa **microarray de expresión** (HuGene 1.0 ST), NO RNA-seq.
> Los datos de **WES** (n=134) no están en GEO; ver sección "Datos de exoma" abajo.

---

## 📥 Descarga directa (desde tu máquina local)

### Opción 1: wget por línea de comandos
```bash
mkdir -p ~/geo_cancer_mama_mx && cd ~/geo_cancer_mama_mx

# Series matrix de expresión (GSE87050)
PREFIX_50="GSE87nnn/GSE87050"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX_50}/matrix/GSE87050_series_matrix.txt.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX_50}/soft/GSE87050_family.soft.gz"

# Series matrix de copy-number (GSE87048)
PREFIX_48="GSE87nnn/GSE87048"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX_48}/matrix/GSE87048_series_matrix.txt.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX_48}/soft/GSE87048_family.soft.gz"

# Archivos suplementarios del SuperSeries
wget -r -l1 --no-parent \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87049/suppl/"
```

### Opción 2: GEOquery en R (recomendado)
```r
BiocManager::install("GEOquery")
library(GEOquery)

# Expresión génica (109 tumores)
gse_expr <- getGEO("GSE87050", GSEMatrix = TRUE)
expr_mat  <- exprs(gse_expr[[1]])
metadata  <- pData(gse_expr[[1]])

# Copy-number
gse_cn <- getGEO("GSE87048", GSEMatrix = TRUE)
```

### Opción 3: GEOparse en Python
```python
import GEOparse
gse = GEOparse.get_GEO("GSE87050", destdir="./geo_data/")
```

---

## 🧬 Estructura del dataset confirmada

### GSE87050 — Expresión (microarray)
- **109 muestras tumorales** nombradas como: `AT*`, `BaT*`, `BT*`, `HT*`, `NT*`, `T175`
  - `AT` = adenocarcinoma tubular / luminal A-like
  - `BaT` = basal-like / triple negativo
  - `BT` = luminal B / HER2+
  - `HT` = HER2-enriquecido
  - `NT` = sin receptor de estrógeno (ER-)
- **Tejido normal adyacente**: muestras `C*` (controles)

### GSE87048 — Copy-number (SNP6)
- **78 tumores** de mama (`T*`)
- **100 muestras de sangre** (buffy coat, `C*`) pareadas como control germinal

---

## 🔬 Datos de WES (Whole Exome Sequencing, n=134)

Los datos de exoma **NO están en GEO** debido a restricciones de privacidad de datos genómicos germinales.

**Dónde buscar:**

### 1. dbGaP (NCBI)
- URL: https://www.ncbi.nlm.nih.gov/gap/
- Buscar: `"Romero Cordoba" OR "Mexican breast cancer" OR PMID:33854067`
- Requiere: cuenta dbGaP + solicitud de acceso (Data Access Request)

### 2. EGA (European Genome-phenome Archive)
- URL: https://ega-archive.org/
- Buscar: DOI `10.1038/s41467-021-22478-5` o autor `Romero-Cordoba`
- Requiere: Data Access Committee (DAC) approval del INMEGEN (México)

### 3. Contacto directo con los autores
```
Sandra L. Romero-Córdoba
Instituto Nacional de Medicina Genómica (INMEGEN)
Email: sromero_cordoba@hotmail.com
```
Los datos de mutaciones somáticas (VCF/MAF) del WES pueden solicitarse directamente.

### 4. Suplementarios del paper
- URL: https://www.nature.com/articles/s41467-021-22478-5
- Buscar: "Supplementary Data" — pueden incluir tablas de mutaciones curadas

---

## 🖥️ Pipeline de análisis recomendado

### En R (script completo adjunto: `analisis_R_cancer_mama_mexicanas.R`)
```r
# 1. Descargar y normalizar
gse <- getGEO("GSE87050")

# 2. Control de calidad (boxplots, PCA)
plotMA(exprs(gse[[1]]))

# 3. Análisis diferencial (limma)
library(limma)
fit <- lmFit(exprs(gse[[1]]), design)
fit2 <- eBayes(contrasts.fit(fit, contrasts))

# 4. Anotación funcional (GO, KEGG)
library(clusterProfiler)
enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db)

# 5. Subtipificación molecular (PAM50)
library(genefu)
pam50_subtypes <- intrinsic.cluster.predict(sbt.model = pam50)
```

### En Python (script completo adjunto: `analisis_python_cancer_mama_mexicanas.py`)
```python
import GEOparse, pandas as pd, numpy as np
from scipy import stats

gse = GEOparse.get_GEO("GSE87050")
# → construir matriz, DEG, volcano plot
```

---

## 📦 Archivos adjuntos en este paquete

| Archivo | Descripción |
|---------|-------------|
| `descarga_GEO_cancer_mama_mexicanas.sh` | Script bash para descarga directa por wget |
| `analisis_R_cancer_mama_mexicanas.R` | Pipeline completo en R (GEOquery + limma + anotación) |
| `analisis_python_cancer_mama_mexicanas.py` | Pipeline equivalente en Python (GEOparse + scipy) |

---

## 📚 Citas relevantes

- Romero-Cordoba et al. (2021). "Loss of function of miR-342-3p results in MCT1 over-expression and contributes to oncogenic metabolic reprogramming in triple negative breast cancer." *Nat Commun* 12, 2402. https://doi.org/10.1038/s41467-021-22478-5
- PubMed ID: 33854067
- GEO SuperSeries: GSE87049 (público desde Aug 31, 2018)
