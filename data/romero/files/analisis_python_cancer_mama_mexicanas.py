#!/usr/bin/env python3
"""
Descarga y análisis de datos GEO: Cáncer de mama en mujeres mexicanas
Romero-Cordoba et al., Nat Commun 2021 | DOI: 10.1038/s41467-021-22478-5

Accession Numbers:
  GSE87049 = SuperSeries (todo el dataset)
  GSE87050 = Microarray de expresión HuGene 1.0 ST (109 tumores + normales)
  GSE87048 = Copy-number SNP6 arrays (78 pares tumor/sangre)

Uso:
  pip install GEOparse pandas numpy scipy matplotlib seaborn
  python3 analisis_python_cancer_mama_mexicanas.py
"""

import os
import sys
import GEOparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ─── Configuración ────────────────────────────────────────────
OUTDIR = Path.home() / "geo_cancer_mama_mx"
OUTDIR.mkdir(exist_ok=True)
os.chdir(OUTDIR)

print("=" * 60)
print("ANÁLISIS GEO: Cáncer de mama - Mujeres Mexicanas")
print("GSE87049 / GSE87050 / GSE87048")
print("=" * 60)

# ═══════════════════════════════════════════════════════════════
# BLOQUE 1: DESCARGA DE DATOS DE EXPRESIÓN (GSE87050)
# ═══════════════════════════════════════════════════════════════
print("\n[1/4] Descargando GSE87050 (microarray de expresión)...")
print("      Nota: primera descarga puede tardar 5-15 min (~50 MB)\n")

try:
    gse = GEOparse.get_GEO(
        geo="GSE87050",
        destdir=str(OUTDIR),
        silent=False
    )
    print(f"\n✓ GSE87050 cargado exitosamente")
    print(f"  Título: {gse.metadata.get('title', ['N/A'])[0]}")
    print(f"  Muestras (GSMs): {len(gse.gsms)}")
    print(f"  Plataformas (GPLs): {list(gse.gpls.keys())}")
except Exception as e:
    print(f"✗ Error al descargar GSE87050: {e}")
    sys.exit(1)

# ─── Extraer metadatos clínicos ────────────────────────────────
print("\n[2/4] Extrayendo metadatos clínicos...")

metadata_rows = []
for gsm_name, gsm in gse.gsms.items():
    row = {"gsm_id": gsm_name}
    # Extraer campos clave del SOFT
    for key in ["title", "source_name_ch1", "organism_ch1",
                "characteristics_ch1", "description"]:
        val = gsm.metadata.get(key, ["N/A"])
        row[key] = " | ".join(val) if isinstance(val, list) else str(val)
    metadata_rows.append(row)

metadata_df = pd.DataFrame(metadata_rows)
metadata_df.to_csv("GSE87050_metadata_clinico.csv", index=False)
print(f"✓ Metadatos guardados: GSE87050_metadata_clinico.csv")
print(f"  Columnas: {list(metadata_df.columns)}")

# Identificar tumores vs normales
tumor_mask  = metadata_df["source_name_ch1"].str.contains(
    "cancer|tumor|Breast cancer", case=False, na=False)
normal_mask = metadata_df["source_name_ch1"].str.contains(
    "normal|adjacent|Normal adjacent", case=False, na=False)

tumor_ids  = metadata_df.loc[tumor_mask, "gsm_id"].tolist()
normal_ids = metadata_df.loc[normal_mask, "gsm_id"].tolist()

print(f"\n  Muestras tumorales:    {len(tumor_ids)}")
print(f"  Tejido normal adyacente: {len(normal_ids)}")

# ─── Construir matriz de expresión ────────────────────────────
print("\n[3/4] Construyendo matriz de expresión...")

expr_data = {}
for gsm_name, gsm in gse.gsms.items():
    if gsm.table is not None and not gsm.table.empty:
        # La columna de expresión generalmente se llama 'VALUE'
        col_val = "VALUE" if "VALUE" in gsm.table.columns else gsm.table.columns[-1]
        expr_data[gsm_name] = gsm.table.set_index("ID_REF")[col_val]

if expr_data:
    expr_matrix = pd.DataFrame(expr_data)
    print(f"✓ Matriz de expresión: {expr_matrix.shape[0]} probes x {expr_matrix.shape[1]} muestras")

    # Separar tumores y normales
    tumor_cols  = [c for c in expr_matrix.columns if c in tumor_ids]
    normal_cols = [c for c in expr_matrix.columns if c in normal_ids]

    expr_tumor  = expr_matrix[tumor_cols]
    expr_normal = expr_matrix[normal_cols]

    print(f"  Tumores:  {expr_tumor.shape}")
    print(f"  Normales: {expr_normal.shape}")

    # Guardar (nota: archivos grandes ~500 MB en formato CSV)
    print("  Guardando matrices (puede tardar unos minutos)...")
    expr_tumor.to_csv("GSE87050_expr_tumores.csv.gz", compression="gzip")
    expr_normal.to_csv("GSE87050_expr_normales.csv.gz", compression="gzip")
    print("✓ Matrices guardadas (formato .csv.gz comprimido)")
else:
    print("✗ No se encontraron datos de expresión en las muestras.")

# ─── Análisis diferencial simple ──────────────────────────────
print("\n[4/4] Calculando expresión diferencial básica (t-test)...")

from scipy import stats

results = []
for probe in expr_matrix.index:
    t_vals = expr_matrix.loc[probe, tumor_cols].dropna().values
    n_vals = expr_matrix.loc[probe, normal_cols].dropna().values

    if len(t_vals) < 3 or len(n_vals) < 3:
        continue

    logfc  = np.mean(t_vals) - np.mean(n_vals)  # ya están en escala log
    t_stat, pval = stats.ttest_ind(t_vals, n_vals, equal_var=False)
    results.append({"probe_id": probe, "logFC": logfc, "p_value": pval})

results_df = pd.DataFrame(results)

# Corrección FDR (Benjamini-Hochberg)
from scipy.stats import rankdata
n = len(results_df)
ranked = rankdata(results_df["p_value"])
results_df["adj_p_BH"] = np.minimum(1, results_df["p_value"] * n / ranked)
results_df = results_df.sort_values("adj_p_BH")

sig_df = results_df[(results_df["adj_p_BH"] < 0.05) & (results_df["logFC"].abs() > 1)]
print(f"✓ Genes diferenciales (FDR<5%, |logFC|>1): {len(sig_df)}")
print(f"\nTop 10 sobreexpresados en tumores:")
print(sig_df[sig_df["logFC"] > 0].head(10)[["probe_id", "logFC", "adj_p_BH"]].to_string())
print(f"\nTop 10 subexpresados en tumores:")
print(sig_df[sig_df["logFC"] < 0].head(10)[["probe_id", "logFC", "adj_p_BH"]].to_string())

results_df.to_csv("GSE87050_DEG_completo.csv", index=False)
sig_df.to_csv("GSE87050_DEG_significativos.csv", index=False)

# ─── Volcano plot ──────────────────────────────────────────────
print("\nGenerando volcano plot...")

fig, ax = plt.subplots(figsize=(10, 8))
plot_df = results_df.copy()
plot_df["-log10_pval"] = -np.log10(plot_df["adj_p_BH"].clip(lower=1e-300))

# Colorear puntos
colors = np.where(
    (plot_df["adj_p_BH"] < 0.05) & (plot_df["logFC"] > 1), "red",
    np.where(
        (plot_df["adj_p_BH"] < 0.05) & (plot_df["logFC"] < -1), "blue",
        "grey"
    )
)

ax.scatter(plot_df["logFC"], plot_df["-log10_pval"],
           c=colors, alpha=0.4, s=5, rasterized=True)
ax.axhline(-np.log10(0.05), ls="--", c="black", lw=0.8, label="FDR=5%")
ax.axvline(-1, ls="--", c="navy", lw=0.8)
ax.axvline(1,  ls="--", c="darkred", lw=0.8)
ax.set_xlabel("log2 Fold Change (Tumor vs Normal)", fontsize=13)
ax.set_ylabel("-log10(FDR-adjusted p-value)", fontsize=13)
ax.set_title("Expresión diferencial: Cáncer de mama en mujeres mexicanas\n"
             "GSE87050 | Romero-Cordoba et al. 2021", fontsize=13)

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="red",  label=f"Sobreexpresado en tumor ({sum(colors=='red'):,})"),
    Patch(facecolor="blue", label=f"Subexpresado en tumor ({sum(colors=='blue'):,})"),
    Patch(facecolor="grey", label=f"No significativo ({sum(colors=='grey'):,})")
]
ax.legend(handles=legend_elements, loc="upper left")

plt.tight_layout()
plt.savefig("GSE87050_volcano_plot.png", dpi=150, bbox_inches="tight")
plt.close()
print("✓ Volcano plot guardado: GSE87050_volcano_plot.png")

# ─── Descargar GSE87048 (Copy-Number) ──────────────────────────
print("\n" + "="*60)
print("DESCARGANDO GSE87048 (copy-number SNP6 arrays)...")
print("="*60)

try:
    gse_cn = GEOparse.get_GEO(
        geo="GSE87048",
        destdir=str(OUTDIR),
        silent=False
    )
    meta_cn = []
    for gsm_name, gsm in gse_cn.gsms.items():
        row = {"gsm_id": gsm_name}
        for key in ["title", "source_name_ch1", "characteristics_ch1"]:
            val = gsm.metadata.get(key, ["N/A"])
            row[key] = " | ".join(val) if isinstance(val, list) else str(val)
        meta_cn.append(row)
    meta_cn_df = pd.DataFrame(meta_cn)
    meta_cn_df.to_csv("GSE87048_metadata_CNV.csv", index=False)
    print(f"✓ GSE87048 descargado: {len(meta_cn_df)} muestras")
    print(f"  Tipos: {meta_cn_df['source_name_ch1'].value_counts().to_dict()}")
except Exception as e:
    print(f"✗ Error GSE87048: {e}")

# ─── Resumen final ─────────────────────────────────────────────
print("\n" + "="*60)
print("RESUMEN: Archivos generados en", OUTDIR)
print("="*60)

output_files = [
    ("GSE87050_metadata_clinico.csv",       "Metadatos clínicos de las 109+ muestras"),
    ("GSE87050_expr_tumores.csv.gz",         "Matriz de expresión de tumores (~109 muestras)"),
    ("GSE87050_expr_normales.csv.gz",        "Matriz de expresión de tejido normal adyacente"),
    ("GSE87050_DEG_completo.csv",            "Resultados completos de expresión diferencial"),
    ("GSE87050_DEG_significativos.csv",      "Genes significativos (FDR<5%, |logFC|>1)"),
    ("GSE87050_volcano_plot.png",            "Volcano plot del análisis diferencial"),
    ("GSE87048_metadata_CNV.csv",            "Metadatos copy-number (78 pares tumor/sangre)"),
]

for fname, desc in output_files:
    fpath = OUTDIR / fname
    exists = "✓" if fpath.exists() else "✗"
    size   = f"{fpath.stat().st_size / 1e6:.1f} MB" if fpath.exists() else "no generado"
    print(f"  {exists} {fname:45s} {size}")
    print(f"    → {desc}")

print("\n⚠️  NOTA sobre WES (Whole Exome Sequencing, n=134):")
print("   Los datos de exoma NO están en GEO por restricciones de privacidad.")
print("   Buscar en:")
print("   • dbGaP: https://www.ncbi.nlm.nih.gov/gap/")
print("     Término: 'Romero-Cordoba breast cancer Mexico'")
print("   • EGA:   https://ega-archive.org/")
print("     Término: EGAD o EGAS asociado al paper DOI:10.1038/s41467-021-22478-5")
print("   • Suplementarios del paper en Nature Communications")
