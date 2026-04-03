# ============================================================
# Análisis de datos GEO: Cáncer de mama en mujeres mexicanas
# Romero-Cordoba et al., Nat Commun 2021
# DOI: 10.1038/s41467-021-22478-5
# ============================================================
# GSE87049 = SuperSeries
# GSE87050 = Expresión microarray (109 tumores + normales)
# GSE87048 = Copy-number SNP6 (78 pares tumor/sangre)
# ============================================================

# --- Instalación de dependencias (ejecutar solo una vez) ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",        # descarga desde GEO
  "limma",           # análisis diferencial de expresión
  "oligo",           # normalización arrays Affymetrix (RMA)
  "hugene10sttranscriptcluster.db", # anotación HuGene 1.0 ST
  "affy",
  "Biobase"
), update = FALSE, ask = FALSE)

install.packages(c("tidyverse", "pheatmap", "RColorBrewer", "ggrepel"))

# ---- Librerías -----------------------------------------------
library(GEOquery)
library(Biobase)
library(tidyverse)
library(pheatmap)

# =============================================================
# BLOQUE 1: DESCARGA Y EXPLORACIÓN DE METADATOS
# =============================================================

# Descargar la SubSeries de EXPRESIÓN (GSE87050)
# Nota: La primera vez tarda ~5-10 min (archivo ~50 MB)
cat("Descargando GSE87050 (microarray de expresión)...\n")
gse_expr <- getGEO("GSE87050", GSEMatrix = TRUE, getGPL = FALSE)

# Si hay varias plataformas, seleccionar la correcta
if (length(gse_expr) > 1) {
  cat("Plataformas disponibles:", names(gse_expr), "\n")
  # Seleccionar GPL6244 = HuGene 1.0 ST
  gse <- gse_expr[["GSE87050-GPL6244_series_matrix.txt.gz"]]
} else {
  gse <- gse_expr[[1]]
}

cat("Dimensiones de la matriz de expresión:", dim(exprs(gse)), "\n")
cat("Número de muestras:", ncol(exprs(gse)), "\n")
cat("Número de probesets:", nrow(exprs(gse)), "\n")

# Extraer metadatos clínicos
pheno <- pData(gse)
cat("\nColumnas de metadatos disponibles:\n")
print(colnames(pheno))

# Ver las primeras filas de las columnas más relevantes
cat("\nResumen de tipos de muestra:\n")
if ("source_name_ch1" %in% colnames(pheno)) {
  print(table(pheno$source_name_ch1))
}

# Guardar metadatos en CSV para explorar
write.csv(pheno, "GSE87050_metadata_clinico.csv", row.names = TRUE)
cat("Metadatos guardados en: GSE87050_metadata_clinico.csv\n")

# =============================================================
# BLOQUE 2: SEPARAR TUMORES VS TEJIDO NORMAL ADYACENTE
# =============================================================

# Identificar tumores y normales por nombre de fuente
tumor_idx  <- grepl("cancer|tumor|Breast cancer", pheno$source_name_ch1,
                    ignore.case = TRUE)
normal_idx <- grepl("normal|adjacent|Normal adjacent", pheno$source_name_ch1,
                    ignore.case = TRUE)

cat("\nMuestras tumorales:   ", sum(tumor_idx), "\n")
cat("Tejido normal adyacente:", sum(normal_idx), "\n")

# Extraer matrices de expresión separadas
expr_all    <- exprs(gse)
expr_tumor  <- expr_all[, tumor_idx]
expr_normal <- expr_all[, normal_idx]

cat("Matriz tumores:", dim(expr_tumor), "\n")
cat("Matriz normales:", dim(expr_normal), "\n")

# Guardar matrices
write.csv(expr_tumor,  "GSE87050_expr_tumores.csv")
write.csv(expr_normal, "GSE87050_expr_normales.csv")

# =============================================================
# BLOQUE 3: ANOTACIÓN DE PROBESETS -> SÍMBOLOS DE GENES
# =============================================================
library(hugene10sttranscriptcluster.db)

# Mapear probe IDs a símbolos de genes
probe_ids <- rownames(expr_all)
gene_symbols <- mapIds(
  hugene10sttranscriptcluster.db,
  keys     = probe_ids,
  column   = "SYMBOL",
  keytype  = "PROBEID",
  multiVals = "first"
)
gene_entrez <- mapIds(
  hugene10sttranscriptcluster.db,
  keys     = probe_ids,
  column   = "ENTREZID",
  keytype  = "PROBEID",
  multiVals = "first"
)

annotation_df <- data.frame(
  probe_id     = probe_ids,
  gene_symbol  = gene_symbols,
  entrez_id    = gene_entrez,
  stringsAsFactors = FALSE
)

write.csv(annotation_df, "GSE87050_anotacion_probes.csv", row.names = FALSE)
cat("Anotación guardada: GSE87050_anotacion_probes.csv\n")
cat("Probes con símbolo génico:", sum(!is.na(gene_symbols)), "/", length(probe_ids), "\n")

# =============================================================
# BLOQUE 4: ANÁLISIS DE EXPRESIÓN DIFERENCIAL (limma)
# =============================================================
library(limma)

# Crear grupo factor
group <- factor(ifelse(tumor_idx, "Tumor", "Normal"))

# Diseño experimental
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Ajustar modelo lineal
fit <- lmFit(expr_all[, tumor_idx | normal_idx], design)
contrast_matrix <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extraer resultados
results <- topTable(fit2, number = Inf, adjust.method = "BH")
results$probe_id    <- rownames(results)
results$gene_symbol <- gene_symbols[rownames(results)]

# Filtrar genes significativos
sig_genes <- results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  arrange(desc(abs(logFC)))

cat("\nGenes significativamente diferenciales (FDR<5%, |logFC|>1):", nrow(sig_genes), "\n")
cat("Top 10 genes sobreexpresados en tumores:\n")
print(head(sig_genes %>% filter(logFC > 0) %>% select(gene_symbol, logFC, adj.P.Val), 10))

cat("Top 10 genes subexpresados en tumores:\n")
print(head(sig_genes %>% filter(logFC < 0) %>% select(gene_symbol, logFC, adj.P.Val), 10))

write.csv(results, "GSE87050_DEG_results_completo.csv", row.names = FALSE)
write.csv(sig_genes, "GSE87050_DEG_significativos.csv", row.names = FALSE)

# =============================================================
# BLOQUE 5: COPY-NUMBER (GSE87048 - SNP6 arrays)
# =============================================================
cat("\n\n=== Descargando GSE87048 (copy-number SNP6 arrays) ===\n")
gse_cnv <- getGEO("GSE87048", GSEMatrix = TRUE, getGPL = FALSE)
gse_cn  <- gse_cnv[[1]]

pheno_cn <- pData(gse_cn)
cat("Muestras en GSE87048:", nrow(pheno_cn), "\n")
cat("Tipos de muestra:\n")
print(table(pheno_cn$source_name_ch1))

write.csv(pData(gse_cn), "GSE87048_metadata_CNV.csv")

# =============================================================
# BLOQUE 6: HEATMAP DE LOS TOP GENES DIFERENCIALES
# =============================================================
library(pheatmap)
library(RColorBrewer)

# Seleccionar top 50 genes diferenciales
top50_probes <- head(sig_genes$probe_id, 50)
top50_symbols <- sig_genes$gene_symbol[sig_genes$probe_id %in% top50_probes]

# Matriz de expresión para heatmap
hm_matrix <- expr_all[top50_probes, tumor_idx | normal_idx]
rownames(hm_matrix) <- ifelse(is.na(top50_symbols), top50_probes, top50_symbols)

# Escalar por fila (z-score)
hm_scaled <- t(scale(t(hm_matrix)))

# Anotación de columnas
col_annotation <- data.frame(
  Tipo = ifelse(colSums(tumor_idx | normal_idx) > 0, group, NA),
  row.names = colnames(hm_matrix)
)
col_annotation$Tipo <- group[tumor_idx | normal_idx]

# Guardar heatmap
png("GSE87050_heatmap_top50_DEG.png", width = 1200, height = 1000, res = 150)
pheatmap(
  hm_scaled,
  annotation_col   = col_annotation,
  color            = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  show_colnames    = FALSE,
  fontsize_row     = 7,
  main             = "Top 50 DEGs: Tumores vs Normal - Mujeres Mexicanas (GSE87050)",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean"
)
dev.off()
cat("Heatmap guardado: GSE87050_heatmap_top50_DEG.png\n")

cat("\n=== ANÁLISIS COMPLETO ===\n")
cat("Archivos generados:\n")
cat("  - GSE87050_metadata_clinico.csv\n")
cat("  - GSE87050_expr_tumores.csv\n")
cat("  - GSE87050_expr_normales.csv\n")
cat("  - GSE87050_anotacion_probes.csv\n")
cat("  - GSE87050_DEG_results_completo.csv\n")
cat("  - GSE87050_DEG_significativos.csv\n")
cat("  - GSE87048_metadata_CNV.csv\n")
cat("  - GSE87050_heatmap_top50_DEG.png\n")
