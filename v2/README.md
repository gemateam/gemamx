# GEMAMX v2

Esta carpeta contiene la nueva implementacion del pipeline basada en el cohorte Romero.

## Objetivo

Construir una version nueva del pipeline sin afectar la version actual (`v1`) basada en TCGA/cBioPortal.

## Principios

- No modificar archivos de salida de `v1`
- No reutilizar `data/` de la raiz como espacio de trabajo de `v2`
- Mantener scripts, datos y documentos de `v2` aislados

## Estructura

- `data/`: salidas y artefactos de `v2`
- `docs/`: planes, notas y validaciones de la nueva implementacion
- `scripts/`: scripts nuevos de ingestion, transformacion y pipeline

## Fuente principal de mutaciones somaticas

La nueva implementacion toma como base:

- `data/romero/41467_2021_22478_MOESM8_ESM.xlsx`
- `data/romero/41467_2021_22478_MOESM9_ESM.xlsx`

En particular:

- hoja `HM maf`
- hoja `Mutsig HM`
- hoja `Driver mutation in HM`

## Estado actual

- carpeta base creada
- pipeline `v1` preservado
- ingestion Romero implementada
- FASTA WT/MUT generados para `v2`
- flujo GPU en Colab preparado para estructuras con Boltz

## Ejecucion con GPU en Colab

Para la fase pesada de estructuras en `v2`, usa:

- `v2/01_romero_structures_colab_gpu.ipynb`
- `v2/scripts/predict_structures_romero_boltz.py`
- `v2/scripts/scan_boltz_v2.py`

Flujo recomendado:
1. Generar localmente `v2/data/top_variants.json` y `v2/data/fasta_sequences.json`.
2. Abrir la libreta de Colab con runtime GPU.
3. Correr Boltz en Colab para WT/MUT.
4. Actualizar manifest y metricas.
5. Regresar a local para docking.
