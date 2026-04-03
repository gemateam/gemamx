# Plan de trabajo v2 Romero

## Objetivo

Recrear el pipeline GEMAMX usando el cohorte Romero como fuente de mutaciones somaticas, preservando el flujo general:

1. ingestión
2. priorización
3. FASTA WT/mutante
4. estructuras
5. docking
6. propuesta de farmacos

## Fuente validada

- `D:\ProyectosIA\GEMAMX\data\romero\41467_2021_22478_MOESM8_ESM.xlsx`

Hojas relevantes:

- `HM maf`
- `Mutsig HM`

## Fases

### Fase 1

Definir columnas y normalizacion desde `HM maf` hacia:

- `hugo_symbol`
- `protein_change`
- `chromosome`
- `start_position`
- `sample_id`
- `variant_classification`

### Fase 2

Generar en `v2/data/`:

- `raw_mutations.json`
- `top_variants.json`
- `mutation_stats.json`

### Fase 3

Adaptar scripts nuevos en `v2/scripts/` para:

- generar FASTA
- preparar estructuras
- preparar docking

### Fase 4

Decidir si la UI de `v2` sera:

- una copia aislada de la app actual
- o una variante configurable por dataset

## Regla de seguridad

Todo cambio nuevo debe vivir primero en `v2/`.
