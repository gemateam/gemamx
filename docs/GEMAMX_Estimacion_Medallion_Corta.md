# GEMAMX: estimacion corta de costos

## Alcance

Esta version corrige el alcance para estimar **solo** la plataforma de datos en **GCP + Databricks** bajo un esquema **medallion**:

- `Bronze`: ingesta de datos crudos
- `Silver`: limpieza, estandarizacion y enriquecimiento
- `Gold`: datasets curados para consulta y exportacion
- Databricks Jobs para ETL/ELT
- Databricks SQL Warehouse para analitica y dashboards
- Cloud Storage para persistencia de capas
- Monitoreo/logging y respaldo operativo minimo

## Lo que NO incluye

Queda fuera de esta estimacion:

- AlphaFold, Boltz, docking, GPUs cientificas
- OpenShift, microservicios, DevOps extendido
- Vertex AI, Gemini, RAG, MLOps avanzado
- soporte funcional, mesa de ayuda, consultoria
- capacitacion, licencias externas, Antigravity
- costos de desarrollo de software

## Supuestos base

- Horizonte: 3 anos
- Casos procesados por ano:
  - Ano 1: 9,000
  - Ano 2: 19,500
  - Ano 3: 30,000
- Tamano promedio persistido en medallion por caso:
  - Bronze: 2.0 GB
  - Silver: 0.8 GB
  - Gold: 0.2 GB
  - Total: 3.0 GB por caso
- Retencion: acumulativa durante los 3 anos
- Almacenamiento: `Cloud Storage Standard`, una sola region
- Databricks: uso orientado a `Jobs Compute` para pipelines y `SQL Warehouse` para consulta

## Supuestos de precios usados

- GCS Standard single-region: **$0.02 USD/GB-mes** como aproximacion de trabajo
- Operaciones/egress de storage: **20%** adicional sobre storage
- Databricks Jobs Compute: **$0.30 USD/DBU-h**
- Databricks SQL Warehouse Pro: **$0.55 USD/DBU-h**

Nota: los precios reales pueden variar por region, contrato, descuentos, tipo de workspace y consumo comprometido. Esta hoja sirve como **estimacion ejecutiva corta**, no como cotizacion contractual.

## Escenario recomendado

### Resumen anual

| Concepto | Ano 1 | Ano 2 | Ano 3 |
|---|---:|---:|---:|
| GCS storage | 3,240 | 13,500 | 31,320 |
| GCS operaciones/transfer | 648 | 2,700 | 6,264 |
| Databricks jobs | 16,200 | 30,600 | 46,800 |
| GCP compute para clusters | 10,000 | 20,000 | 32,000 |
| Databricks SQL warehouse | 6,970 | 10,454 | 13,939 |
| Monitoreo/logging/backups | 4,000 | 5,000 | 7,000 |
| **Total anual** | **41,058** | **82,254** | **137,323** |

### Total 3 anos

- **Total 3 anos estimado: $260,635 USD**

## Lectura ejecutiva

Si el objetivo es **solo** operar un lakehouse medallion en GCP con Databricks para ingesta, transformacion y consulta, el orden de magnitud ya no es de millones de dolares. Bajo estos supuestos, el proyecto cae mas cerca de:

- **Ano 1:** ~`$41K`
- **Ano 2:** ~`$82K`
- **Ano 3:** ~`$137K`
- **Total 3 anos:** ~`$261K`

## Recomendacion de uso

Usa esta version corta para:

- presentacion ejecutiva inicial
- comparativo rapido contra la version inflada
- pedir cotizacion formal a Google/Databricks
- base de un presupuesto de plataforma de datos, separado del presupuesto cientifico

La recomendacion es manejar **dos presupuestos separados**:

1. **Plataforma de datos medallion**
2. **Plataforma cientifica/IA/bioinformatica**

Mezclarlos en una sola hoja distorsiona la decision.

## Fuentes usadas para acotar el modelo

- Databricks pricing overview: https://www.databricks.com/product/pricing
- Serverless DBU consumption notes: https://docs.databricks.com/gcp/en/resources/pricing
- Google Cloud Storage pricing: https://cloud.google.com/storage/pricing
- BigQuery pricing referencial para warehouse alterno: https://cloud.google.com/bigquery/pricing.html
