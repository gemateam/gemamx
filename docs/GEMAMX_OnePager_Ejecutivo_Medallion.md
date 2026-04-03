# GEMAMX
## Estimacion Ejecutiva de Plataforma de Datos

### Objetivo

Implementar en **Google Cloud Platform + Databricks** una plataforma de datos tipo **Lakehouse Medallion** para GEMAMX, enfocada exclusivamente en:

- ingesta de datos
- almacenamiento por capas `Bronze`, `Silver` y `Gold`
- transformacion y curacion de datos
- consulta analitica y exportacion de resultados
- monitoreo y operacion minima de la plataforma

### Alcance considerado

Esta estimacion contempla solo la **plataforma de datos**.

Incluye:

- `Cloud Storage` para capas medallion
- `Databricks Jobs` para pipelines de ingesta y transformacion
- `Databricks SQL Warehouse` para consulta y dashboards
- capacidad base de computo en GCP para ejecucion de cargas
- monitoreo, logging y respaldo operativo minimo

No incluye:

- prediccion estructural de proteinas
- docking molecular
- GPUs cientificas
- OpenShift o microservicios
- IA generativa, RAG o MLOps avanzado
- soporte funcional, capacitacion o licencias externas

### Supuestos de capacidad

- horizonte de planeacion: **3 anos**
- adopcion gradual:
  - **Ano 1:** 9,000 casos
  - **Ano 2:** 19,500 casos
  - **Ano 3:** 30,000 casos
- almacenamiento estimado por caso en medallion:
  - **3 GB por caso** acumulando `Bronze + Silver + Gold`

### Inversion estimada

| Concepto | Ano 1 | Ano 2 | Ano 3 | Total 3 anos |
|---|---:|---:|---:|---:|
| Storage GCP | 3,888 | 16,200 | 37,584 | 57,672 |
| Databricks Jobs + compute | 26,200 | 50,600 | 78,800 | 155,600 |
| SQL Warehouse | 6,970 | 10,454 | 13,939 | 31,363 |
| Monitoreo, logging y backups | 4,000 | 5,000 | 7,000 | 16,000 |
| **Total estimado** | **41,058** | **82,254** | **137,323** | **260,635** |

### Lectura ejecutiva

La plataforma de datos requerida para un esquema `medallion` en `GCP + Databricks` tiene un costo estimado de:

- **USD 41 mil** en el primer ano
- **USD 82 mil** en el segundo ano
- **USD 137 mil** en el tercer ano
- **USD 261 mil** acumulados a tres anos

Esto ubica el costo de la **plataforma de datos** en un rango manejable para una estrategia de crecimiento gradual, y claramente separado de los costos cientificos o de investigacion avanzada.

### Recomendacion

Para fines presupuestales y de gobernanza, se recomienda manejar **dos bolsas separadas**:

1. **Plataforma de datos**: GCP + Databricks medallion
2. **Plataforma cientifica**: bioinformatica, IA, docking, GPUs y servicios especializados

La estimacion presentada en esta hoja corresponde unicamente a la **primera bolsa**.

### Mensaje clave para presentacion

**El costo de montar la plataforma de datos medallion para GEMAMX no es de millones de dolares; bajo un despliegue gradual de tres anos, se estima en aproximadamente USD 261 mil.**
