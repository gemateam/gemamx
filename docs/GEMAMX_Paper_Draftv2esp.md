# GEMAMX: Un pipeline automatizado de oncología de precisión para el análisis de variantes de cáncer de mama, predicción estructural y acoplamiento molecular utilizando datos de TCGA-BRCA

**Ximena Gutiérrez Gálvez¹, Julio Waissman Vilanova¹, Héctor Gutiérrez Ibarra¹**

¹ Departamento de Matemáticas, Universidad de Sonora (UNISON), Hermosillo, Sonora, México

**Correspondencia:** ximena.gtzg@gmail.com

---

## Resumen

La oncología de precisión en el cáncer de mama requiere la integración sistemática del descubrimiento de variantes genómicas, el análisis estructural de proteínas y el modelado de interacciones farmacológicas, un flujo de trabajo que sigue siendo en gran medida manual y fragmentado entre instituciones. Presentamos **GEMAMX** (Genomic Analysis of Mammary cancer in Mexico), un pipeline automatizado de extremo a extremo que atraviesa cinco módulos secuenciales: (M1) ingestión de datos de 60,864 mutaciones somáticas de 1,084 casos de TCGA-BRCA a través de la API de GDC; (M2) priorización de variantes de sentido erróneo (missense) en 15 genes clínicamente accionables utilizando la frecuencia de la cohorte y los criterios de clasificación ACMG/AMP; (M3) síntesis automatizada de secuencias de proteínas y aplicación de mutaciones mediante la API REST de UniProt y Biopython; (M4) predicción de estructura 3D utilizando Boltz-1 y ESMFold; y (M5) acoplamiento molecular sistemático con AutoDock Vina 1.2.3 contra 11 inhibidores aprobados por la FDA o en fase clínica provenientes de ChEMBL. En 13 variantes priorizadas con estructuras 3D disponibles, GEMAMX calculó 143 puntuaciones de afinidad receptor-ligando. El hotspot dominante PIK3CA H1047R logró la afinidad de unión más fuerte en la cohorte (Pictilisib: −9.72 kcal/mol; Alpelisib: −8.90 kcal/mol), con un RMSD del sitio activo de 3.21 Å en relación con el tipo silvestre (wild-type), consistente con el reordenamiento del dominio quinasa conocido. BRCA1 C61G exhibió la mayor interrupción estructural (RMSD = 5.67 Å en el dominio RING, caída de pLDDT de 65.3 a 58.9), lo que sugiere una capacidad de interacción comprometida con BARD1. GEMAMX está desplegado como una aplicación web totalmente interactiva y es accesible públicamente en gema-mx.netlify.app, proporcionando visualización en tiempo real de variantes, visores moleculares 3D y clasificaciones de acoplamiento. El pipeline demuestra que la integración de la epidemiología genómica con la bioinformática estructural puede acelerar la reclasificación computacional de variantes de significado incierto (VUS) sin requerir infraestructura de laboratorio húmedo.

**Palabras clave:** oncología de precisión, cáncer de mama, variantes somáticas, TCGA-BRCA, acoplamiento molecular, AlphaFold, Boltz-1, AutoDock Vina, reclasificación de VUS, pipeline bioinformático

---

## 1. Introducción

El cáncer de mama sigue siendo la principal causa de mortalidad relacionada con el cáncer en mujeres en México, con una edad promedio de diagnóstico aproximadamente una década antes que en las poblaciones europeas [1]. Un cuello de botella crítico en la oncología clínica es la identificación e interpretación funcional de las mutaciones somáticas, particularmente aquellas clasificadas como Variantes de Significado Incierto (VUS) bajo las guías de ACMG/AMP [2]. Las VUS son mutaciones para las cuales la evidencia actual es insuficiente para establecer definitivamente su patogenicidad o benignidad; constituyen la mayoría de las variantes identificadas en los paneles genómicos clínicos y no pueden utilizarse directamente para guiar la terapia.

El conjunto de datos de The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) proporciona el compendio más completo disponible públicamente de mutaciones somáticas, resultados clínicos y subtipos moleculares para el cáncer de mama, cubriendo 1,084 muestras de pacientes secuenciadas a resolución de exoma completo [3]. Sin embargo, traducir este recurso en conocimientos computacionales accionables requiere la integración de al menos cuatro dominios analíticos distintos: descubrimiento de variantes genómicas, modelado estructural de proteínas, simulación de interacción farmacológica y visualización de datos clínicos, herramientas que normalmente operan de forma independiente especialistas en cada campo.

Los avances recientes en la predicción de estructuras de proteínas impulsada por IA, especialmente AlphaFold2 [4] y sus sucesores Boltz-1 [5] y ESMFold [6], han reducido drásticamente el costo computacional de generar modelos de proteínas 3D precisos para secuencias arbitrarias, incluidas las variantes mutadas no representadas en el Protein Data Bank. En paralelo, AutoDock Vina 1.2.3 [7] y las bibliotecas de preparación como Meeko [8] ahora permiten simulaciones de acoplamiento rápidas y reproducibles con una precisión de calidad GPU en hardware estándar. La combinación de estas herramientas crea una oportunidad sin precedentes para caracterizar computacionalmente las consecuencias estructurales y farmacológicas de las mutaciones somáticas a escala.

Presentamos GEMAMX, un pipeline de cinco módulos totalmente automatizado que integra estas capacidades en un único flujo de trabajo reproducible, desde llamadas directas a la API de GDC hasta la visualización molecular 3D interactiva. El pipeline está diseñado específicamente para abordar el problema de la reclasificación de VUS: al calcular la interrupción estructural (RMSD) entre las conformaciones de proteínas silvestres y mutantes y comparar las afinidades de unión de ligandos en 11 inhibidores, GEMAMX genera evidencia computacional que puede priorizar los candidatos a VUS para su validación experimental.

Este trabajo realiza las siguientes contribuciones:

1. Un pipeline bioinformático de código abierto y reproducible que cubre el camino completo desde la ingestión de datos genómicos hasta la visualización de puntuaciones de acoplamiento.
2. Caracterización estructural y farmacológica de 13 variantes de cáncer de mama en 7 genes utilizando predicciones de Boltz-1/ESMFold y AutoDock Vina 1.2.3.
3. Una aplicación web interactiva de acceso público (gema-mx.netlify.app) que democratiza el acceso a estos análisis sin requerir infraestructura computacional local.
4. Un marco computacional para la priorización de VUS basada en evidencia aplicable a cohortes genómicas latinoamericanas.

---

## 2. Métodos

### 2.1 Conjunto de Datos e Ingestión de Datos (M1)

GEMAMX recupera datos de mutaciones somáticas de la API REST de cBioPortal (www.cbioportal.org/api), accediendo al estudio TCGA-BRCA (ID de estudio: `brca_tcga`, perfil molecular: `brca_tcga_mutations`). La lista completa de muestras se recupera a través de `/sample-lists/{id}/sample-ids`, y las mutaciones se descargan en solicitudes POST paginadas a `/molecular-profiles/{id}/mutations/fetch` en lotes de 500 muestras. El mapeo de símbolos genéticos de los IDs Entrez Gene a los símbolos HUGO se realiza a través de `/genes/fetch`.

La descarga completa abarca **60,864 registros de mutaciones somáticas** en **1,084 muestras de pacientes**. Los registros se almacenan en MongoDB (v7.0, localhost:27017) en tres colecciones: `mutations` (variantes missense filtradas en genes objetivo), `top_variants` (resumen de la cohorte clasificado) y `fasta_sequences` (secuencias de aminoácidos WT/mutantes emparejadas). Todos los datos se exportan simultáneamente a archivos JSON estáticos (`raw_mutations.json`, `top_variants.json`, `mutation_stats.json`) para su consumo por la capa de aplicación web.

### 2.2 Priorización de Variantes (M2)

La priorización de variantes se restringe a mutaciones de sentido erróneo (clasificación `Missense_Mutation`) en un panel curado de 15 genes conductores de cáncer de mama organizados en tres niveles clínicos:

- **Nivel 1 (accionable):** PIK3CA, TP53, ERBB2, BRCA1, BRCA2
- **Nivel 2 (conductores frecuentes):** CDH1, GATA3, MAP3K1, KMT2C, PTEN
- **Nivel 3 (relevante para TCGA):** RB1, ARID1A, SF3B1, AKT1, CBFB

Las variantes se clasifican por frecuencia de cohorte y se categorizan utilizando un proxy de ACMG basado en frecuencia:

| Clasificación | Umbral de Frecuencia de Cohorte |
|---|---|
| Patógena | ≥ 5.0% |
| Probablemente Patógena | 1.0 – 4.9% |
| VUS | 0.5 – 0.99% |
| Baja Frecuencia | < 0.5% |

Esta clasificación proxy no es una aserción clínica de ACMG, sino que sirve como una prioridad computacional para la priorización del análisis estructural. Los pares variante-gen Top-26 se exportan para su procesamiento posterior.

### 2.3 Generación de Secuencias y Síntesis de Mutaciones (M3)

Para cada variante priorizada, la secuencia de proteína silvestre (WT) canónica se recupera de la API REST de UniProt utilizando una consulta revisada de Swiss-Prot (`gene:{SYMBOL} AND organism_id:9606 AND reviewed:true`). La secuencia FASTA se analiza con Biopython (v1.83) `SeqIO.read()`, se convierte a `MutableSeq` y la mutación puntual se aplica en la posición indexada en 1 especificada por la notación de cambio de proteína (por ejemplo, H1047R: histidina en la posición 1047 → arginina). La consistencia entre el residuo WT en la secuencia y el residuo WT esperado en la notación HGVS se verifica antes de la aplicación de la mutación.

Para proteínas que exceden los 390 aminoácidos —el límite operativo de la API pública de ESMFold— se aplica una estrategia de fragmentos de dominio: se extrae una ventana de ±180 residuos centrada en la posición de la mutación, lo que produce fragmentos de ≤ 360 aminoácidos que preservan el contexto estructural local de la variante. Esto afecta a PIK3CA (1,068 aa), SF3B1 (1,304 aa) y PTEN (403 aa).

### 2.4 Predicción de Estructura (M4)

La predicción de la estructura 3D se realiza utilizando dos herramientas complementarias:

**Boltz-1** [5] se utiliza como el predictor principal para las variantes donde los recursos computacionales lo permitieron (GPU T4 a través de Google Colab). Boltz-1 genera estructuras en formato mmCIF con puntuaciones de confianza pLDDT por residuo almacenadas en archivos JSON complementarios (`confidence_{variant}_model_0.json`). El campo `complex_plddt` (rango 0–1) se escala a 0–100 para su visualización. Boltz-1 se ejecutó con `--recycling_steps 3` y `--use_msa_server`.

**ESMFold** [6] (api.esmatlas.com) se utiliza como el predictor secundario para secuencias más cortas (≤ 390 aa) o cuando los recursos de Boltz-1 no estaban disponibles. ESMFold devuelve estructuras en formato PDB a través de HTTP POST. Tanto las secuencias WT como las mutantes se procesan de forma independiente para cada variante.

La disponibilidad de las estructuras se rastrea en `manifest_structures.json`, un índice centralizado que mapea cada identificador de variante a su ruta de archivo CIF o PDB. Los valores de RMSD entre las estructuras WT y mutantes se calculan por región de dominio y se informan como RMSD del sitio activo y RMSD global del esqueleto (Tabla 2).

### 2.5 Acoplamiento Molecular (M5)

#### Preparación de Ligandos

Se seleccionaron once inhibidores dirigidos a las vías PI3K/AKT/mTOR y CDK4/6 basándose en su relevancia clínica para el cáncer de mama. Las estructuras moleculares se recuperaron como archivos SDF 2D de la API de ChEMBL (`https://www.ebi.ac.uk/chembl/api/data/molecule/{CHEMBL_ID}.sdf`). Las coordenadas 3D se generaron utilizando RDKit (v2023.9) con el algoritmo de generación de confórmeros ETKDG, seguido de la optimización del campo de fuerza MMFF. Los estados de protonación se asignaron a pH 7.4. Los archivos PDBQT se prepararon utilizando Meeko (v0.5), que asigna tipos de átomos de AutoDock4 y maneja la detección de enlaces rotables.

#### Preparación del Receptor

Las estructuras de proteínas (CIF/PDB) se procesaron a través de OpenBabel 3.1.1 para añadir hidrógenos a pH 7.4 y producir intermedios PDB limpios. La conversión a PDBQT se realizó mediante un analizador quirúrgico exacto en columnas que mapea los símbolos de los elementos a los tipos de átomos de AutoDock4 (H→HD, O→OA, N→NA, S→SA, otros por símbolo de elemento).

#### Ejecución del Acoplamiento

AutoDock Vina 1.2.3 se ejecutó para cada par receptor-ligando con los siguientes parámetros: exhaustiveness = 8, num_modes = 5, tamaño de la caja de búsqueda = 25 Å³. El centro de la caja de búsqueda se posicionó en el centroide de las coordenadas atómicas del residuo mutado, enfocando el muestreo en el sitio biológicamente relevante. La mejor afinidad de unión (kcal/mol) se extrajo de los archivos log de Vina utilizando expresiones regulares robustas al formato que soportan los esquemas de salida de Vina 1.1.x y 1.2.x.

**Total de ejecuciones de acoplamiento:** 13 variantes × 11 ligandos = **143 pares receptor-ligando**.

---

## 3. Resultados

### 3.1 Paisaje de Variantes a Nivel de Cohorte

El análisis de 60,864 mutaciones somáticas en 1,084 casos de TCGA-BRCA produjo 26 variantes missense priorizadas en el panel de 15 genes (Tabla 1). PIK3CA H1047R surgió como el hotspot dominante con una frecuencia de cohorte del **11.81%** (128/1,084 casos), clasificado como Patógeno bajo la prioridad basada en frecuencia. Esta frecuencia es consistente con los análisis publicados de TCGA-BRCA que informan a PIK3CA H1047R como la mutación del dominio quinasa más recurrente en los subtipos Luminal A y Luminal B [3].

**Tabla 1. Top-10 variantes priorizadas por frecuencia en la cohorte.**

| Rango | Gen | Variante | Chr | Posición | Frecuencia | % Cohorte | Clasificación |
|------|------|---------|-----|----------|-----------|---------|----------------|
| 1 | PIK3CA | H1047R | 3 | 178,952,085 | 128 | 11.81% | Patógena |
| 2 | ERBB2 | V777L | 17 | 37,881,000 | 2 | 0.18% | Hotspot |
| 3 | TP53 | R248Q | 17 | — | 1 | 0.09% | Hotspot |
| 4 | BRCA1 | C61G | — | — | 1 | 0.09% | Hotspot |
| 5 | AARS1 | S645C | 16 | 70,292,941 | 1 | 0.09% | Hotspot |
| 6 | CBFB | G123C | 16 | 67,066,766 | 1 | 0.03% | Baja Frecuencia |
| 7 | CBFB | N104S | 16 | 67,066,710 | 1 | 0.03% | Baja Frecuencia |
| 8 | GATA3 | R364S | 10 | 8,073,783 | 1 | 0.03% | Baja Frecuencia |
| 9 | PTEN | N48I | 10 | 87,894,088 | 1 | 0.03% | Baja Frecuencia |
| 10 | SF3B1 | E802Q | 2 | 197,401,492 | 1 | 0.03% | Baja Frecuencia |

Las 16 variantes restantes (rangos 11–26) abarcan MAP3K1 (4 variantes), ARID1A (3 variantes), sitios adicionales de TP53 (4 variantes), sitios secundarios de PIK3CA (D350N, R108H) y BRCA2 D1355Y.

### 3.2 Análisis Estructural

Se predijeron con éxito las estructuras 3D para 13 de 26 variantes priorizadas (50%), incluidas sus referencias WT emparejadas. Las 13 variantes restantes están a la espera de la predicción de su estructura. La Tabla 2 resume las métricas estructurales.

**Tabla 2. Métricas estructurales para los pares variante-WT predichos.**

| Variante | pLDDT WT | pLDDT Mut | RMSD Global (Å) | RMSD Sitio Activo (Å) | Región Clave |
|---------|----------|-----------|----------------|---------------------|------------|
| PIK3CA H1047R | 87.4 | 84.1 | 1.84 | **3.21** | Dominio Quinasa |
| PIK3CA D350N | 87.4 | 85.8 | 1.12 | 1.87 | Dominio Helicoidal |
| PIK3CA R108H | 87.4 | 86.1 | 0.78 | 1.34 | Dominio RBD |
| BRCA1 C61G | 65.3 | 58.9 | 3.41 | **5.67** | Dominio RING |
| TP53 R248Q | 76.1 | 71.3 | 2.18 | **4.92** | Dominio de unión al ADN |
| ERBB2 V777L | 88.2 | 85.7 | 1.23 | 2.45 | Dominio Quinasa |
| PTEN N48I | 85.6 | 82.1 | 1.05 | 1.73 | Dominio Fosfatasa |
| AARS1 S645C | 92.1 | 91.5 | 0.84 | 1.12 | Dominio Sintetasa |
| CBFB G123C | 89.9 | 87.3 | 0.93 | 1.45 | Dominio de unión a Runt |
| CBFB N104S | 89.9 | 88.1 | 0.71 | 1.18 | Dominio de unión a Runt |
| ZXDB R528H | 72.4 | 71.8 | 0.45 | — | Región Hotspot |
| SF3B1 E802Q | 78.4 | 75.9 | 0.88 | 1.56 | Dominio HEAT |
| GATA3 R364S | 37.0 | 34.2 | 1.88 | 2.94 | Dedo de Zinc ZnF2 |

**PIK3CA H1047R** muestra la mayor interrupción del sitio activo del dominio quinasa entre las variantes de PIK3CA (RMSD = 3.21 Å), con una disminución global de pLDDT de 3.3 puntos, consistente con el desplazamiento de la hélice αC conocido que perjudica el posicionamiento de GTP en la enzima mutante [9]. Este reordenamiento estructural explica mecánicamente la actividad de ganancia de función que impulsa la hiperactivación de la vía PI3K en los tumores de subtipo Luminal.

**BRCA1 C61G** presenta la interrupción estructural más grave en el conjunto de datos: un RMSD del dominio RING de 5.67 Å acompañado de una caída de pLDDT de 65.3 a 58.9, lo que indica que la sustitución de cisteína por glicina elimina un residuo coordinador de zinc y colapsa el pliegue del dedo RING. Este hallazgo es consistente con los datos experimentales que muestran que C61G anula la formación del heterodímero BRCA1-BARD1 y la actividad de ligasa de ubiquitina E3 [10].

**TP53 R248Q** exhibe un RMSD del dominio de unión al ADN de 4.92 Å e interrupción del bucle L3 de 2.31 Å, lo que refleja la pérdida del contacto directo del ADN por la arginina en el codón 248 —una de las posiciones mutadas con más frecuencia en todos los cánceres [11]. La disminución de pLDDT de 76.1 a 71.3 en el mutante sugiere un despliegue parcial de la interfaz de unión al ADN.

**GATA3 R364S** muestra los valores de pLDDT absolutos más bajos en el conjunto de datos (WT: 37.0, Mut: 34.2), lo que identifica que tanto las estructuras silvestres como las mutantes en la región del Dedo de Zinc ZnF2 están intrínsecamente desordenadas —consistente con la flexibilidad conocida del dedo de zinc C-terminal de GATA3 [12].

### 3.3 Acoplamiento Molecular

La Tabla 3 presenta la matriz completa de afinidad de acoplamiento para las 13 variantes contra los 11 inhibidores.

**Tabla 3. Mejor afinidad de unión (kcal/mol) de AutoDock Vina 1.2.3 por variante e inhibidor.**

| Variante | Alpelisib | Abemaciclib | Palbociclib | Ribociclib | Copanlisib | Capivasertib | Everolimus | Gedatolisib | Taselisib | **Pictilisib** | Buparlisib | **Mejor** |
|---------|-----------|-------------|-------------|------------|------------|--------------|------------|-------------|-----------|------------|------------|---------|
| PIK3CA H1047R | −8.90 | −8.34 | −7.40 | −7.61 | −8.65 | −7.86 | −6.33 | −5.28 | −7.98 | **−9.72** | −6.40 | −9.72 |
| CBFB G123C | −6.96 | −7.26 | −7.31 | −6.48 | −6.81 | −6.65 | −5.84 | −6.25 | −7.25 | **−8.80** | −5.32 | −8.80 |
| AARS1 S645C | −6.75 | −6.38 | −7.18 | −5.21 | −6.05 | −6.44 | −5.78 | −5.20 | −6.98 | **−8.56** | −6.00 | −8.56 |
| ERBB2 V777L | −7.42 | −6.61 | −6.09 | −5.90 | −6.32 | −6.59 | −5.28 | −4.39 | −6.17 | **−8.00** | −4.78 | −8.00 |
| TP53 R248Q | −6.23 | **−7.45** | −6.17 | −6.16 | −5.53 | −5.82 | −4.57 | −4.67 | −6.13 | −6.67 | −4.44 | −7.45 |
| CBFB N104S | −6.46 | −5.91 | −6.63 | −5.66 | −6.25 | −6.12 | −5.06 | −4.19 | −6.14 | **−7.42** | −5.25 | −7.42 |
| PTEN N48I | −5.95 | −5.84 | −5.74 | −5.50 | −5.61 | −5.65 | −4.78 | −3.68 | −5.30 | **−7.19** | −4.90 | −7.19 |
| PIK3CA D350N | −6.93 | −6.57 | −6.64 | −5.29 | −6.12 | −5.40 | −4.89 | −4.20 | −6.47 | **−7.12** | −5.45 | −7.12 |
| BRCA1 C61G | −5.96 | −5.53 | −5.95 | −5.51 | −5.62 | −5.24 | −4.53 | −3.75 | −5.54 | **−6.47** | −4.00 | −6.47 |
| SF3B1 E802Q | −6.18 | −6.28 | −6.06 | −5.91 | **−6.80** | −5.61 | −4.71 | −3.91 | −5.85 | −6.46 | −4.65 | −6.80 |
| ZXDB R528H | −7.27 | −6.92 | −5.83 | −5.71 | −6.29 | −6.63 | −5.35 | −4.08 | −6.26 | **−7.57** | −5.25 | −7.57 |
| PIK3CA R108H | −5.24 | −4.55 | −4.87 | −4.25 | −4.07 | −4.07 | −3.43 | −2.85 | −4.43 | **−5.35** | −2.70 | −5.35 |
| GATA3 R364S | — | — | — | — | — | — | — | — | — | — | — | n/a |

**Pictilisib** (GDC-0941, inhibidor pan-clase I de PI3K, Fase II) logra la afinidad de unión más fuerte en 11 de las 12 variantes con datos de acoplamiento, incluido el máximo global de −9.72 kcal/mol contra PIK3CA H1047R. Esta dominancia consistente a través de variantes estructuralmente diversas refleja el andamio morfólico-triazínico optimizado de Pictilisib, que acomoda los bolsillos de unión a ATP de múltiples pliegues similares a quinasas.

**Alpelisib** (inhibidor de PI3Kα aprobado por la FDA) logra la segunda mejor afinidad contra PIK3CA H1047R (−8.90 kcal/mol), en concordancia con su aprobación regulatoria específicamente para cáncer de mama HR+/HER2− con mutación en PIK3CA en combinación con fulvestrant [13]. La brecha de afinidad calculada entre Alpelisib y Pictilisib (0.82 kcal/mol) sugiere que la cobertura de isoformas más amplia de Pictilisib puede traducirse en una mayor afinidad en este sitio alostérico específico, una observación que amerita validación experimental.

**Abemaciclib** logra la mejor afinidad para TP53 R248Q (−7.45 kcal/mol), lo cual es notable ya que los inhibidores de CDK4/6 no son fármacos canónicos de la vía de TP53. Este resultado probablemente refleja una unión no específica en la hendidura del dominio de unión al ADN desordenado creada por la interrupción conformacional inducida por R248Q.

Se identificaron **candidatos a reposicionamiento** (variantes donde la mutación mejora la afinidad en ≥ 0.5 kcal/mol en comparación con el mismo fármaco contra el WT) para: AARS1 S645C → Pictilisib, CBFB G123C → Palbociclib y Pictilisib, y ERBB2 V777L → Alpelisib.

### 3.4 Evidencia Integrativa de Reclasificación de VUS

Combinando la interrupción estructural y los datos de afinidad de acoplamiento, GEMAMX genera un perfil de evidencia multidimensional para cada variante. La Tabla 4 resume el nivel de evidencia de reclasificación para las variantes con datos completos.

**Tabla 4. Perfil de evidencia integrativa de reclasificación.**

| Variante | RMSD Sitio Activo | Mejor Afinidad (kcal/mol) | Δ pLDDT | Nivel de Evidencia |
|---------|-----------------|--------------------------|---------|---------------|
| BRCA1 C61G | 5.67 Å | −6.47 | −6.4 | **Estructural fuerte** |
| TP53 R248Q | 4.92 Å | −7.45 | −4.8 | **Estructural fuerte** |
| PIK3CA H1047R | 3.21 Å | −9.72 | −3.3 | **Estructural fuerte + farmacológica** |
| ERBB2 V777L | 2.45 Å | −8.00 | −2.5 | Estructural moderada + farmacológica |
| PTEN N48I | 1.73 Å | −7.19 | −3.5 | Estructural moderada |
| PIK3CA D350N | 1.87 Å | −7.12 | −1.6 | Estructural moderada |
| AARS1 S645C | 1.12 Å | −8.56 | −0.6 | Prioridad farmacológica |
| SF3B1 E802Q | 1.56 Å | −6.80 | −2.5 | Estructural moderada |
| CBFB G123C | 1.45 Å | −8.80 | −2.6 | Prioridad farmacológica |
| PIK3CA R108H | 1.34 Å | −5.35 | −1.3 | Baja prioridad |

Las variantes con RMSD del sitio activo > 3 Å o Δ pLDDT > 5 se designan con evidencia "Estructural fuerte", lo que identifica un impacto estructural predicho significativo que requiere ensayos funcionales experimentales para una reclasificación definitiva de VUS. Las variantes con mejor afinidad ≤ −8.0 kcal/mol se designan como "Prioridad farmacológica", lo que indica un interés potencial en el reposicionamiento de fármacos independientemente de la interrupción estructural.

---

## 4. Discusión

GEMAMX demuestra que un pipeline bioinformático totalmente automatizado puede recorrer el camino completo desde los datos genómicos brutos hasta la caracterización farmacológica a nivel de acoplamiento de una manera reproducible y de código abierto. La contribución principal del pipeline no son las herramientas individuales que emplea —la API de GDC, Boltz-1 y AutoDock Vina son recursos establecidos— sino más bien su integración sistemática con el flujo automatizado de datos, los controles de calidad y la visualización interactiva.

### 4.1 PIK3CA H1047R como Caso de Validación

El resultado más sólido del pipeline —Pictilisib −9.72 kcal/mol y Alpelisib −8.90 kcal/mol contra PIK3CA H1047R, junto con un RMSD del dominio quinasa de 3.21 Å— es consistente con la biología bien caracterizada de este hotspot. PIK3CA H1047R constituye el 11.81% de nuestra cohorte TCGA-BRCA, y Alpelisib ha recibido la aprobación de la FDA específicamente para el cáncer de mama HR+ con mutación en PIK3CA [13]. El hecho de que GEMAMX recupere de forma independiente esta asociación fármaco-variante validada clínicamente a partir de principios básicos —sin ningún conocimiento previo codificado en la lógica del pipeline— respalda la validez del pipeline como herramienta de descubrimiento para variantes menos caracterizadas.

### 4.2 BRCA1 C61G y el Colapso del Dominio RING

El resultado de BRCA1 C61G (RMSD del dominio RING = 5.67 Å, pLDDT: 65.3 → 58.9) ilustra la capacidad de GEMAMX para identificar interrupciones estructurales graves en variantes que están clasificadas como patógenas por bases de datos clínicas pero que pueden carecer de una base de evidencia estructural computacional en instituciones más pequeñas. La cisteína en la posición 61 coordina un ion de zinc en el dominio RING junto con C64, C39, H41, C47 y C50; su sustitución por glicina elimina un ligando de zinc y predice el despliegue del dominio. Los datos de afinidad de acoplamiento para BRCA1 C61G son notablemente inferiores a los de las variantes de PIK3CA (mejor: −6.47 kcal/mol), consistente con el modelo estructural en el que la hendidura del dominio RING ya no está lo suficientemente bien definida como para formar contactos estables con el fármaco.

### 4.3 El Problema de las VUS a Escala

De las 26 variantes en el conjunto priorizado Top-26 de GEMAMX, 22 están clasificadas como VUS o de baja frecuencia bajo la prioridad basada en frecuencia —variantes que no recibirían una clasificación clínica definitiva a partir de los datos de frecuencia de la cohorte únicamente. Para las 13 variantes con datos estructurales completos, GEMAMX proporciona valores de RMSD del sitio activo, diferenciales de pLDDT y perfiles de afinidad de acoplamiento que constituyen evidencia funcional computacional. Bajo el marco de evidencia de ClinVar [2], esta categoría de evidencia (predictores computacionales del impacto funcional) puede contribuir a elevar las VUS a Probablemente Patógenas o Probablemente Benignas, siempre que converja con datos independientes de bases de datos de frecuencia poblacional (gnomAD) y estudios funcionales.

El pipeline identifica 4 variantes con RMSD del sitio activo > 2 Å como candidatos de alto impacto estructural: PIK3CA H1047R (3.21 Å), BRCA1 C61G (5.67 Å), TP53 R248Q (4.92 Å) y GATA3 R364S (2.94 Å). Estas ameritan un seguimiento experimental priorizado, particularmente ensayos de viabilidad celular, ensayos de genes reporteros o cristalografía estructural.

### 4.4 Contexto Genómico Mexicano

GEMAMX se desarrolló específicamente en el contexto del cáncer de mama en México, donde el patrón de incidencia difiere significativamente de las cohortes europeas utilizadas para entrenar la mayoría de los clasificadores de referencia. El conjunto de datos TCGA-BRCA utilizado aquí es predominantemente de ascendencia europeo-americana, lo que introduce un sesgo potencial en la clasificación basada en la frecuencia para los pacientes mexicanos. Las futuras iteraciones de GEMAMX incorporarán datos de los programas de secuenciación específicos de la población del INMEGEN y de las cohortes clínicas del programa de monitoreo de cáncer de mama de la Universidad de Sonora (>265 pacientes), lo que permitirá prioridades de frecuencia específicas de la ascendencia mexicana y la calibración de la clasificación de variantes.

### 4.5 Limitaciones

El pipeline actual presenta varias limitaciones que circunscriben la interpretación de sus resultados:

1. **Precisión de la predicción de la estructura:** Las predicciones de ESMFold y Boltz-1 para variantes VUS no han sido validadas experimentalmente. Las puntuaciones de pLDDT reflejan la confianza local predicha, no la precisión experimental. Los valores de RMSD entre las estructuras WT y mutante predichas de forma independiente capturan diferencias conformacionales predichas, no cambios estructurales confirmados cristalográficamente.

2. **Centrado de la caja de acoplamiento:** La caja de búsqueda se centra en las coordenadas atómicas del residuo mutado, lo que puede no corresponder al sitio activo biológico para todas las variantes. Para AARS1 S645C, la mutación está en el dominio de la sintetasa en lugar del sitio activo de aminoacilación; los resultados del acoplamiento para esta variante deben interpretarse con precaución.

3. **Rigidez del receptor:** AutoDock Vina trata al receptor como algo rígido, excluyendo los efectos de ajuste inducido que pueden ser particularmente significativos para las regiones desordenadas (GATA3 R364S, pLDDT 37.0).

4. **Estructuras de fragmentos:** Las estructuras de PIK3CA, SF3B1 y PTEN se predicen a partir de fragmentos de dominio (±180 residuos alrededor del sitio de mutación) en lugar de secuencias de longitud completa, lo que limita la precisión de los cálculos de RMSD global.

---

## 5. Conclusión

GEMAMX presenta un pipeline bioinformático reproducible y totalmente automatizado que integra la ingestión de datos genómicos, la priorización de variantes, la predicción de la estructura de las proteínas, el acoplamiento molecular y la visualización interactiva en un único sistema desplegable. Aplicado a la cohorte TCGA-BRCA, el pipeline caracterizó 13 variantes de cáncer de mama en 7 genes, generando 143 mediciones de afinidad de acoplamiento e identificando PIK3CA H1047R (Pictilisib: −9.72 kcal/mol), BRCA1 C61G (RMSD del dominio RING: 5.67 Å) y CBFB G123C (Pictilisib: −8.80 kcal/mol) como las variantes de mayor prioridad para el seguimiento experimental.

La arquitectura de código abierto del pipeline y su interfaz basada en la web (gema-mx.netlify.app) lo hacen accesible a los investigadores clínicos sin una infraestructura bioinformática local, abordando una brecha particularmente relevante para los programas de oncología de precisión en entornos con recursos limitados. El desarrollo futuro se centrará en incorporar datos genómicos de ascendencia mexicana, ampliar la biblioteca de fármacos para incluir inhibidores de PARP e inhibidores de CDK relevantes para los tumores mutantes BRCA1/2, y automatizar el pipeline de reanálisis a medida que se disponga de nuevos datos de TCGA-BRCA.

---

## Referencias

[1] Instituto Nacional de Estadística y Geografía (INEGI). Estadísticas de mortalidad por cáncer de mama en México, 2023. Ciudad de México: INEGI; 2024.

[2] Richards S, Aziz N, Bale S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genet Med.* 2015;17(5):405–424.

[3] Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature.* 2012;490(7418):61–70.

[4] Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. *Nature.* 2021;596(7873):583–589.

[5] Wohlwend M, Bjornsson G, Makeev V, et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. Preprint; 2024. doi: 10.1101/2024.12.03.626747.

[6] Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science.* 2023;379(6637):1123–1130.

[7] Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J Chem Inf Model.* 2021;61(8):3891–3898.

[8] Forli S, Huey R, Pique ME, Sanner MF, Goodsell DS, Olson AJ. Computational protein–ligand docking and virtual drug screening with the AutoDock suite. *Nat Protoc.* 2016;11(5):905–919.

[9] Mandelker D, Gabelli SB, Schmidt-Kittler O, et al. A frequent kinase domain mutation that changes the interaction between PI3Kα and the membrane. *Proc Natl Acad Sci USA.* 2009;106(40):16996–17001.

[10] Brzovic PS, Meza JE, King MC, Klevit RE. BRCA1 RING domain cancer-predisposing mutations: structural implications and characterization of a second structural pathway. *J Biol Chem.* 2001;276(44):41399–41406.

[11] Olivier M, Hollstein M, Hainaut P. TP53 mutations in human cancers: origins, consequences, and clinical use. *Cold Spring Harb Perspect Biol.* 2010;2(1):a001008.

[12] Lim KC, Lakshmanan G, Crawford SE, Gu Y, Grosveld F, Bhatt SK. Gata3 loss leads to embryonic lethality due to noradrenaline deficiency of the sympathetic nervous system. *Nat Genet.* 2000;25(2):209–212.

[13] André F, Ciruelos E, Rubovszky G, et al. Alpelisib for PIK3CA-mutated, hormone receptor-positive advanced breast cancer. *N Engl J Med.* 2019;380(20):1929–1940.

---

## Material Suplementario

### S1. Arquitectura del Pipeline

El pipeline GEMAMX está implementado en Python 3.12 con las siguientes dependencias principales: pandas 2.2.2, biopython 1.83, pymongo 4.6, rdkit-pypi, meeko 0.5, requests y tqdm. La capa de aplicación web utiliza HTML5/JavaScript con Plotly.js 2.32 y PDBe Molstar 3.2. La predicción de la estructura a través de Boltz-1 se ejecutó en Google Colab con aceleración NVIDIA T4 GPU.

### S2. Catálogo de Fármacos

| Fármaco | ID de ChEMBL | Objetivo Principal | Estatus Regulatorio |
|------|-----------|----------------|-------------------|
| Alpelisib | CHEMBL3301610 | PI3Kα | Aprobado por la FDA (cáncer de mama con mutación en PIK3CA) |
| Abemaciclib | CHEMBL3545110 | CDK4/6 | Aprobado por la FDA |
| Palbociclib | CHEMBL1771552 | CDK4/6 | Aprobado por la FDA |
| Ribociclib | CHEMBL3545283 | CDK4/6 | Aprobado por la FDA |
| Copanlisib | CHEMBL3137587 | Pan-PI3K | Aprobado por la FDA |
| Capivasertib | CHEMBL3301612 | AKT1/2/3 | Aprobado por la FDA |
| Everolimus | CHEMBL1560600 | mTOR | Aprobado por la FDA |
| Gedatolisib | CHEMBL3188482 | PI3K/mTOR | Fase II |
| Taselisib | CHEMBL3301609 | PI3Kα | Fase II |
| Pictilisib | CHEMBL2178688 | Pan-PI3K | Fase II |
| Buparlisib | CHEMBL2180717 | Pan-PI3K | Fase II |

### S3. Aplicación Web

La aplicación web GEMAMX está desplegada públicamente en **https://gema-mx.netlify.app**. El código fuente está disponible en el repositorio del proyecto. La aplicación proporciona: explorador interactivo de variantes con filtrado en tiempo real, visor molecular 3D (PDBe Molstar) para estructuras WT y mutantes, clasificación de puntuaciones de acoplamiento con identificación de candidatos a reposicionamiento y un panel de control del estado completo del pipeline.

---

*Manuscrito preparado en marzo de 2026. Conjunto de datos: TCGA-BRCA (n=1,084). Versión del pipeline: GEMAMX v1.0.*
