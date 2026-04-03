# =============================================================================
# GEMA MX — Guía de Instalación Completa
# Genomic Analysis of Mammary cancer in Mexico
# =============================================================================
# Sistema operativo: Windows 10/11 x64
# Directorio del proyecto: D:\ProyectosIA\GEMAMX
# =============================================================================


## ÍNDICE
  1. Python 3.12
  2. MongoDB Community
  3. OpenBabel 3.1.1
  4. AutoDock Vina 1.2.3
  5. Entorno virtual de Python
  6. Librerías Python (requirements.txt)
  7. Google Colab (para predicción estructural con Boltz-1)
  8. Verificación final del entorno


═══════════════════════════════════════════════════════════════════════════════
## PASO 1 — Python 3.12
═══════════════════════════════════════════════════════════════════════════════

1. Descarga Python 3.12 desde:
   https://www.python.org/downloads/release/python-3120/
   → Archivo: python-3.12.x-amd64.exe

2. Durante la instalación:
   ✅ Marca "Add Python 3.12 to PATH"
   ✅ Marca "Install for all users"
   → Click en "Customize installation"
   ✅ Marca "pip"
   ✅ Marca "py launcher"

3. Verifica la instalación:
   python --version
   pip --version

   Resultado esperado:
   Python 3.12.x
   pip 24.x.x


═══════════════════════════════════════════════════════════════════════════════
## PASO 2 — MongoDB Community 7.0
═══════════════════════════════════════════════════════════════════════════════

MongoDB se usa durante la fase de ingesta (M1) para almacenar mutaciones
y el ranking de variantes. No se necesita en tiempo de ejecución de la app web.

1. Descarga MongoDB Community Server desde:
   https://www.mongodb.com/try/download/community
   → Versión: 7.0.x  → Platform: Windows  → Package: msi

2. Ejecuta el instalador .msi:
   ✅ Selecciona "Complete"
   ✅ Marca "Install MongoDB as a Service"
   ✅ Marca "Run service as Network Service user"
   ✅ Marca "Install MongoDB Compass" (GUI opcional pero útil)

3. Verifica que el servicio está corriendo:
   Abre PowerShell como administrador y ejecuta:

   Get-Service -Name MongoDB

   Resultado esperado:
   Status: Running

4. Si el servicio no inicia automáticamente:
   net start MongoDB

5. Verifica la conexión:
   mongosh
   → Debes ver el prompt: test>
   → Escribe: exit


═══════════════════════════════════════════════════════════════════════════════
## PASO 3 — OpenBabel 3.1.1
═══════════════════════════════════════════════════════════════════════════════

OpenBabel se usa para convertir estructuras PDB de receptores a formato PDBQT
que requiere AutoDock Vina.

1. Descarga OpenBabel 3.1.1 para Windows desde:
   https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1
   → Archivo: OpenBabel-3.1.1-x64.exe

2. Ejecuta el instalador:
   → Directorio de instalación: D:\Program Files\OpenBabel-3.1.1
   ✅ Marca "Add to PATH" durante la instalación

3. Verifica la instalación:
   obabel --version

   Resultado esperado:
   Open Babel 3.1.1

NOTA: El path del ejecutable está hardcodeado en prepare_docking.py como:
      OBABEL_EXE = r"D:\Program Files\OpenBabel-3.1.1\obabel.exe"
      Si instalas en otro directorio, actualiza esa variable.


═══════════════════════════════════════════════════════════════════════════════
## PASO 4 — AutoDock Vina 1.2.3
═══════════════════════════════════════════════════════════════════════════════

AutoDock Vina es el motor de docking molecular. Se ejecuta localmente
para todas las simulaciones de unión proteína-fármaco.

1. Descarga AutoDock Vina 1.2.3 para Windows desde:
   https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.3
   → Archivo: vina_1.2.3_windows_x86_64.exe

2. Renombra el archivo descargado a:
   vina123.exe

3. Copia vina123.exe al directorio del proyecto:
   D:\ProyectosIA\GEMAMX\vina123.exe

   (El ejecutable ya debe estar en ese path — verifica que existe)

4. Verifica que funciona desde la carpeta del proyecto:
   cd D:\ProyectosIA\GEMAMX
   .\vina123.exe --version

   Resultado esperado:
   AutoDock Vina 1.2.3


═══════════════════════════════════════════════════════════════════════════════
## PASO 5 — Entorno Virtual de Python
═══════════════════════════════════════════════════════════════════════════════

Se recomienda un entorno virtual para aislar las dependencias del proyecto
y evitar conflictos con otras versiones de librerías.

1. Abre una terminal en el directorio del proyecto:
   cd D:\ProyectosIA\GEMAMX

2. Crea el entorno virtual:
   python -m venv venv

3. Activa el entorno virtual:
   .\venv\Scripts\activate

   Resultado esperado: verás (venv) al inicio del prompt.

4. Actualiza pip dentro del entorno:
   python -m pip install --upgrade pip

IMPORTANTE: Siempre activa el entorno antes de ejecutar cualquier
script del proyecto:
   cd D:\ProyectosIA\GEMAMX
   .\venv\Scripts\activate


═══════════════════════════════════════════════════════════════════════════════
## PASO 6 — Librerías Python
═══════════════════════════════════════════════════════════════════════════════

Con el entorno virtual activado, instala todas las dependencias:

   pip install -r requirements.txt

Esto instalará:
  • numpy 1.26.4           — cómputo numérico
  • pandas 2.2.2           — manejo de datos tabulares
  • requests               — llamadas HTTP a APIs externas
  • tqdm                   — barras de progreso
  • biopython 1.83         — SeqIO, MutableSeq (mutaciones en FASTA)
  • pymongo                — cliente MongoDB
  • rdkit-pypi             — quimioinformática (generación 3D de ligandos)
  • meeko                  — preparación de PDBQT para Vina 1.2+
  • python-dotenv          — gestión de variables de entorno

Verificación individual de instalaciones críticas:

   python -c "import numpy; print('numpy', numpy.__version__)"
   python -c "import pandas; print('pandas', pandas.__version__)"
   python -c "from Bio import SeqIO; print('biopython OK')"
   python -c "from pymongo import MongoClient; print('pymongo OK')"
   python -c "from rdkit import Chem; print('rdkit OK')"
   python -c "from meeko import MoleculePreparation; print('meeko OK')"


═══════════════════════════════════════════════════════════════════════════════
## PASO 7 — Google Colab (Predicción Estructural con Boltz-1)
═══════════════════════════════════════════════════════════════════════════════

Boltz-1 requiere GPU para tiempos razonables de predicción. Se ejecuta
en Google Colab con GPU T4 gratuita.

PREREQUISITO: Cuenta de Google con acceso a Google Drive.

─── Configuración inicial (solo se hace una vez) ────────────────────────────

1. Crea la carpeta del proyecto en tu Google Drive:
   Mi unidad / GEMAMX /

2. Sube el notebook al Drive o ábrelo en Colab:
   Archivo: 01_tcga_exploration.ipynb
   → Abre Google Colab: https://colab.research.google.com
   → File → Open notebook → Google Drive → selecciona el .ipynb

3. Activa la GPU:
   Entorno de ejecución → Cambiar tipo de entorno de ejecución
   → Acelerador de hardware: T4 GPU
   → Guardar

─── Ejecución del notebook ──────────────────────────────────────────────────

4. Ejecuta la Celda 0 (Configuración del Entorno):
   → Instala: numpy 1.26.4, pandas 2.2.2, pytorch-lightning 2.1.4, boltz
   → Monta Google Drive
   → Crea los symlinks /content/data → Drive/GEMAMX/data

   ⚠️  IMPORTANTE: Después de esta celda, REINICIA el entorno:
   Entorno de ejecución → Reiniciar entorno

5. Ejecuta la Celda 0b (Reconexión post-reinicio):
   → Restaura el montaje de Drive y los symlinks
   → Verifica versiones instaladas

6. Ejecuta las celdas en orden:
   Celda 1-2  → Descarga mutaciones TCGA-BRCA (GDC API)
   Celda 3-4  → Genera FASTA WT y mutados (UniProt API)
   Celda 5    → Predicción estructural con Boltz-1 (≈ 5-15 min por variante)
   Celda 6.7  → Docking molecular con AutoDock Vina (versión Linux)
   Celda final → Exporta gemamx_data.zip

─── Sincronización con el proyecto local ────────────────────────────────────

7. Descarga gemamx_data.zip desde Colab o Google Drive

8. Descomprime en el directorio del proyecto:
   → Extrae en: D:\ProyectosIA\GEMAMX\
   → Esto actualizará la carpeta data/ con los archivos generados

9. Actualiza el manifest de estructuras localmente:
   cd D:\ProyectosIA\GEMAMX
   .\venv\Scripts\activate
   python scan_boltz.py

10. Abre index.html en el navegador para ver los resultados:
    → Abre D:\ProyectosIA\GEMAMX\index.html con Chrome o Edge
    → Los datos se cargan automáticamente desde data/


═══════════════════════════════════════════════════════════════════════════════
## PASO 8 — Verificación Final del Entorno
═══════════════════════════════════════════════════════════════════════════════

Ejecuta este checklist en orden para confirmar que todo está instalado
correctamente:

─── Verificaciones de sistema ───────────────────────────────────────────────

   python --version
   # Esperado: Python 3.12.x

   pip --version
   # Esperado: pip 24.x.x

   obabel --version
   # Esperado: Open Babel 3.1.1

   cd D:\ProyectosIA\GEMAMX && .\vina123.exe --version
   # Esperado: AutoDock Vina 1.2.3

   Get-Service -Name MongoDB  (PowerShell)
   # Esperado: Status Running

─── Verificaciones de Python (con venv activo) ──────────────────────────────

   cd D:\ProyectosIA\GEMAMX
   .\venv\Scripts\activate

   python -c "
   import numpy, pandas, requests, tqdm
   from Bio import SeqIO
   from pymongo import MongoClient
   from rdkit import Chem
   from meeko import MoleculePreparation
   print('✓ Todas las librerías importadas correctamente')
   print(f'  numpy   {numpy.__version__}')
   print(f'  pandas  {pandas.__version__}')
   "

─── Verificación de MongoDB ─────────────────────────────────────────────────

   python -c "
   from pymongo import MongoClient
   c = MongoClient('mongodb://localhost:27017/', serverSelectionTimeoutMS=3000)
   c.server_info()
   print('✓ MongoDB conectado correctamente')
   "

─── Verificación de estructura de datos ─────────────────────────────────────

   # Verifica que existen los archivos JSON necesarios para la app web:
   dir D:\ProyectosIA\GEMAMX\data\top_variants.json
   dir D:\ProyectosIA\GEMAMX\data\raw_mutations.json
   dir D:\ProyectosIA\GEMAMX\data\mutation_stats.json
   dir D:\ProyectosIA\GEMAMX\data\structural_metrics.json
   dir D:\ProyectosIA\GEMAMX\data\manifest_structures.json

   # Si algún archivo falta, ejecuta primero:
   python ingest_data.py       # genera top_variants.json y raw_mutations.json
   python rebuild_stats.py     # genera mutation_stats.json
   python scan_boltz.py        # genera/actualiza manifest_structures.json


═══════════════════════════════════════════════════════════════════════════════
## RESUMEN DE VERSIONES
═══════════════════════════════════════════════════════════════════════════════

  Software externo:
  ┌─────────────────────┬──────────────────────────────────────────────────┐
  │ Componente          │ Versión / Fuente                                 │
  ├─────────────────────┼──────────────────────────────────────────────────┤
  │ Python              │ 3.12.x — python.org                              │
  │ MongoDB Community   │ 7.0.x — mongodb.com                              │
  │ OpenBabel           │ 3.1.1 — github.com/openbabel                     │
  │ AutoDock Vina       │ 1.2.3 — github.com/ccsb-scripps/AutoDock-Vina    │
  │ Google Colab        │ Gratuito — colab.research.google.com             │
  └─────────────────────┴──────────────────────────────────────────────────┘

  Librerías Python (ver requirements.txt para versiones exactas):
  ┌─────────────────────┬──────────────────────────────────────────────────┐
  │ Librería            │ Uso en el pipeline                               │
  ├─────────────────────┼──────────────────────────────────────────────────┤
  │ numpy 1.26.4        │ Cómputo numérico (compatible con Boltz-1)        │
  │ pandas 2.2.2        │ Procesamiento del dataset TCGA-BRCA              │
  │ requests            │ APIs: cBioPortal, GDC, UniProt, ChEMBL          │
  │ tqdm                │ Barras de progreso en descarga masiva            │
  │ biopython 1.83      │ Lectura FASTA, sustitución de aminoácidos        │
  │ pymongo             │ Almacén temporal de mutaciones e índices         │
  │ rdkit-pypi          │ Generación de coordenadas 3D de fármacos        │
  │ meeko               │ Preparación de ligandos PDBQT para Vina 1.2+    │
  │ python-dotenv       │ Variables de entorno                             │
  │ boltz (Colab only)  │ Predicción estructural 3D con GPU T4             │
  └─────────────────────┴──────────────────────────────────────────────────┘

  Dependencias de la app web (cargadas desde CDN, sin instalación):
  ┌─────────────────────┬──────────────────────────────────────────────────┐
  │ Librería            │ Versión                                          │
  ├─────────────────────┼──────────────────────────────────────────────────┤
  │ Plotly.js           │ 2.32.0                                           │
  │ PDBe Molstar        │ 3.2.0                                            │
  └─────────────────────┴──────────────────────────────────────────────────┘


═══════════════════════════════════════════════════════════════════════════════
## SOLUCIÓN DE PROBLEMAS FRECUENTES
═══════════════════════════════════════════════════════════════════════════════

PROBLEMA: "numpy binary incompatibility" al importar pandas
SOLUCIÓN:
   pip install numpy==1.26.4 pandas==2.2.2 --force-reinstall

───────────────────────────────────────────────────────────────────────────────

PROBLEMA: "ModuleNotFoundError: No module named 'rdkit'"
SOLUCIÓN:
   pip install rdkit-pypi
   # Si falla, instalar con conda:
   conda install -c conda-forge rdkit

───────────────────────────────────────────────────────────────────────────────

PROBLEMA: "EPERM: operation not permitted" al escribir archivos con MongoDB
SOLUCIÓN: Asegúrate de correr la terminal como administrador o verificar
          que el directorio del proyecto tiene permisos de escritura.

───────────────────────────────────────────────────────────────────────────────

PROBLEMA: Vina devuelve error "Could not find receptor"
SOLUCIÓN: Los paths deben ser absolutos. run_docking.py usa
          Path(receptor_pdbqt).absolute() automáticamente.
          Verifica que los archivos .pdbqt existen en:
          D:\ProyectosIA\GEMAMX\data\docking\receptors\

───────────────────────────────────────────────────────────────────────────────

PROBLEMA: Boltz-1 falla con "LeafSpec is deprecated" en Colab
SOLUCIÓN: Reinstala pytorch-lightning con la versión compatible:
   !pip install pytorch-lightning==2.1.4 lightning==2.1.4 --force-reinstall
   Luego reinicia el entorno de ejecución.

───────────────────────────────────────────────────────────────────────────────

PROBLEMA: La app web no carga datos (muestra "—" en todos los KPIs)
SOLUCIÓN: Verifica que existen los archivos JSON en data/:
   dir D:\ProyectosIA\GEMAMX\data\*.json
   Si faltan, ejecuta el pipeline en orden:
   1. python ingest_data.py
   2. python fetch_fasta.py
   3. python scan_boltz.py
   Luego recarga index.html en el navegador.


═══════════════════════════════════════════════════════════════════════════════
Desarrollado en la Universidad de Sonora (UNISON)
Google.org Impact Challenge: AI for Science
═══════════════════════════════════════════════════════════════════════════════
