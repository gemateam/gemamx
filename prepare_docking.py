"""
prepare_docking.py
==================
Pipeline completo para docking molecular real con AutoDock Vina.

Pasos:
  1. Descarga SDF de fármacos desde ChEMBL API (por nombre/SMILES)
  2. Convierte SDF → PDBQT con OpenBabel
  3. Convierte receptores PDB/CIF → PDBQT con OpenBabel
  4. Calcula centro de masa de cada receptor (proxy del sitio activo)
  5. Corre AutoDock Vina para cada par receptor-ligando
  6. Parsea resultados y guarda vina_affinities.json
  7. Actualiza structural_metrics.json con scores reales

Requisitos:
  pip install requests biopython rdkit-pypi meeko
  conda install -c conda-forge autodock-vina openbabel  (o apt-get en Linux)

Uso:
  python prepare_docking.py
"""

import os, json, subprocess, glob, re, time, requests
import sys
import shutil
import warnings
from pathlib import Path

# Silenciar advertencias de Meeko/Deprecation
warnings.filterwarnings("ignore", category=DeprecationWarning)

# ── Directorios ────────────────────────────────────────────────────────────
DOCKING_DIR   = Path("data/docking")
LIGANDS_DIR   = DOCKING_DIR / "ligands"
RECEPTORS_DIR = DOCKING_DIR / "receptors"
RESULTS_DIR   = DOCKING_DIR / "results"
STRUCTURES_DIR = Path("data/structures")

for d in [LIGANDS_DIR, RECEPTORS_DIR, RESULTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Configuración de Binarios ──────────────────────────────────────────────
# Usar vina123.exe (v1.2.3) para mayor estabilidad
VINA_EXE = "vina123.exe"
OBABEL_EXE = r"D:\Program Files\OpenBabel-3.1.1\obabel.exe"

# ── Catálogo de fármacos con ChEMBL IDs ───────────────────────────────────
# ChEMBL IDs verificados para cada fármaco
DRUG_CATALOG = {
    "Alpelisib":    {"chembl_id": "CHEMBL3301610", "target": "PI3Kα",     "status": "FDA"},
    "Abemaciclib":  {"chembl_id": "CHEMBL3545110", "target": "CDK4/6",    "status": "FDA"},
    "Palbociclib":  {"chembl_id": "CHEMBL1771552", "target": "CDK4/6",    "status": "FDA"},
    "Ribociclib":   {"chembl_id": "CHEMBL3545283", "target": "CDK4/6",    "status": "FDA"},
    "Copanlisib":   {"chembl_id": "CHEMBL3137587", "target": "Pan-PI3K",  "status": "FDA"},
    "Capivasertib": {"chembl_id": "CHEMBL3301612", "target": "AKT1/2/3",  "status": "FDA"},
    "Everolimus":   {"chembl_id": "CHEMBL1560600", "target": "mTOR",      "status": "FDA"},
    "Gedatolisib":  {"chembl_id": "CHEMBL3188482", "target": "PI3K/mTOR", "status": "Trial"},
    "Taselisib":    {"chembl_id": "CHEMBL3301609", "target": "PI3Kα",     "status": "Trial"},
    "Pictilisib":   {"chembl_id": "CHEMBL2178688", "target": "PI3K",      "status": "Trial"},
    "Buparlisib":   {"chembl_id": "CHEMBL2180717", "target": "Pan-PI3K",  "status": "Trial"},
}

# ── Variantes con estructura disponible para docking ──────────────────────
# Solo las que tienen PDB/CIF en disco
VARIANTS_FOR_DOCKING = [
    "PIK3CA_H1047R",
    "PIK3CA_D350N",
    "PIK3CA_R108H",
    "ERBB2_V777L",
    "TP53_R248Q",
    "BRCA1_C61G",
    "PTEN_N48I",
    "AARS1_S645C",
    "CBFB_G123C",
    "CBFB_N104S",
    "GATA3_R364S",
    "SF3B1_E802Q",
    "ZXDB_R528H",
]

# ── Parámetros de Vina ─────────────────────────────────────────────────────
VINA_EXHAUSTIVENESS = 1    # 1 para screening rápido, 8 para producción
VINA_NUM_MODES      = 5    # Poses a generar por ligando
VINA_BOX_SIZE       = 25   # Å — caja de búsqueda centrada en centro de masa


# ==========================================================================
# PASO 1 — Descargar SDF de ChEMBL
# ==========================================================================
def download_sdf_from_chembl(drug_name, chembl_id):
    """Descarga el SDF 2D del compuesto desde ChEMBL API."""
    out_path = LIGANDS_DIR / f"{drug_name}.sdf"
    if out_path.exists() and out_path.stat().st_size > 100:
        print(f"  ↳ {drug_name}.sdf ya existe")
        return str(out_path)

    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.sdf"
    print(f"  Descargando {drug_name} ({chembl_id})...")
    try:
        r = requests.get(url, timeout=30)
        if r.ok and len(r.text) > 50:
            out_path.write_text(r.text, encoding="utf-8")
            print(f"  ✓ {drug_name}.sdf ({len(r.text)//1024} KB)")
            return str(out_path)
        else:
            print(f"  ✗ Error {r.status_code} descargando {drug_name}")
    except Exception as e:
        print(f"  ✗ Excepción: {e}")
    return None


# ==========================================================================
# PASO 2 — Convertir SDF → PDBQT (ligando)
# ==========================================================================
def sdf_to_pdbqt_ligand(sdf_path, drug_name):
    """Convierte SDF 2D → PDBQT 3D con RDKit y Meeko."""
    out_path = LIGANDS_DIR / f"{drug_name}.pdbqt"
    if out_path.exists() and out_path.stat().st_size > 100:
        print(f"  ↳ {drug_name}.pdbqt ya existe")
        return str(out_path)

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.MolStandardize import rdMolStandardize
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
        if mol is None:
            print(f"  ✗ RDKit no pudo leer {sdf_path}")
            return None
            
        # Seleccionar fragmento más grande (eliminar contraiones/sales)
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)
            
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)

        # Preparar con Meeko (estándar para Vina 1.2+)
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        
        # Guardar usando el writer de Meeko (maneja ramas correctamente)
        pdbqt_string = preparator.write_pdbqt_string()
        
        # Asegurar saltos de línea consistentes
        with open(out_path, "w", encoding="ascii") as f:
            f.write(pdbqt_string)
        
        if out_path.exists() and out_path.stat().st_size > 100:
            print(f"  ✓ {drug_name}.pdbqt (Meeko)")
            return str(out_path)
    except Exception as e:
        print(f"  ✗ Error Meeko: {e}")
    return None


# ==========================================================================
# PASO 3 — Preparar receptor PDBQT
# ==========================================================================
def write_pdbqt_line(record, serial, name, resname, chain, resseq, x, y, z, occ, temp, charge, type_ad4):
    """Construye una línea PDBQT concatenando strings por columnas exactas."""
    # Columnas:
    # 1-6   : Record (ATOM/HETATM)
    # 7-11  : Serial
    # 13-16 : Name
    # 18-20 : Resname
    # 22    : Chain
    # 23-26 : ResSeq
    # 31-38 : X
    # 39-46 : Y
    # 47-54 : Z
    # 55-60 : Occ
    # 61-66 : Temp
    # 71-76 : Charge
    # 78-79 : Type
    
    line =  f"{record:<6}"[:6]                   # 1-6
    line += f"{serial:>5d}"[:5]                  # 7-11
    line += " "                                    # 12
    line += f"{name:<4}"[:4]                     # 13-16
    line += " "                                    # 17
    line += f"{resname:<3}"[:3]                  # 18-20
    line += " "                                    # 21
    line += f"{chain:1}"[:1]                     # 22
    line += f"{resseq:>4d}"[:4]                  # 23-26
    line += "    "                                 # 27-30
    line += f"{x:8.3f}"[:8]                      # 31-38
    line += f"{y:8.3f}"[:8]                      # 39-46
    line += f"{z:8.3f}"[:8]                      # 47-54
    line += f"{occ:6.2f}"[:6]                    # 55-60
    line += f"{temp:6.2f}"[:6]                   # 61-66
    line += "    "                                 # 67-70
    line += f"{charge:6.3f}"[:6]                 # 71-76
    line += " "                                    # 77
    line += f"{type_ad4.upper():<2}"[:2]         # 78-79
    return line

def pdb_to_pdbqt_simple(pdb_path, out_path):
    """Convierte PDB a PDBQT de forma quirúrgica (especial para receptores rígidos)."""
    try:
        output = []
        with open(pdb_path, "r", encoding="ascii", errors="ignore") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        name = line[12:16].strip()
                        resname = line[17:20].strip()
                        chain = line[21].strip() or "A"
                        resseq = int(line[22:26] or 1)
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        
                        element = line[76:78].strip().upper()
                        if not element: element = name[:1].upper() or "C"
                        
                        # Mapeo simple de elementos a tipos AD4
                        if element == "H": ad_type = "HD"
                        else: ad_type = element[:2]
                        
                        output.append(write_pdbqt_line("ATOM", len(output)+1, name, resname, chain, resseq, x, y, z, 1.0, 0.0, 0.0, ad_type))
                    except: continue
        
        with open(out_path, "wb") as f:
            for l in output: f.write((l + "\n").encode("ascii"))
            f.write(b"TER\nEND\n")
        
        if out_path.exists() and out_path.stat().st_size > 10000:
            return True
    except Exception as e:
        print(f"      ✗ Error quirúrgico: {e}")
    return False

def structure_to_pdbqt_receptor(variant_id, manifest):
    """Convierte PDB del receptor → PDBQT usando Meeko."""
    out_path = RECEPTORS_DIR / f"{variant_id}.pdbqt"
    if out_path.exists() and out_path.stat().st_size > 5000:
        print(f"  ↳ {variant_id}.pdbqt ya existe")
        return str(out_path)

    struct_path = manifest.get(variant_id)
    if not struct_path or not Path(struct_path).exists():
        print(f"  ✗ Sin estructura para {variant_id}")
        return None

    print(f"  Preparando receptor {variant_id}...")
    try:
        # 1. Limpiar con obabel a PDB (esto funciona en OB 3.1.1)
        tmp_pdb = str(RECEPTORS_DIR / f"{variant_id}_clean.pdb")
        subprocess.run([
            OBABEL_EXE, struct_path,
            "-O", tmp_pdb,
            "-xr", "--addh", "-p", "7.4"
        ], capture_output=True, timeout=60)
        
        if not Path(tmp_pdb).exists(): return None

        # 2. Convertir PDB limpio a PDBQT de forma quirúrgica
        if pdb_to_pdbqt_simple(tmp_pdb, out_path):
            print(f"  ✓ {variant_id}.pdbqt (Surgical, {out_path.stat().st_size//1024} KB)")
            return str(out_path)
        else:
            print(f"  ✗ Error en conversión quirúrgica")
            
    except Exception as e:
        print(f"  ✗ Error en receptor: {e}")
    return None


# ==========================================================================
# PASO 4 — Calcular centro de masa del receptor
# ==========================================================================
def get_mutation_center(pdb_path, variant_id):
    """Calcula el centro basado en el residuo mutado (ej. H1047)."""
    try:
        # Extraer posición del ID (ej. PIK3CA_H1047R -> 1047)
        match = re.search(r'([A-Z])(\d+)[A-Z]', variant_id)
        if not match:
            return get_center_of_mass(pdb_path)
            
        pos = int(match.group(2))
        coords = []
        
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    # Offset de ResSeq en PDB es 22-26
                    try:
                        res_seq = int(line[22:26].strip())
                        if res_seq == pos:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            coords.append([x, y, z])
                    except: continue
        
        if coords:
            import numpy as np
            cx, cy, cz = np.mean(coords, axis=0)
            return float(cx), float(cy), float(cz)
            
    except Exception as e:
        print(f"      ✗ Error buscando mutación en {variant_id}: {e}")
        
    return get_center_of_mass(pdb_path)

def get_center_of_mass(pdb_path):
    """Calcula el centro de masa de toda la proteína (fallback)."""
    try:
        import numpy as np
        coords = []
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except: continue
        if coords:
            cx, cy, cz = np.mean(coords, axis=0)
            return float(cx), float(cy), float(cz)
    except: pass
    return 0.0, 0.0, 0.0


# ==========================================================================
# PASO 5 — Correr AutoDock Vina
# ==========================================================================
def run_vina(variant_id, receptor_pdbqt, ligand_pdbqt, drug_name, cx, cy, cz):
    """Ejecuta AutoDock Vina y retorna la mejor afinidad (kcal/mol)."""
    drug_safe = drug_name.replace(" ", "_")
    out_pdbqt = RESULTS_DIR / f"{variant_id}_{drug_safe}_out.pdbqt"
    log_file  = RESULTS_DIR / f"{variant_id}_{drug_safe}_log.txt"
    config_file = RESULTS_DIR / f"{variant_id}_{drug_safe}_config.txt"

    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 100:
        # Ya existe — parsear resultado previo
        return parse_vina_log(str(log_file))

    try:
        # Usar argumentos directos para Vina 1.2.3 (más robusto en Windows)
        cmd = [
            str(Path(VINA_EXE).absolute()),
            "--receptor",     str(Path(receptor_pdbqt).absolute()),
            "--ligand",       str(Path(ligand_pdbqt).absolute()),
            "--center_x",     f"{cx:.3f}",
            "--center_y",     f"{cy:.3f}",
            "--center_z",     f"{cz:.3f}",
            "--size_x",       str(VINA_BOX_SIZE),
            "--size_y",       str(VINA_BOX_SIZE),
            "--size_z",       str(VINA_BOX_SIZE),
            "--exhaustiveness", str(VINA_EXHAUSTIVENESS),
            "--num_modes",    str(VINA_NUM_MODES),
            "--out",          str(out_pdbqt)
        ]
        
        # Guardar comando para debug
        config_file.write_text(" ".join(cmd), encoding="utf-8")
        
        # Ejecutar capturando el output en el log manual
        with open(log_file, "w") as f_log:
            result = subprocess.run(cmd, stdout=f_log, stderr=subprocess.PIPE, text=True, timeout=600)

        if out_pdbqt.exists():
            affinity = parse_vina_log(str(log_file))
            return affinity
        else:
            print(f"\n      ✗ Vina falló para {drug_name} vs {variant_id}")
            if result.stderr: print(f"      ERR: {result.stderr.strip()[:200]}")
            return None
    except FileNotFoundError:
        print("    ✗ AutoDock Vina no encontrado. Instala con: conda install -c conda-forge autodock-vina")
    except subprocess.TimeoutExpired:
        print(f"    ✗ Timeout en {variant_id} vs {drug_name}")
    except Exception as e:
        print(f"    ✗ Error: {e}")
    return None


def parse_vina_log(log_path):
    """Extrae el mejor score de afinidad (kcal/mol) de un log de Vina."""
    if not Path(log_path).exists(): return None
    try:
        content = Path(log_path).read_text(encoding="utf-8", errors="ignore")
        # Patrón robusto para tabla de Vina 1.1.x y 1.2.x
        # Ejemplo: "   1       -8.903          0          0"
        match = re.search(r'^\s*1\s+([-+]?\d+\.\d+)', content, re.MULTILINE)
        if match:
            return float(match.group(1))
        
        # Fallback: buscar "Estimated Free Energy of Binding"
        match = re.search(r'Estimated Free Energy of Binding\s*:\s*([-+]?\d+\.\d+)', content)
        if match:
            return float(match.group(1))
    except Exception as e:
        print(f"      ✗ Error parseando log: {e}")
    return None


# ==========================================================================
# MAIN
# ==========================================================================
def main():
    print("=" * 60)
    print("PIPELINE DE DOCKING MOLECULAR — GEN-MAMA-MX")
    print("=" * 60)

    # Cargar manifest
    manifest_path = "data/manifest_structures.json"
    with open(manifest_path) as f:
        manifest = json.load(f)

    # ── PASO 1+2: Descargar y convertir ligandos ───────────────────────────
    print("\n📦 PASO 1-2: Preparando ligandos (SDF → PDBQT)")
    print("-" * 60)
    ligand_pdbqts = {}
    for drug_name, info in DRUG_CATALOG.items():
        print(f"\n[{drug_name}]")
        sdf = download_sdf_from_chembl(drug_name, info["chembl_id"])
        if sdf:
            pdbqt = sdf_to_pdbqt_ligand(sdf, drug_name)
            if pdbqt:
                ligand_pdbqts[drug_name] = pdbqt
        time.sleep(0.5)  # Rate limit ChEMBL

    print(f"\n✓ Ligandos preparados: {len(ligand_pdbqts)}/{len(DRUG_CATALOG)}")
    if not ligand_pdbqts:
        print("✗ Sin ligandos disponibles. Verifica OpenBabel y conexión a internet.")
        return

    # ── PASO 3: Preparar receptores ────────────────────────────────────────
    print("\n\n🧬 PASO 3: Preparando receptores (PDB/CIF → PDBQT)")
    print("-" * 60)
    receptor_pdbqts = {}
    for variant_id in VARIANTS_FOR_DOCKING:
        print(f"\n[{variant_id}]")
        pdbqt = structure_to_pdbqt_receptor(variant_id, manifest)
        if pdbqt:
            receptor_pdbqts[variant_id] = pdbqt

    print(f"\n✓ Receptores preparados: {len(receptor_pdbqts)}/{len(VARIANTS_FOR_DOCKING)}")

    # ── PASO 4+5: Correr Vina ─────────────────────────────────────────────
    print("\n\n⚗️  PASO 4-5: Ejecutando AutoDock Vina")
    print("-" * 60)

    vina_results  = {}   # { variant_id: { drug_name: affinity } }
    total_runs    = len(receptor_pdbqts) * len(ligand_pdbqts)
    completed     = 0

    for variant_id, receptor_pdbqt in receptor_pdbqts.items():
        # Centrar en sitio de mutación si es posible
        print(f"\n[{variant_id}]")
        cx, cy, cz = get_mutation_center(receptor_pdbqt, variant_id)
        print(f"  Centro (sitio mutación): ({cx:.1f}, {cy:.1f}, {cz:.1f})")
        vina_results[variant_id] = {}

        for drug_name, ligand_pdbqt in ligand_pdbqts.items():
            print(f"  vs {drug_name}...", end=" ", flush=True)
            affinity = run_vina(
                variant_id, receptor_pdbqt, ligand_pdbqt, drug_name, cx, cy, cz
            )
            vina_results[variant_id][drug_name] = affinity
            completed += 1
            if affinity is not None:
                print(f"{affinity} kcal/mol")
            else:
                print("N/A")

        # Guardar resultados parciales (por si se interrumpe)
        with open("data/docking/vina_affinities.json", "w") as f:
            json.dump(vina_results, f, indent=4)

    # ── PASO 6: Guardar resultados finales ────────────────────────────────
    print("\n\n💾 PASO 6: Guardando resultados")
    print("-" * 60)
    with open("data/docking/vina_affinities.json", "w") as f:
        json.dump(vina_results, f, indent=4)
    print(f"✓ data/docking/vina_affinities.json ({completed} pares evaluados)")

    # ── PASO 7: Actualizar structural_metrics.json ────────────────────────
    print("\n📊 PASO 7: Actualizando structural_metrics.json con scores reales")
    print("-" * 60)

    with open("data/structural_metrics.json") as f:
        metrics = json.load(f)

    for variant_id, drug_scores in vina_results.items():
        valid_scores = {k: v for k, v in drug_scores.items() if v is not None}
        if not valid_scores:
            continue
        best_drug  = min(valid_scores, key=valid_scores.get)
        best_score = valid_scores[best_drug]
        worst_drug = max(valid_scores, key=valid_scores.get)

        if variant_id not in metrics:
            metrics[variant_id] = {}

        metrics[variant_id]["best_vina_drug"]    = best_drug
        metrics[variant_id]["best_vina_score"]   = best_score
        metrics[variant_id]["vina_scores"]        = valid_scores
        print(f"  {variant_id}: mejor={best_drug} ({best_score} kcal/mol)")

    with open("data/structural_metrics.json", "w") as f:
        json.dump(metrics, f, indent=4)
    print("✓ structural_metrics.json actualizado con scores Vina reales")

    # ── Resumen final ─────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("RESUMEN FINAL")
    print("=" * 60)
    print(f"  Ligandos procesados:  {len(ligand_pdbqts)}")
    print(f"  Receptores procesados: {len(receptor_pdbqts)}")
    print(f"  Pares evaluados:      {completed}")
    print(f"\nPoses guardadas en: data/docking/results/")
    print("Para ver las poses en la app: click en '🔬 Pose' en el ranking")
    print("\nNOTA: Actualiza docking.js para leer vina_scores reales desde")
    print("      structural_metrics.json en lugar de calcular valores simulados.")


if __name__ == "__main__":
    main()
