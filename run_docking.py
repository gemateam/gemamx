"""
run_docking.py
==============
Rerun limpio del docking con diagnóstico previo.
- Verifica version de Vina y formato de log
- Calcula centro de masa REAL desde coordenadas del receptor
- Maneja ambos formatos de log (Vina 1.1.x y 1.2.x)
- Genera vina_affinities.json y actualiza structural_metrics.json

Uso:
  (venv) python run_docking.py
"""

import os, json, subprocess, re, glob
from pathlib import Path

RECEPTORS_DIR = Path("data/docking/receptors")
LIGANDS_DIR   = Path("data/docking/ligands")
RESULTS_DIR   = Path("data/docking/results")
RESULTS_DIR.mkdir(exist_ok=True)

# ── Verificar Vina disponible ──────────────────────────────────────────────
def check_vina():
    for exe in ["vina", "vina.exe", "vina_1.2.3_windows_x86_64", "vina_1.2.5_win.exe"]:
        try:
            r = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=5)
            version_line = (r.stdout + r.stderr).strip().split("\n")[0]
            print(f"✓ Vina encontrado: {exe} — {version_line}")
            return exe
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    print("✗ AutoDock Vina no encontrado en PATH")
    print("  Descarga desde: https://github.com/ccsb-scripps/AutoDock-Vina/releases")
    print("  Extrae vina.exe a D:\\ProyectosIA\\GEMAMX\\ o agrégalo al PATH")
    return None

# ── Calcular centro de masa real del receptor ──────────────────────────────
def center_of_mass(pdbqt_path):
    xs, ys, zs = [], [], []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except (ValueError, IndexError):
                    continue
    if not xs:
        return 0.0, 0.0, 0.0
    n = len(xs)
    return sum(xs)/n, sum(ys)/n, sum(zs)/n

# ── Parsear afinidad del log (soporta Vina 1.1.x y 1.2.x) ─────────────────
def parse_affinity(log_path):
    """
    Vina 1.1.x formato:
        -----+------------+----------+----------
           1       -8.2      0.000      0.000
    Vina 1.2.x formato:
        mode |   affinity | dist from best mode
           1        -8.3          0          0
    """
    if not Path(log_path).exists():
        return None
    content = Path(log_path).read_text()

    # Buscar cualquier linea que empiece con espacios + "1" + numero negativo
    patterns = [
        r'^\s+1\s+([-\d]+\.\d+)',           # Vina 1.2.x
        r'^\s+1\s+\|\s+([-\d]+\.\d+)',      # Vina 1.2.x con pipe
        r'^   1\s+([-\d]+\.\d+)',            # Vina 1.1.x
    ]
    for line in content.split('\n'):
        for pat in patterns:
            m = re.match(pat, line)
            if m:
                try:
                    val = float(m.group(1))
                    if -20 < val < 0:   # Rango razonable de afinidades
                        return val
                except ValueError:
                    pass
    return None

# ── Correr Vina para un par receptor-ligando ──────────────────────────────
def run_vina_pair(vina_exe, variant_id, receptor_pdbqt, ligand_pdbqt, drug_name, cx, cy, cz):
    drug_safe = drug_name.replace(" ", "_")
    out_pdbqt = RESULTS_DIR / f"{variant_id}_{drug_safe}_out.pdbqt"
    log_file  = RESULTS_DIR / f"{variant_id}_{drug_safe}_log.txt"

    # Si ya existe con resultado valido, no recorrer
    if log_file.exists():
        existing = parse_affinity(str(log_file))
        if existing is not None:
            return existing  # Ya tiene resultado real

    cmd = [
        vina_exe,
        "--receptor",       str(receptor_pdbqt),
        "--ligand",         str(ligand_pdbqt),
        "--center_x",       f"{cx:.3f}",
        "--center_y",       f"{cy:.3f}",
        "--center_z",       f"{cz:.3f}",
        "--size_x",         "25",
        "--size_y",         "25",
        "--size_z",         "25",
        "--exhaustiveness", "8",
        "--num_modes",      "5",
        "--out",            str(out_pdbqt),
        "--log",            str(log_file),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        print(f"    ✗ Timeout")
        return None
    except FileNotFoundError:
        print(f"    ✗ Vina no encontrado: {vina_exe}")
        return None

    affinity = parse_affinity(str(log_file))
    return affinity

# ── MAIN ──────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("RUN_DOCKING — Rerun limpio con diagnóstico")
    print("=" * 60)

    # 1. Verificar Vina
    vina_exe = check_vina()
    if not vina_exe:
        return

    # 2. Detectar receptores y ligandos disponibles
    receptors = {}
    for pdbqt in RECEPTORS_DIR.glob("*.pdbqt"):
        # Solo receptores con contenido real (>1 KB)
        if pdbqt.stat().st_size > 1000 and "_tmp" not in pdbqt.stem and "_clean" not in pdbqt.stem:
            receptors[pdbqt.stem] = str(pdbqt)

    ligands = {}
    for pdbqt in LIGANDS_DIR.glob("*.pdbqt"):
        if pdbqt.stat().st_size > 100:
            ligands[pdbqt.stem] = str(pdbqt)

    print(f"\nReceptores válidos ({len(receptors)}): {list(receptors.keys())}")
    print(f"Ligandos válidos  ({len(ligands)}):   {list(ligands.keys())}")

    if not receptors:
        print("\n✗ Sin receptores PDBQT válidos. Corre prepare_docking.py primero.")
        return
    if not ligands:
        print("\n✗ Sin ligandos PDBQT. Corre prepare_docking.py primero.")
        return

    # 3. Docking
    print(f"\n{'='*60}")
    print("EJECUTANDO DOCKING")
    print(f"{'='*60}")

    results = {}
    total = len(receptors) * len(ligands)
    done  = 0

    for variant_id, receptor_path in sorted(receptors.items()):
        cx, cy, cz = center_of_mass(receptor_path)
        print(f"\n[{variant_id}] centro=({cx:.1f}, {cy:.1f}, {cz:.1f})")
        results[variant_id] = {}

        for drug_name, ligand_path in sorted(ligands.items()):
            print(f"  {drug_name:15s} ...", end=" ", flush=True)
            aff = run_vina_pair(vina_exe, variant_id, receptor_path, ligand_path, drug_name, cx, cy, cz)
            results[variant_id][drug_name] = aff
            done += 1
            print(f"{aff if aff else 'FALLÓ'} kcal/mol  [{done}/{total}]")

        # Guardar parcial
        with open("data/docking/vina_affinities.json", "w") as f:
            json.dump(results, f, indent=4)

    # 4. Resumen y actualizar structural_metrics
    print(f"\n{'='*60}")
    print("RESULTADOS")
    print(f"{'='*60}")

    with open("data/structural_metrics.json") as f:
        metrics = json.load(f)

    success = 0
    for variant_id, drug_scores in results.items():
        valid = {k: v for k, v in drug_scores.items() if v is not None}
        success += len(valid)
        if valid:
            best  = min(valid, key=valid.get)
            score = valid[best]
            print(f"  {variant_id:20s} → mejor: {best} ({score} kcal/mol)  [{len(valid)}/{len(drug_scores)} OK]")
            if variant_id not in metrics:
                metrics[variant_id] = {}
            metrics[variant_id]["best_vina_drug"]  = best
            metrics[variant_id]["best_vina_score"] = score
            metrics[variant_id]["vina_scores"]      = valid
        else:
            print(f"  {variant_id:20s} → TODOS FALLARON")

    with open("data/structural_metrics.json", "w") as f:
        json.dump(metrics, f, indent=4)

    print(f"\n✓ {success}/{total} pares con score real")
    print(f"✓ vina_affinities.json actualizado")
    print(f"✓ structural_metrics.json actualizado")
    print(f"\nHaz Ctrl+Shift+R en el navegador para ver los scores reales en la app.")

if __name__ == "__main__":
    main()
