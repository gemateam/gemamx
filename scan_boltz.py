"""
scan_boltz.py — Escanea todas las carpetas boltz_results_* y reporta cuáles tienen CIF.
También actualiza manifest_structures.json con todo lo que existe en disco.
"""
import os, json, glob

STRUCTURES_DIR = "data/structures"
MANIFEST_PATH  = "data/manifest_structures.json"

print("=" * 55)
print("SCAN DE ESTRUCTURAS BOLTZ")
print("=" * 55)

# 1. Escanear CIF de Boltz
boltz_found = {}
for cif in glob.glob(f"{STRUCTURES_DIR}/**/predictions/**/*.cif", recursive=True):
    name = os.path.basename(os.path.dirname(cif))
    size_kb = os.path.getsize(cif) / 1024
    # Leer pLDDT del confidence_*.json si existe
    conf_path = cif.replace("_model_0.cif", "").replace(name, f"confidence_{name}") + "_model_0.json"
    plddt = "--"
    if os.path.exists(conf_path):
        try:
            c = json.load(open(conf_path))
            plddt = round(c.get("complex_plddt", 0) * 100, 1)
        except: pass
    boltz_found[name] = {"path": cif.replace("\\", "/"), "size_kb": round(size_kb, 1), "plddt": plddt}
    print(f"  ✓ {name:35s} {size_kb:7.1f} KB   pLDDT: {plddt}")

print(f"\n  Total CIF Boltz: {len(boltz_found)}")

# 2. Escanear PDB directos (AlphaFold / ESMFold)
print("\n" + "=" * 55)
print("PDB DIRECTOS (AlphaFold / ESMFold)")
print("=" * 55)
pdb_found = {}
for pdb in glob.glob(f"{STRUCTURES_DIR}/*.pdb"):
    name = os.path.basename(pdb).replace(".pdb", "")
    size_kb = os.path.getsize(pdb) / 1024
    if "placeholder" in name or "_MUT" in name or size_kb < 1:
        continue
    pdb_found[name] = {"path": pdb.replace("\\", "/"), "size_kb": round(size_kb, 1)}
    print(f"  ✓ {name:35s} {size_kb:7.1f} KB")

for cif in glob.glob(f"{STRUCTURES_DIR}/*.cif"):
    name = os.path.basename(cif).replace(".cif", "")
    size_kb = os.path.getsize(cif) / 1024
    pdb_found[name] = {"path": cif.replace("\\", "/"), "size_kb": round(size_kb, 1)}
    print(f"  ✓ {name:35s} {size_kb:7.1f} KB  (AlphaFold)")

print(f"\n  Total PDB/CIF directos: {len(pdb_found)}")

# 3. Cargar manifest actual y actualizarlo
manifest = {}
if os.path.exists(MANIFEST_PATH):
    with open(MANIFEST_PATH) as f:
        manifest = json.load(f)

# Agregar Boltz
for name, info in boltz_found.items():
    manifest[name] = info["path"]

# Agregar PDB/CIF directos
for name, info in pdb_found.items():
    manifest[name] = info["path"]

# Guardar
with open(MANIFEST_PATH, "w") as f:
    json.dump(manifest, f, indent=4)

print("\n" + "=" * 55)
print("MANIFEST ACTUALIZADO")
print("=" * 55)
total = len(manifest)
con_archivo = sum(1 for v in manifest.values() if v)
sin_archivo = total - con_archivo
print(f"  Total entradas:  {total}")
print(f"  Con archivo:     {con_archivo}")
print(f"  Sin archivo (null): {sin_archivo}")
print(f"\n  Guardado en: {MANIFEST_PATH}")

# 4. Resumen de lo que falta
print("\n" + "=" * 55)
print("VARIANTES SIN ESTRUCTURA (null en manifest)")
print("=" * 55)
for key, val in sorted(manifest.items()):
    if val is None:
        print(f"  ✗ {key}")
