import json, os
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO

def fix_receptors():
    manifest_path = "data/manifest_structures.json"
    with open(manifest_path) as f:
        manifest = json.load(f)

    problematic = ["CBFB_G123C", "CBFB_N104S", "GATA3_R364S", "ZXDB_R528H"]
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    for variant in problematic:
        cif_path = manifest.get(variant)
        if not cif_path or not os.path.exists(cif_path):
            # Buscar en boltz_results si el manifest tiene una ruta vieja o nula
            search_dir = Path("data/structures")
            cif_files = list(search_dir.glob(f"boltz_results_{variant}/**/*.cif"))
            if cif_files:
                cif_path = str(cif_files[0])
            else:
                print(f"  ✗ No se encontró CIF para {variant}")
                continue

        out_pdb = f"data/structures/{variant}.pdb"
        print(f"  Convirtiendo {variant}: {cif_path} -> {out_pdb}")
        
        try:
            structure = parser.get_structure(variant, cif_path)
            io.set_structure(structure)
            io.save(out_pdb)
            if os.path.exists(out_pdb) and os.path.getsize(out_pdb) > 10000:
                print(f"  ✓ Exito: {variant}.pdb pronto ({os.path.getsize(out_pdb)//1024} KB)")
                # Actualizar manifest
                manifest[variant] = out_pdb
            else:
                print(f"  ✗ Error: PDB generado es sospechosamente pequeño")
        except Exception as e:
            print(f"  ✗ Falllo en {variant}: {e}")

    # Guardar manifest actualizado
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=4)
    print("\n✓ Manifest actualizado.")

if __name__ == "__main__":
    fix_receptors()
