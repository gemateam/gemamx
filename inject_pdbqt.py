import json

with open('01_tcga_exploration.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

new_markdown_cell = {
    'cell_type': 'markdown',
    'metadata': {},
    'source': [
        '## 6.6 Módulo Chemoinformático: Preparación para Vina (.pdbqt)\n',
        '\n',
        'Convertimos automáticamente los ligandos (.sdf) y receptores (.pdb/.cif) al formato `.pdbqt` requerido por AutoDock Vina utilizando `Meeko` y `OpenBabel`.'
    ]
}

new_code_cell = {
    'cell_type': 'code',
    'execution_count': None,
    'metadata': {},
    'outputs': [],
    'source': [
        '!sudo apt-get install openbabel -y --quiet\n',
        '!pip install meeko --quiet\n',
        '\n',
        'import os, glob\n',
        'import subprocess\n',
        'from meeko import MoleculePreparation\n',
        '\n',
        'LIGANDS_DIR = "data/docking/ligands"\n',
        'RECEPTORS_DIR = "data/docking/receptors"\n',
        'STRUCTURES_DIR = "data/structures"\n',
        'os.makedirs(RECEPTORS_DIR, exist_ok=True)\n',
        '\n',
        'print("--- PREPARANDO LIGANDOS (Meeko) ---")\n',
        'sdf_files = glob.glob(f"{LIGANDS_DIR}/*.sdf")\n',
        'for sdf in sdf_files:\n',
        '    try:\n',
        '        from rdkit import Chem\n',
        '        mol = next(Chem.SDMolSupplier(sdf, removeHs=False))\n',
        '        if mol:\n',
        '            preparator = MoleculePreparation()\n',
        '            preparator.prepare(mol)\n',
        '            pdbqt_string = preparator.write_pdbqt_string()\n',
        '            \n',
        '            out_path = sdf.replace(".sdf", ".pdbqt")\n',
        '            with open(out_path, "w") as f:\n',
        '                f.write(pdbqt_string)\n',
        '            print(f"Ligando preparado: {os.path.basename(out_path)}")\n',
        '    except Exception as e:\n',
        '        print(f"Error al preparar {os.path.basename(sdf)}: {e}")\n',
        '\n',
        'print("\\n--- PREPARANDO RECEPTORES (OpenBabel) ---")\n',
        'receptor_files = glob.glob(f"{STRUCTURES_DIR}/**/*.cif", recursive=True) + glob.glob(f"{STRUCTURES_DIR}/**/*.pdb", recursive=True)\n',
        'for filepath in receptor_files:\n',
        '    name = os.path.basename(filepath).split(".")[0]\n',
        '    if name == "placeholder": continue\n',
        '    \n',
        '    ext = filepath.split(".")[-1]\n',
        '    out_path = f"{RECEPTORS_DIR}/{name}.pdbqt"\n',
        '    \n',
        '    try:\n',
        '        subprocess.run(["obabel", "-i", ext, filepath, "-o", "pdbqt", "-O", out_path, "-xr"], check=False, capture_output=True)\n',
        '        if os.path.exists(out_path):\n',
        '            print(f"Receptor preparado: {os.path.basename(out_path)}")\n',
        '    except Exception as e:\n',
        '        print(f"Error en obabel: {e}")\n',
        '\n',
        'print("\\n--- PREPARACIÓN COMPLETADA ---")\n'
    ]
}

idx = 16
if '7. Exportación' in nb['cells'][idx]['source'][0]:
    nb['cells'].insert(idx, new_code_cell)
    nb['cells'].insert(idx, new_markdown_cell)

    with open('01_tcga_exploration.ipynb', 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
    print('Notebook successfuly updated using Python.')
else:
    print('Failed to find insertion point.')
