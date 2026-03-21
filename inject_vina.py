import json

with open('01_tcga_exploration.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

new_markdown_cell = {
    'cell_type': 'markdown',
    'metadata': {},
    'source': [
        '## 6.7 Módulo Computacional: AutoDock Vina (Batch)\n',
        '\n',
        'Ejecución masiva de simulaciones de docking molecular entre cada variante priorizada y los fármacos correspondientes descargados de ChEMBL. Extraemos la mejor afinidad (kcal/mol) de cada pareo.'
    ]
}

new_code_cell = {
    'cell_type': 'code',
    'execution_count': None,
    'metadata': {},
    'outputs': [],
    'source': [
        '!sudo apt-get install autodock-vina -y --quiet\n',
        '\n',
        'import os, glob, subprocess, json\n',
        '\n',
        'LIGANDS_DIR = "data/docking/ligands"\n',
        'RECEPTORS_DIR = "data/docking/receptors"\n',
        'OUT_DIR = "data/docking/results"\n',
        'os.makedirs(OUT_DIR, exist_ok=True)\n',
        '\n',
        'ligands = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")\n',
        'receptors = glob.glob(f"{RECEPTORS_DIR}/*.pdbqt")\n',
        '\n',
        'results = {}\n',
        '\n',
        'print("--- EJECUTANDO AUTODOCK VINA ---")\n',
        'for rec in receptors:\n',
        '    rec_name = os.path.basename(rec).replace(".pdbqt", "")\n',
        '    gene = rec_name.split("_")[0]\n',
        '    \n',
        '    if rec_name not in results: results[rec_name] = {}\n',
        '    \n',
        '    for lig in ligands:\n',
        '        lig_name = os.path.basename(lig).replace(".pdbqt", "")\n',
        '        if not lig_name.startswith(f"{gene}_"): continue\n',
        '        \n',
        '        drug_name = lig_name.replace(f"{gene}_", "")\n',
        '        out_pdbqt = f"{OUT_DIR}/{rec_name}_{drug_name}_out.pdbqt"\n',
        '        log_txt = f"{OUT_DIR}/{rec_name}_{drug_name}_log.txt"\n',
        '        \n',
        '        print(f"Docking: {rec_name} vs {drug_name}...")\n',
        '        \n',
        '        cmd = [\n',
        '            "vina", "--receptor", rec, "--ligand", lig,\n',
        '            "--center_x", "0", "--center_y", "0", "--center_z", "0",\n',
        '            "--size_x", "40", "--size_y", "40", "--size_z", "40",\n',
        '            "--exhaustiveness", "8", "--out", out_pdbqt, "--log", log_txt\n',
        '        ]\n',
        '        \n',
        '        subprocess.run(cmd, capture_output=True)\n',
        '        \n',
        '        affinity = "N/A"\n',
        '        if os.path.exists(log_txt):\n',
        '            with open(log_txt, "r") as f:\n',
        '                for line in f:\n',
        '                    if "   1 " in line:\n',
        '                        parts = [p for p in line.split(" ") if p]\n',
        '                        if len(parts) >= 2:\n',
        '                            affinity = float(parts[2])\n',
        '                        break\n',
        '        results[rec_name][drug_name] = affinity\n',
        '        print(f"  -> Score: {affinity} kcal/mol")\n',
        '\n',
        'with open("data/docking/vina_affinities.json", "w") as f:\n',
        '    json.dump(results, f, indent=4)\n',
        '\n',
        'print("\\n--- DOCKING COMPLETADO ---")\n'
    ]
}

idx = -1
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'markdown' and len(cell['source']) > 0 and '7. Exportación' in cell['source'][0]:
        idx = i
        break

if idx != -1:
    nb['cells'].insert(idx, new_code_cell)
    nb['cells'].insert(idx, new_markdown_cell)

    with open('01_tcga_exploration.ipynb', 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
    print('Notebook successfuly updated using Python dynamic insertion.')
else:
    print('Failed to find dynamic insertion point.')
