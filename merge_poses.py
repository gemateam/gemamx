"""
merge_poses.py
==============
Genera PDB combinados (proteína + ligando) con formato estricto PDB 3.30
compatible con PDBe Molstar.

Columnas PDB:
 1- 6  Record name
 7-11  Serial number
12     space
13-16  Atom name (alineado según longitud)
17     Alternate location
18-20  Residue name
21     space
22     Chain ID
23-26  Residue sequence number
27     Code for insertion of residues
28-30  spaces
31-38  X orthogonal coordinate
39-46  Y orthogonal coordinate
47-54  Z orthogonal coordinate
55-60  Occupancy
61-66  Temperature factor
77-78  Element symbol (right-justified)
"""

import os, json, re
from pathlib import Path

RECEPTORS_DIR = Path("data/docking/receptors")
RESULTS_DIR   = Path("data/docking/results")
COMPLEX_DIR   = Path("data/docking/complex")
COMPLEX_DIR.mkdir(parents=True, exist_ok=True)

# Mapa tipos AD4 → elemento PDB
AD4_TO_ELEM = {
    'C': 'C', 'A': 'C', 'CA': 'C',
    'N': 'N', 'NA': 'N', 'NS': 'N',
    'O': 'O', 'OA': 'O', 'OS': 'O',
    'H': 'H', 'HD': 'H', 'HS': 'H',
    'S': 'S', 'SA': 'S',
    'F': 'F', 'Cl': 'CL', 'CL': 'CL',
    'Br': 'BR', 'BR': 'BR', 'I': 'I',
    'P': 'P', 'FE': 'FE', 'ZN': 'ZN',
    'MG': 'MG', 'MN': 'MN',
}

def elem_from_ad4(ad4):
    return AD4_TO_ELEM.get(ad4.strip(), ad4.strip()[:1].upper() or 'C')

def format_atom_name(name, element):
    """
    Formatea el nombre del átomo en 4 caracteres según PDB 3.30:
    - Elemento de 1 letra: _ nombre _ _ (ej ' CA ')
    - Elemento de 2 letras: nombre _ _ _ (ej 'FE  ')
    """
    name = name.strip()
    elem = element.strip()
    if len(elem) == 1:
        # Átomos normales: nombre de 1-3 chars alineado a partir de col 14
        if len(name) <= 3:
            return f' {name:<3}'
        else:
            return name[:4]
    else:
        # Elementos de 2 letras (FE, ZN, CL, etc): col 13-14
        return f'{name:<4}'

def make_pdb_line(record, serial, atom_name, resname, chain, resseq, x, y, z, element):
    """Construye una línea PDB de 80 caracteres con columnas exactas."""
    aname = format_atom_name(atom_name, element)
    elem  = f'{element.strip():>2}'
    line = (
        f'{record:<6}'        # 1-6   record
        f'{serial:5d}'        # 7-11  serial
        f' '                  # 12    space
        f'{aname}'            # 13-16 atom name
        f' '                  # 17    altLoc
        f'{resname:>3}'       # 18-20 resname
        f' '                  # 21    space
        f'{chain:1}'          # 22    chain
        f'{resseq:4d}'        # 23-26 resseq
        f'    '               # 27-30 (iCode + spaces)
        f'{x:8.3f}'           # 31-38 x
        f'{y:8.3f}'           # 39-46 y
        f'{z:8.3f}'           # 47-54 z
        f'  1.00'             # 55-60 occupancy
        f'  0.00'             # 61-66 tempFactor
        f'          '         # 67-76 spaces
        f'{elem}'             # 77-78 element
    )
    return line[:80]  # Truncar a 80 caracteres exactos


def read_receptor_atoms(pdbqt_path):
    """Lee receptor PDBQT → líneas PDB estrictas en cadena A."""
    lines  = []
    serial = 1
    with open(pdbqt_path, encoding='ascii', errors='ignore') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                atom_name = line[12:16].strip()
                resname   = line[17:20].strip() or 'UNK'
                chain     = 'A'
                resseq    = int(line[22:26].strip() or 1)
                x         = float(line[30:38])
                y         = float(line[38:46])
                z         = float(line[46:54])
                ad4       = line[77:79].strip() if len(line) > 77 else atom_name[:1]
                element   = elem_from_ad4(ad4)

                lines.append(make_pdb_line('ATOM', serial, atom_name, resname, chain, resseq, x, y, z, element) + '\n')
                serial += 1
            except (ValueError, IndexError):
                continue
    return lines


def read_ligand_model1(pdbqt_path):
    """Extrae MODEL 1 del PDBQT Vina → líneas HETATM estrictas en cadena L."""
    lines      = []
    in_model1  = False
    base_serial = 10000

    with open(pdbqt_path, encoding='ascii', errors='ignore') as f:
        for line in f:
            if line.startswith('MODEL 1'):
                in_model1 = True
                continue
            if line.startswith('ENDMDL') and in_model1:
                break
            if not in_model1 or not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                atom_name = line[12:16].strip()
                x         = float(line[30:38])
                y         = float(line[38:46])
                z         = float(line[46:54])
                ad4       = line[77:79].strip() if len(line) > 77 else atom_name[:1]
                element   = elem_from_ad4(ad4)
                serial    = base_serial + len(lines) + 1

                lines.append(make_pdb_line('HETATM', serial, atom_name, 'LIG', 'L', 1, x, y, z, element) + '\n')
            except (ValueError, IndexError):
                continue
    return lines


def merge_complex(receptor_pdbqt, ligand_out_pdbqt, out_pdb):
    rec_lines = read_receptor_atoms(receptor_pdbqt)
    lig_lines = read_ligand_model1(ligand_out_pdbqt)

    if not rec_lines:
        print(f'  ✗ Receptor vacío: {receptor_pdbqt}')
        return False
    if not lig_lines:
        print(f'  ✗ Sin poses en: {ligand_out_pdbqt}')
        return False

    with open(out_pdb, 'w', encoding='ascii') as f:
        f.write('REMARK GEN-MAMA-MX Complex\n')
        f.write(f'REMARK Receptor: {Path(receptor_pdbqt).name}\n')
        f.write(f'REMARK Ligand:   {Path(ligand_out_pdbqt).name}\n')
        f.write('REMARK Chain A = Protein | Chain L = Ligand (Vina pose 1)\n')
        f.writelines(rec_lines)
        f.write('TER\n')
        f.writelines(lig_lines)
        f.write('END\n')

    size_kb = Path(out_pdb).stat().st_size / 1024
    print(f'  ✓ {Path(out_pdb).name} ({size_kb:.0f} KB, {len(rec_lines)} rec + {len(lig_lines)} lig)')
    return True


def main():
    print('=' * 55)
    print('MERGE POSES — PDB formato estricto para PDBe Molstar')
    print('=' * 55)

    with open('data/docking/vina_affinities.json') as f:
        vina_data = json.load(f)

    ok = skip = fail = 0

    for variant_id, drug_scores in vina_data.items():
        receptor_pdbqt = RECEPTORS_DIR / f'{variant_id}.pdbqt'
        if not receptor_pdbqt.exists() or receptor_pdbqt.stat().st_size < 1000:
            print(f'\n[{variant_id}] ⚠ Sin receptor válido')
            skip += len(drug_scores)
            continue

        print(f'\n[{variant_id}]')
        for drug_name, affinity in drug_scores.items():
            if affinity is None:
                skip += 1
                continue

            drug_safe   = drug_name.replace(' ', '_')
            ligand_out  = RESULTS_DIR / f'{variant_id}_{drug_safe}_out.pdbqt'
            complex_out = COMPLEX_DIR / f'{variant_id}_{drug_safe}_complex.pdb'

            if not ligand_out.exists():
                print(f'  ⚠ {drug_name}: sin poses')
                skip += 1
                continue

            success = merge_complex(str(receptor_pdbqt), str(ligand_out), str(complex_out))
            if success: ok += 1
            else: fail += 1

    print(f'\n{"="*55}')
    print(f'  ✓ Generados: {ok}  ⚠ Omitidos: {skip}  ✗ Fallidos: {fail}')


if __name__ == '__main__':
    main()
