from __future__ import annotations

import argparse
from pathlib import Path


def read_gro(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Archivo GRO invalido: {path}")
    title = lines[0]
    natoms = int(lines[1].strip())
    atoms = lines[2 : 2 + natoms]
    box = lines[2 + natoms]
    return title, natoms, atoms, box


def write_gro(path: Path, title: str, atoms: list[str], box: str):
    payload = [title, f"{len(atoms):5d}", *atoms, box]
    path.write_text("\n".join(payload) + "\n", encoding="utf-8")


def merge_gro(protein_gro: Path, ligand_gro: Path, output_gro: Path):
    p_title, _, p_atoms, p_box = read_gro(protein_gro)
    l_title, _, l_atoms, l_box = read_gro(ligand_gro)
    box = p_box if p_box.strip() else l_box
    write_gro(output_gro, f"{p_title} + ligand", p_atoms + l_atoms, box)


def update_topology(
    topol_path: Path,
    ligand_itp: Path,
    ligand_posre: Path,
    ligand_name: str,
    output_topol: Path | None = None,
):
    target = output_topol or topol_path
    text = topol_path.read_text(encoding="utf-8")

    include_block = (
        f'\n; Include ligand topology\n#include "{ligand_itp.name}"\n'
        f'#ifdef POSRES_LIG\n#include "{ligand_posre.name}"\n#endif\n'
    )

    if ligand_itp.name not in text:
        if "[ system ]" in text:
            text = text.replace("[ system ]", include_block + "\n[ system ]", 1)
        else:
            text += include_block

    lines = text.splitlines()
    out_lines = []
    in_molecules = False
    ligand_added = False

    for line in lines:
        out_lines.append(line)
        stripped = line.strip()
        if stripped.startswith("[ molecules ]"):
            in_molecules = True
            continue
        if in_molecules:
            if stripped and not stripped.startswith(";"):
                pass
            if stripped.lower().startswith(ligand_name.lower()):
                ligand_added = True

    if "[ molecules ]" not in text:
        out_lines.extend(["", "[ molecules ]", "; Compound        nmols", f"{ligand_name:<16} 1"])
    elif not ligand_added:
        out_lines.append(f"{ligand_name:<16} 1")

    target.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def copy_if_needed(src: Path, dst: Path):
    if not dst.exists():
        dst.write_bytes(src.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a protein-ligand GROMACS system from ACPYPE outputs")
    parser.add_argument("--run_dir", required=True)
    parser.add_argument("--protein_gro", default="receptor_processed.gro")
    parser.add_argument("--protein_top", default="topol.top")
    parser.add_argument("--ligand_gro", required=True)
    parser.add_argument("--ligand_itp", required=True)
    parser.add_argument("--ligand_posre", required=True)
    parser.add_argument("--ligand_name", required=True)
    parser.add_argument("--output_gro", default="complex_init.gro")
    parser.add_argument("--output_top", default="topol.top")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    protein_gro = run_dir / args.protein_gro
    protein_top = run_dir / args.protein_top
    ligand_gro = run_dir / args.ligand_gro
    ligand_itp = run_dir / args.ligand_itp
    ligand_posre = run_dir / args.ligand_posre
    output_gro = run_dir / args.output_gro
    output_top = run_dir / args.output_top

    if not protein_gro.exists():
        raise SystemExit(f"No existe {protein_gro}")
    if not protein_top.exists():
        raise SystemExit(f"No existe {protein_top}")
    if not ligand_gro.exists():
        raise SystemExit(f"No existe {ligand_gro}")
    if not ligand_itp.exists():
        raise SystemExit(f"No existe {ligand_itp}")
    if not ligand_posre.exists():
        raise SystemExit(f"No existe {ligand_posre}")

    copy_if_needed(ligand_itp, run_dir / ligand_itp.name)
    copy_if_needed(ligand_posre, run_dir / ligand_posre.name)

    merge_gro(protein_gro, ligand_gro, output_gro)
    update_topology(protein_top, run_dir / ligand_itp.name, run_dir / ligand_posre.name, args.ligand_name, output_top)

    print(f"complex_gro={output_gro}")
    print(f"topology={output_top}")


if __name__ == "__main__":
    main()
