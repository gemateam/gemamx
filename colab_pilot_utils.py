from __future__ import annotations

import json
import shutil
from pathlib import Path


EM_MDP = """integrator  = steep
nsteps      = 5000
emtol       = 1000.0
emstep      = 0.01
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""

NVT_MDP = """title       = NVT equilibration
define      = -DPOSRES
integrator  = md
nsteps      = 50000
dt          = 0.002
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300
pcoupl      = no
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""

NPT_MDP = """title       = NPT equilibration
define      = -DPOSRES
integrator  = md
nsteps      = 50000
dt          = 0.002
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
"""

MD_MDP = """title       = Production MD
integrator  = md
nsteps      = 10000000
dt          = 0.002
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 300 300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
constraints = h-bonds
constraint_algorithm = lincs
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
nstxout-compressed = 5000
nstenergy   = 1000
nstlog      = 1000
"""


def load_manifest(path: str | Path) -> dict:
    path = Path(path)
    return json.loads(path.read_text(encoding="utf-8"))


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def copy_inputs(pair: dict, run_dir: str | Path) -> dict:
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    copied = {}
    for key, src in pair.get("input_paths", {}).items():
        src_path = Path(src)
        if src_path.exists():
            dst = run_dir / src_path.name
            shutil.copy2(src_path, dst)
            copied[key] = str(dst)
        else:
            copied[key] = None
    return copied


def write_mdp_bundle(run_dir: str | Path, md_time_ns: float = 20.0) -> None:
    run_dir = Path(run_dir)
    md_steps = int((md_time_ns * 1000) / 0.002)
    md_content = MD_MDP.replace("nsteps      = 10000000", f"nsteps      = {md_steps}")
    write_text(run_dir / "em.mdp", EM_MDP)
    write_text(run_dir / "nvt.mdp", NVT_MDP)
    write_text(run_dir / "npt.mdp", NPT_MDP)
    write_text(run_dir / "md.mdp", md_content)


def write_run_script(pair: dict, run_dir: str | Path, md_time_ns: float = 20.0) -> Path:
    run_dir = Path(run_dir)
    safe_drug = pair["drug_name"].replace(" ", "_")
    script = f"""#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

echo "[1/8] Convertir pose y ligando"
obabel "{pair['variant_id']}_{safe_drug}_out.pdbqt" -O ligand_pose.pdb || true
obabel "{safe_drug}.sdf" -O ligand_input.mol2 || true

echo "[2/8] Preparar topologia del ligando"
echo "TODO: generar topologia del ligando con ACPYPE o CGenFF"
echo "Ejemplo ACPYPE:"
echo "  acpype -i ligand_input.mol2 -c bcc -n 0"

echo "[3/8] Preparar receptor"
gmx pdb2gmx -f "{pair['variant_id']}_clean.pdb" -o receptor_processed.gro -water tip3p -ff charmm36-jul2022

echo "[4/8] Combinar proteina + ligando"
python /content/drive/MyDrive/GEMAMX/gromacs/prepare_protein_ligand_system.py \\
  --run_dir "$(pwd)" \\
  --protein_gro receptor_processed.gro \\
  --protein_top topol.top \\
  --ligand_gro {safe_drug.lower()}.acpype/{safe_drug.lower()}_GMX.gro \\
  --ligand_itp {safe_drug.lower()}.acpype/{safe_drug.lower()}_GMX.itp \\
  --ligand_posre {safe_drug.lower()}.acpype/posre_{safe_drug.lower()}.itp \\
  --ligand_name {safe_drug.lower()} \\
  --output_gro complex_init.gro \\
  --output_top topol.top

echo "[5/9] Construir caja"
gmx editconf -f complex_init.gro -o boxed.gro -c -d 1.0 -bt cubic

echo "[6/9] Solvatar"
gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top

echo "[7/9] Iones"
gmx grompp -f em.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 5
printf "SOL\\n" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

echo "[8/9] Minimización + equilibración"
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 5
gmx mdrun -deffnm em
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 5
gmx mdrun -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 5
gmx mdrun -deffnm npt

echo "[9/9] Producción"
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 5
gmx mdrun -deffnm md

echo "Producción completada para {pair['tag']}"
echo "TODO post-proceso:"
echo "  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg"
echo "  gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg"
echo "  gmx hbond -s md.tpr -f md.xtc -num hbonds.xvg"
echo "  MM-GBSA con gmx_MMPBSA si el runtime lo soporta"
"""
    script_path = run_dir / "run_gromacs.sh"
    write_text(script_path, script)
    return script_path


def write_pair_metadata(pair: dict, run_dir: str | Path, copied_inputs: dict, md_time_ns: float = 20.0) -> Path:
    run_dir = Path(run_dir)
    payload = {
        "tag": pair["tag"],
        "variant_id": pair["variant_id"],
        "drug_name": pair["drug_name"],
        "vina_score_kcalmol": pair["vina_score_kcalmol"],
        "run_dir": str(run_dir),
        "md_time_ns": md_time_ns,
        "copied_inputs": copied_inputs,
        "status": "scaffolded"
    }
    path = run_dir / "pair_metadata.json"
    write_text(path, json.dumps(payload, indent=2, ensure_ascii=False))
    return path


def scaffold_pair(pair: dict, base_drive_root: str | Path, md_time_ns: float = 20.0) -> dict:
    base_drive_root = Path(base_drive_root)
    run_dir = base_drive_root / pair["recommended_colab_run_dir"]
    copied_inputs = copy_inputs(pair, run_dir)
    write_mdp_bundle(run_dir, md_time_ns=md_time_ns)
    script_path = write_run_script(pair, run_dir, md_time_ns=md_time_ns)
    metadata_path = write_pair_metadata(pair, run_dir, copied_inputs, md_time_ns=md_time_ns)
    return {
        "tag": pair["tag"],
        "run_dir": str(run_dir),
        "script_path": str(script_path),
        "metadata_path": str(metadata_path)
    }


def aggregate_pair_metrics(results: list[dict], output_path: str | Path, parameters: dict | None = None, pipeline_version: str = "GEMAMX-GROMACS-v2-pilot") -> Path:
    output_path = Path(output_path)
    payload = {
        "generated_at": None,
        "pipeline_version": pipeline_version,
        "parameters": parameters or {},
        "summary": {
            "pairs_total": len(results),
            "pairs_completed": sum(1 for row in results if row.get("status") == "completed"),
            "pairs_failed": sum(1 for row in results if row.get("status") in ("failed", "error"))
        },
        "pairs": results
    }
    write_text(output_path, json.dumps(payload, indent=2, ensure_ascii=False))
    return output_path
