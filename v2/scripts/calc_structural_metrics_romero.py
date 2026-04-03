#!/usr/bin/env python3
"""
GEMAMX v2 - Structural metrics for Romero cohort.

Inputs:
- v2/data/structure_jobs.json
- v2/data/manifest_structures.json

Output:
- v2/data/structural_metrics.json

Notes:
- Uses PDB B-factor as pLDDT proxy for ESMFold outputs.
- Computes CA RMSD globally and in a local mutation window.
- Works on fragment-based WT/MUT jobs generated in v2.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
JOBS_PATH = ROOT / "v2" / "data" / "structure_jobs.json"
OUT_PATH = ROOT / "v2" / "data" / "structural_metrics.json"

LOCAL_WINDOW = 12


def parse_pdb_ca(path: Path):
    atoms = []
    with path.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            resseq = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            b = float(line[60:66]) if line[60:66].strip() else 0.0
            atoms.append({"resseq": resseq, "coord": (x, y, z), "plddt": b})
    return atoms


def rmsd(coords_a, coords_b):
    n = min(len(coords_a), len(coords_b))
    if n == 0:
        return None
    s = 0.0
    for i in range(n):
        ax, ay, az = coords_a[i]
        bx, by, bz = coords_b[i]
        s += (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
    return math.sqrt(s / n)


def mean(values):
    vals = [v for v in values if v is not None]
    return (sum(vals) / len(vals)) if vals else None


def scale_plddt(value):
    if value is None:
        return None
    return value * 100.0 if value <= 1.5 else value


def main():
    jobs = json.loads(JOBS_PATH.read_text(encoding="utf-8"))
    metrics = {}

    for job in jobs:
        wt_path = ROOT / job["wt_out"]
        mut_path = ROOT / job["mut_out"]
        if not wt_path.exists() or not mut_path.exists():
            continue

        wt_atoms = parse_pdb_ca(wt_path)
        mut_atoms = parse_pdb_ca(mut_path)
        n = min(len(wt_atoms), len(mut_atoms))
        if n == 0:
            continue

        wt_atoms = wt_atoms[:n]
        mut_atoms = mut_atoms[:n]
        wt_coords = [a["coord"] for a in wt_atoms]
        mut_coords = [a["coord"] for a in mut_atoms]

        mut_pos = int(job["mutation_pos"])
        frag_start = int(job["mut_fragment"]["start"])
        local_start = max(frag_start, mut_pos - LOCAL_WINDOW)
        local_end = min(int(job["mut_fragment"]["end"]), mut_pos + LOCAL_WINDOW)

        local_wt = []
        local_mut = []
        for wa, ma in zip(wt_atoms, mut_atoms):
            resseq = wa["resseq"]
            genomic_like = frag_start + (resseq - 1)
            if local_start <= genomic_like <= local_end:
                local_wt.append(wa["coord"])
                local_mut.append(ma["coord"])

        wt_plddt = scale_plddt(mean([a["plddt"] for a in wt_atoms]))
        mut_plddt = scale_plddt(mean([a["plddt"] for a in mut_atoms]))
        rmsd_global = rmsd(wt_coords, mut_coords)
        rmsd_active = rmsd(local_wt, local_mut)

        region_label = f"Mutation window +/-{LOCAL_WINDOW} aa"
        metrics[job["variant_id"]] = {
            "wt_plddt": round(wt_plddt, 1) if wt_plddt is not None else "--",
            "mut_plddt": round(mut_plddt, 1) if mut_plddt is not None else "--",
            "rmsd_global": round(rmsd_global, 2) if rmsd_global is not None else "--",
            "rmsd_active": round(rmsd_active, 2) if rmsd_active is not None else "--",
            "rmsd_regions": [
                {
                    "region": region_label,
                    "value": round(rmsd_active, 2) if rmsd_active is not None else 0,
                    "color": "#f59e0b" if (rmsd_active or 0) < 2 else "#f87171",
                }
            ],
            "structure_mode": job["mut_fragment"]["mode"],
            "fragment_start": job["mut_fragment"]["start"],
            "fragment_end": job["mut_fragment"]["end"],
            "mutation_pos": mut_pos,
        }

    OUT_PATH.write_text(json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8")
    print("GEMAMX v2 Romero structural metrics complete")
    print(f"  variants with metrics: {len(metrics)}")
    print(f"  output: {OUT_PATH}")
    for vid, m in list(metrics.items())[:10]:
        print(f"    {vid}: global={m['rmsd_global']} active={m['rmsd_active']} pLDDT={m['wt_plddt']}->{m['mut_plddt']}")


if __name__ == "__main__":
    main()
