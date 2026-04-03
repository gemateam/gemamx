#!/usr/bin/env python3
"""
GEMAMX v2 - Structure generation for Romero cohort using ESMFold API.

Input:
- v2/data/fasta_sequences.json

Outputs:
- v2/data/structures/*.pdb
- v2/data/manifest_structures.json
- v2/data/structure_jobs.json

Notes:
- For proteins <= 390 aa, uses full WT and full mutant sequence.
- For proteins > 390 aa, uses a mutation-centered fragment (+/- 180 aa).
- Mutant structures are always variant-specific.
- WT structures are reused per gene for short proteins and are variant-specific for long proteins.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
FASTA_INDEX = ROOT / "v2" / "data" / "fasta_sequences.json"
STRUCTURES_DIR = ROOT / "v2" / "data" / "structures"
STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
MANIFEST_PATH = ROOT / "v2" / "data" / "manifest_structures.json"
JOBS_PATH = ROOT / "v2" / "data" / "structure_jobs.json"

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_LIMIT = 390
FRAGMENT_WINDOW = 180
REQUEST_TIMEOUT = 240
SLEEP_SECONDS = 1.0


def parse_fasta_sequence(path: Path) -> str:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return "".join(line for line in lines if not line.startswith(">"))


def get_fragment(seq: str, mut_pos: int, window: int = FRAGMENT_WINDOW):
    start = max(0, mut_pos - window - 1)
    end = min(len(seq), mut_pos + window)
    return seq[start:end], start + 1, end


def predict_structure(seq: str, out_path: Path) -> bool:
    response = requests.post(
        ESMFOLD_URL,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data=seq,
        timeout=REQUEST_TIMEOUT,
    )
    if not response.ok or "ATOM" not in response.text:
        return False
    out_path.write_text(response.text, encoding="utf-8")
    return True


def build_jobs(limit: int | None = None):
    entries = json.loads(FASTA_INDEX.read_text(encoding="utf-8"))
    if limit is not None:
        entries = entries[:limit]

    jobs = []
    for entry in entries:
        gene = entry["gene"]
        change = entry["protein_change"]
        variant_id = f"{gene}_{change}"
        wt_fasta = ROOT / entry["wt_fasta"]
        mut_fasta = ROOT / entry["mut_fasta"]
        mut_pos = int(entry["mutation_pos"])

        wt_seq = parse_fasta_sequence(wt_fasta)
        mut_seq = parse_fasta_sequence(mut_fasta)

        short_protein = len(wt_seq) <= ESMFOLD_LIMIT and len(mut_seq) <= ESMFOLD_LIMIT
        if short_protein:
            wt_input_seq = wt_seq
            mut_input_seq = mut_seq
            wt_key = f"{gene}_WT"
            wt_out = STRUCTURES_DIR / f"{gene}_WT.pdb"
            wt_fragment = {"start": 1, "end": len(wt_seq), "full_length": len(wt_seq), "mode": "full"}
            mut_fragment = {"start": 1, "end": len(mut_seq), "full_length": len(mut_seq), "mode": "full"}
        else:
            wt_input_seq, wt_start, wt_end = get_fragment(wt_seq, mut_pos)
            mut_input_seq, mut_start, mut_end = get_fragment(mut_seq, mut_pos)
            wt_key = f"WT__{variant_id}"
            wt_out = STRUCTURES_DIR / f"WT__{variant_id}.pdb"
            wt_fragment = {"start": wt_start, "end": wt_end, "full_length": len(wt_seq), "mode": "fragment"}
            mut_fragment = {"start": mut_start, "end": mut_end, "full_length": len(mut_seq), "mode": "fragment"}

        mut_key = variant_id
        mut_out = STRUCTURES_DIR / f"{variant_id}.pdb"

        jobs.append({
            "variant_id": variant_id,
            "gene": gene,
            "protein_change": change,
            "mutation_pos": mut_pos,
            "wt_key": wt_key,
            "wt_out": str(wt_out.relative_to(ROOT)).replace('\\', '/'),
            "wt_sequence_length": len(wt_input_seq),
            "wt_fragment": wt_fragment,
            "mut_key": mut_key,
            "mut_out": str(mut_out.relative_to(ROOT)).replace('\\', '/'),
            "mut_sequence_length": len(mut_input_seq),
            "mut_fragment": mut_fragment,
            "wt_input_seq": wt_input_seq,
            "mut_input_seq": mut_input_seq,
        })
    return jobs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="Process only the first N variants")
    parser.add_argument("--force", action="store_true", help="Regenerate structures even if files already exist")
    args = parser.parse_args()

    jobs = build_jobs(limit=args.limit)
    manifest = {}
    completed = []

    generated_wt = set()

    for idx, job in enumerate(jobs, start=1):
        variant_id = job["variant_id"]
        wt_path = ROOT / job["wt_out"]
        mut_path = ROOT / job["mut_out"]
        print(f"[{idx}/{len(jobs)}] {variant_id}")

        wt_generated = wt_path.exists() and wt_path.stat().st_size > 1000 and not args.force
        if not wt_generated and job["wt_key"] not in generated_wt:
            print(f"  WT -> {job['wt_key']} ({job['wt_fragment']['mode']}, {job['wt_sequence_length']} aa)")
            ok = predict_structure(job["wt_input_seq"], wt_path)
            if not ok:
                print("    WT prediction failed")
            else:
                generated_wt.add(job["wt_key"])
                time.sleep(SLEEP_SECONDS)
        else:
            generated_wt.add(job["wt_key"])

        mut_generated = mut_path.exists() and mut_path.stat().st_size > 1000 and not args.force
        if not mut_generated:
            print(f"  MUT -> {variant_id} ({job['mut_fragment']['mode']}, {job['mut_sequence_length']} aa)")
            ok = predict_structure(job["mut_input_seq"], mut_path)
            if not ok:
                print("    MUT prediction failed")
            else:
                time.sleep(SLEEP_SECONDS)

        wt_ok = wt_path.exists() and wt_path.stat().st_size > 1000
        mut_ok = mut_path.exists() and mut_path.stat().st_size > 1000

        if wt_ok:
            manifest[job["wt_key"]] = job["wt_out"]
        if mut_ok:
            manifest[job["mut_key"]] = job["mut_out"]
        if wt_ok and mut_ok:
            completed.append({k: v for k, v in job.items() if not k.endswith('_input_seq')})

    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    JOBS_PATH.write_text(json.dumps(completed, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\nGEMAMX v2 Romero structure generation complete")
    print(f"  manifest entries: {len(manifest)}")
    print(f"  completed variant pairs: {len(completed)}")
    print(f"  manifest: {MANIFEST_PATH}")
    print(f"  jobs: {JOBS_PATH}")


if __name__ == "__main__":
    main()
