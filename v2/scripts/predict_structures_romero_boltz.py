#!/usr/bin/env python3
"""
GEMAMX v2 - GPU structure prediction for Romero cohort using Boltz in Colab.

Inputs:
- v2/data/fasta_sequences.json

Outputs:
- v2/data/structures/*.cif
- v2/data/manifest_structures.json
- v2/data/structure_jobs.json
- v2/data/boltz_runs/* (raw Boltz outputs)

Notes:
- Uses full sequence for proteins <= 390 aa.
- Uses mutation-centered fragments (+/- 180 aa) for larger proteins.
- Designed to run in Google Colab with GPU and boltz installed.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FASTA_INDEX = ROOT / "v2" / "data" / "fasta_sequences.json"
STRUCTURES_DIR = ROOT / "v2" / "data" / "structures"
RUNS_DIR = ROOT / "v2" / "data" / "boltz_runs"
META_DIR = STRUCTURES_DIR / "meta"
MANIFEST_PATH = ROOT / "v2" / "data" / "manifest_structures.json"
JOBS_PATH = ROOT / "v2" / "data" / "structure_jobs.json"

STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
RUNS_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)

BOLTZ_LIMIT = 390
FRAGMENT_WINDOW = 180


def parse_fasta_sequence(path: Path) -> str:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return "".join(line for line in lines if not line.startswith(">"))


def get_fragment(seq: str, mut_pos: int, window: int = FRAGMENT_WINDOW):
    start = max(0, mut_pos - window - 1)
    end = min(len(seq), mut_pos + window)
    return seq[start:end], start + 1, end


def write_job_fasta(name: str, seq: str, run_dir: Path) -> Path:
    fasta_path = run_dir / f"{name}.fasta"
    fasta_path.write_text(f">{name}\n{seq}\n", encoding="utf-8")
    return fasta_path


def find_prediction_file(run_dir: Path, name: str) -> Path | None:
    exact_candidates = [
        run_dir / name / "predictions" / name / f"{name}_model_0.cif",
        run_dir / "predictions" / name / f"{name}_model_0.cif",
        run_dir / f"{name}_model_0.cif",
    ]
    for candidate in exact_candidates:
        if candidate.exists() and candidate.stat().st_size > 1000:
            return candidate

    files = sorted(
        [p for p in run_dir.rglob("*.cif") if p.stat().st_size > 1000],
        key=lambda p: p.stat().st_size,
        reverse=True,
    )
    return files[0] if files else None


def find_confidence_file(run_dir: Path, name: str) -> Path | None:
    exact_candidates = [
        run_dir / name / "predictions" / name / f"confidence_{name}_model_0.json",
        run_dir / "predictions" / name / f"confidence_{name}_model_0.json",
    ]
    for candidate in exact_candidates:
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate

    files = sorted(run_dir.rglob("confidence_*.json"))
    return files[0] if files else None


def run_boltz(name: str, seq: str, force: bool = False) -> tuple[bool, str]:
    run_dir = RUNS_DIR / name
    if force and run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = write_job_fasta(name, seq, run_dir)
    cmd = ["boltz", "predict", str(fasta_path), "--out_dir", str(run_dir), "--name", name]
    print("  Running:", " ".join(cmd))
    completed = subprocess.run(cmd, cwd=ROOT)
    if completed.returncode != 0:
        return False, f"boltz exited with code {completed.returncode}"

    pred = find_prediction_file(run_dir, name)
    if not pred:
        return False, "prediction file not found"

    out_path = STRUCTURES_DIR / f"{name}{pred.suffix.lower()}"
    shutil.copy2(pred, out_path)

    conf = find_confidence_file(run_dir, name)
    if conf:
        shutil.copy2(conf, META_DIR / f"{name}_confidence.json")

    return True, str(out_path.relative_to(ROOT)).replace('\\', '/')


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

        short_protein = len(wt_seq) <= BOLTZ_LIMIT and len(mut_seq) <= BOLTZ_LIMIT
        if short_protein:
            wt_input_seq = wt_seq
            mut_input_seq = mut_seq
            wt_key = f"{gene}_WT"
            wt_fragment = {"start": 1, "end": len(wt_seq), "full_length": len(wt_seq), "mode": "full"}
            mut_fragment = {"start": 1, "end": len(mut_seq), "full_length": len(mut_seq), "mode": "full"}
        else:
            wt_input_seq, wt_start, wt_end = get_fragment(wt_seq, mut_pos)
            mut_input_seq, mut_start, mut_end = get_fragment(mut_seq, mut_pos)
            wt_key = f"WT__{variant_id}"
            wt_fragment = {"start": wt_start, "end": wt_end, "full_length": len(wt_seq), "mode": "fragment"}
            mut_fragment = {"start": mut_start, "end": mut_end, "full_length": len(mut_seq), "mode": "fragment"}

        jobs.append({
            "variant_id": variant_id,
            "gene": gene,
            "protein_change": change,
            "mutation_pos": mut_pos,
            "wt_key": wt_key,
            "wt_sequence_length": len(wt_input_seq),
            "wt_fragment": wt_fragment,
            "mut_key": variant_id,
            "mut_sequence_length": len(mut_input_seq),
            "mut_fragment": mut_fragment,
            "wt_input_seq": wt_input_seq,
            "mut_input_seq": mut_input_seq,
        })
    return jobs


def load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        try:
            return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def existing_structure_path(key: str) -> str | None:
    for suffix in (".cif", ".pdb"):
        candidate = STRUCTURES_DIR / f"{key}{suffix}"
        if candidate.exists() and candidate.stat().st_size > 1000:
            return str(candidate.relative_to(ROOT)).replace('\\', '/')
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=None, help="Process only the first N variants")
    parser.add_argument("--force", action="store_true", help="Regenerate structures even if they already exist")
    parser.add_argument("--only-missing", action="store_true", help="Skip WT/MUT that already exist in v2/data/structures")
    args = parser.parse_args()

    try:
        subprocess.run(["boltz", "--version"], check=True, capture_output=True)
    except Exception as exc:
        raise SystemExit("boltz is not available in PATH. Install it in Colab first.") from exc

    jobs = build_jobs(limit=args.limit)
    manifest = load_manifest()
    completed = []
    generated_wt = set()

    for idx, job in enumerate(jobs, start=1):
        variant_id = job["variant_id"]
        print(f"[{idx}/{len(jobs)}] {variant_id}")

        wt_path = existing_structure_path(job["wt_key"])
        if args.force:
            wt_path = None
        if job["wt_key"] not in generated_wt:
            if args.only_missing and wt_path:
                print(f"  WT already exists -> {job['wt_key']}")
            elif wt_path and not args.force:
                print(f"  WT cached -> {job['wt_key']}")
            else:
                print(f"  WT -> {job['wt_key']} ({job['wt_fragment']['mode']}, {job['wt_sequence_length']} aa)")
                ok, result = run_boltz(job["wt_key"], job["wt_input_seq"], force=args.force)
                if not ok:
                    print(f"    WT failed: {result}")
                else:
                    wt_path = result
                    print(f"    WT saved: {result}")
            generated_wt.add(job["wt_key"])

        mut_path = existing_structure_path(job["mut_key"])
        if args.force:
            mut_path = None
        if args.only_missing and mut_path:
            print(f"  MUT already exists -> {job['mut_key']}")
        elif mut_path and not args.force:
            print(f"  MUT cached -> {job['mut_key']}")
        else:
            print(f"  MUT -> {job['mut_key']} ({job['mut_fragment']['mode']}, {job['mut_sequence_length']} aa)")
            ok, result = run_boltz(job["mut_key"], job["mut_input_seq"], force=args.force)
            if not ok:
                print(f"    MUT failed: {result}")
            else:
                mut_path = result
                print(f"    MUT saved: {result}")

        wt_path = wt_path or existing_structure_path(job["wt_key"])
        mut_path = mut_path or existing_structure_path(job["mut_key"])

        if wt_path:
            manifest[job["wt_key"]] = wt_path
        if mut_path:
            manifest[job["mut_key"]] = mut_path
        if wt_path and mut_path:
            completed.append({k: v for k, v in job.items() if not k.endswith("_input_seq")})

    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    JOBS_PATH.write_text(json.dumps(completed, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\nGEMAMX v2 Romero Boltz run complete")
    print(f"  manifest entries: {len(manifest)}")
    print(f"  completed variant pairs: {len(completed)}")
    print(f"  manifest: {MANIFEST_PATH}")
    print(f"  jobs: {JOBS_PATH}")


if __name__ == "__main__":
    main()
