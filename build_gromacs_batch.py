from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCKING_DIR = ROOT / "data" / "docking"
VINA_JSON = DOCKING_DIR / "vina_affinities.json"
RESULTS_DIR = ROOT / "gromacs" / "results"
MANIFEST_JSON = RESULTS_DIR / "gromacs_batch_manifest.json"


def load_json(path: Path):
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def build_pair_record(variant_id: str, drug_name: str, score: float) -> dict:
    safe_drug = drug_name.replace(" ", "_")
    tag = f"{variant_id}__{safe_drug}"
    return {
        "tag": tag,
        "variant_id": variant_id,
        "drug_name": drug_name,
        "vina_score_kcalmol": score,
        "input_paths": {
            "vina_pose_pdbqt": str(DOCKING_DIR / "results" / f"{variant_id}_{safe_drug}_out.pdbqt"),
            "vina_log_txt": str(DOCKING_DIR / "results" / f"{variant_id}_{safe_drug}_log.txt"),
            "complex_pdb": str(DOCKING_DIR / "complex" / f"{variant_id}_{safe_drug}_complex.pdb"),
            "receptor_pdbqt": str(DOCKING_DIR / "receptors" / f"{variant_id}.pdbqt"),
            "receptor_clean_pdb": str(DOCKING_DIR / "receptors" / f"{variant_id}_clean.pdb"),
            "ligand_pdbqt": str(DOCKING_DIR / "ligands" / f"{safe_drug}.pdbqt"),
            "ligand_sdf": str(DOCKING_DIR / "ligands" / f"{safe_drug}.sdf")
        },
        "recommended_colab_run_dir": f"gromacs/runs/{tag}"
    }


def existing_inputs(pair: dict) -> dict:
    availability = {}
    for key, path_str in pair["input_paths"].items():
        availability[key] = Path(path_str).exists()
    return availability


def build_manifest(top_n_per_variant: int) -> dict:
    vina = load_json(VINA_JSON)
    pairs = []
    per_variant_counts = {}

    for variant_id, drug_scores in sorted(vina.items()):
        ranked = sorted(
            ((drug, float(score)) for drug, score in drug_scores.items() if score is not None),
            key=lambda item: item[1]
        )[:top_n_per_variant]
        per_variant_counts[variant_id] = len(ranked)

        for drug_name, score in ranked:
            pair = build_pair_record(variant_id, drug_name, score)
            pair["available_inputs"] = existing_inputs(pair)
            pairs.append(pair)

    return {
        "generated_at": datetime.now().isoformat(),
        "pipeline_version": "GEMAMX-GROMACS-v2",
        "selection_policy": {
            "source": str(VINA_JSON),
            "top_n_per_variant": top_n_per_variant,
            "sort_order": "ascending vina score"
        },
        "summary": {
            "variants_considered": len(vina),
            "pairs_selected": len(pairs),
            "per_variant_counts": per_variant_counts
        },
        "pairs": pairs
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build GROMACS batch manifest from Vina results")
    parser.add_argument("--top_n_per_variant", type=int, default=5)
    parser.add_argument("--variants", nargs="*", default=None)
    parser.add_argument("--drugs", nargs="*", default=None)
    parser.add_argument("--output", default=str(MANIFEST_JSON))
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    manifest = build_manifest(args.top_n_per_variant)
    if args.variants:
        allowed_variants = set(args.variants)
        manifest["pairs"] = [row for row in manifest["pairs"] if row["variant_id"] in allowed_variants]
    if args.drugs:
        allowed_drugs = set(args.drugs)
        manifest["pairs"] = [row for row in manifest["pairs"] if row["drug_name"] in allowed_drugs]

    per_variant_counts = {}
    for row in manifest["pairs"]:
        per_variant_counts[row["variant_id"]] = per_variant_counts.get(row["variant_id"], 0) + 1
    manifest["summary"]["pairs_selected"] = len(manifest["pairs"])
    manifest["summary"]["variants_considered"] = len(per_variant_counts) if manifest["pairs"] else 0
    manifest["summary"]["per_variant_counts"] = per_variant_counts
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    print(json.dumps({
        "generated_at": manifest["generated_at"],
        "pairs_selected": manifest["summary"]["pairs_selected"],
        "variants_considered": manifest["summary"]["variants_considered"],
        "output": str(output_path)
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
