from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GROMACS_DIR = ROOT / "gromacs"
RESULTS_DIR = GROMACS_DIR / "results"
DEFAULT_INPUT = RESULTS_DIR / "gromacs_metrics.json"
LOCAL_COPY = GROMACS_DIR / "gromacs_metrics.json"
STRUCTURAL_JSON = ROOT / "data" / "structural_metrics.json"


def load_json(path: Path):
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def save_json(data, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def normalize_pairs(payload: dict) -> list[dict]:
    raw_pairs = payload.get("pairs", payload)
    pairs = []
    if isinstance(raw_pairs, dict):
        for tag, row in raw_pairs.items():
            item = dict(row or {})
            item.setdefault("tag", tag)
            pairs.append(item)
    elif isinstance(raw_pairs, list):
        for row in raw_pairs:
            if isinstance(row, dict):
                item = dict(row)
                item.setdefault("tag", f'{item.get("variant_id","unknown")}__{item.get("drug_name","unknown")}')
                pairs.append(item)
    return pairs


def pick_best_pair(rows: list[dict]) -> dict | None:
    valid = [row for row in rows if row.get("mmgbsa_delta_g") is not None]
    if valid:
        return min(valid, key=lambda row: float(row["mmgbsa_delta_g"]))
    completed = [row for row in rows if row.get("status") == "completed"]
    if completed:
        return min(completed, key=lambda row: float(row.get("vina_score_kcalmol", 999)))
    return None


def ingest_results(input_path: Path) -> dict:
    payload = load_json(input_path)
    pairs = normalize_pairs(payload)
    structural = load_json(STRUCTURAL_JSON) if STRUCTURAL_JSON.exists() else {}
    by_variant: dict[str, list[dict]] = {}

    for row in pairs:
        variant_id = str(row.get("variant_id") or row.get("tag", "unknown").split("__")[0])
        drug_name = str(row.get("drug_name") or row.get("tag", "unknown").split("__")[-1]).replace("_", " ")
        row["variant_id"] = variant_id
        row["drug_name"] = drug_name
        row.setdefault("status", "unknown")
        by_variant.setdefault(variant_id, []).append(row)

    updated_variants = 0
    updated_pairs = 0
    for variant_id, rows in by_variant.items():
        structural.setdefault(variant_id, {})
        structural[variant_id].setdefault("gromacs", {})

        for row in rows:
            structural[variant_id]["gromacs"][row["drug_name"]] = {
                "tag": row.get("tag"),
                "drug_name": row["drug_name"],
                "status": row.get("status"),
                "vina_score_kcalmol": row.get("vina_score_kcalmol"),
                "mmgbsa_delta_g": row.get("mmgbsa_delta_g"),
                "interaction_energy": row.get("interaction_energy"),
                "rmsd_mean_nm": row.get("rmsd_mean_nm"),
                "rmsd_std_nm": row.get("rmsd_std_nm"),
                "rmsf_mean_nm": row.get("rmsf_mean_nm"),
                "hbonds_mean": row.get("hbonds_mean"),
                "hbonds_max": row.get("hbonds_max"),
                "run_dir": row.get("run_dir"),
                "trajectory_path": row.get("trajectory_path"),
                "plot_paths": row.get("plot_paths", {}),
                "notes": row.get("notes"),
                "ingested_at": datetime.now().isoformat()
            }
            updated_pairs += 1

        best = pick_best_pair(rows)
        if best:
            structural[variant_id]["gromacs_best_drug"] = best.get("drug_name")
            structural[variant_id]["gromacs_best_mmgbsa"] = best.get("mmgbsa_delta_g")
            structural[variant_id]["gromacs_best_rmsd_nm"] = best.get("rmsd_mean_nm")
            structural[variant_id]["gromacs_best_hbonds"] = best.get("hbonds_mean")
            updated_variants += 1

    normalized = {
        "generated_at": payload.get("generated_at"),
        "ingested_at": datetime.now().isoformat(),
        "pipeline_version": payload.get("pipeline_version", "GEMAMX-GROMACS-v2"),
        "parameters": payload.get("parameters", {}),
        "summary": payload.get("summary", {}),
        "pairs": pairs
    }

    save_json(structural, STRUCTURAL_JSON)
    save_json(normalized, LOCAL_COPY)
    if input_path.resolve() != DEFAULT_INPUT.resolve():
        shutil.copy2(input_path, DEFAULT_INPUT)

    return {
        "input": str(input_path),
        "updated_variants": updated_variants,
        "updated_pairs": updated_pairs,
        "structural_metrics": str(STRUCTURAL_JSON),
        "local_copy": str(LOCAL_COPY),
        "results_copy": str(DEFAULT_INPUT)
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Ingest GROMACS metrics from Colab or Google Drive")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise SystemExit(f"No se encontro el archivo de entrada: {input_path}")

    summary = ingest_results(input_path)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
