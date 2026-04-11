from __future__ import annotations

import json
from pathlib import Path

from build_gromacs_batch import MANIFEST_JSON, RESULTS_DIR, build_manifest


PILOT_VARIANTS = [
    "PIK3CA_H1047R",
    "ERBB2_V777L",
    "TP53_R248Q",
]

PILOT_PRIORITY = {
    "PIK3CA_H1047R": ["Alpelisib", "Pictilisib"],
    "ERBB2_V777L": ["Pictilisib", "Alpelisib"],
    "TP53_R248Q": ["Abemaciclib", "Pictilisib"],
}

PILOT_OUTPUT = RESULTS_DIR / "gromacs_pilot_manifest.json"


def main() -> None:
    manifest = build_manifest(top_n_per_variant=5)
    pilot_pairs = []

    for row in manifest["pairs"]:
        variant_id = row["variant_id"]
        drug_name = row["drug_name"]
        if variant_id not in PILOT_VARIANTS:
            continue
        if drug_name not in PILOT_PRIORITY.get(variant_id, []):
            continue
        pilot_pairs.append(row)

    pilot_pairs.sort(
        key=lambda row: (
            PILOT_VARIANTS.index(row["variant_id"]),
            PILOT_PRIORITY[row["variant_id"]].index(row["drug_name"]),
        )
    )

    payload = {
        "generated_at": manifest["generated_at"],
        "pipeline_version": "GEMAMX-GROMACS-v2-pilot",
        "selection_policy": {
            "mode": "pilot curated subset",
            "source_manifest": str(MANIFEST_JSON),
            "variants": PILOT_VARIANTS,
            "priority_by_variant": PILOT_PRIORITY,
        },
        "summary": {
            "variants_considered": len(PILOT_VARIANTS),
            "pairs_selected": len(pilot_pairs),
            "per_variant_counts": {
                variant: sum(1 for row in pilot_pairs if row["variant_id"] == variant)
                for variant in PILOT_VARIANTS
            },
        },
        "pairs": pilot_pairs,
    }

    PILOT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    PILOT_OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
