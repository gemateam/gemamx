#!/usr/bin/env python3
"""
GEMAMX v2 - FASTA generation for Romero cohort.

Input:
- v2/data/top_variants.json

Outputs:
- v2/data/fasta/*.fasta
- v2/data/fasta_sequences.json
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
TOP_VARIANTS_PATH = ROOT / "v2" / "data" / "top_variants.json"
OUT_DIR = ROOT / "v2" / "data" / "fasta"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_INDEX = ROOT / "v2" / "data" / "fasta_sequences.json"

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def get_uniprot_id(gene_symbol: str) -> str | None:
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,gene_primary,protein_name",
        "format": "json",
        "size": 3,
    }
    resp = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    results = data.get("results", [])
    if not results:
        return None
    return results[0].get("primaryAccession")


def fetch_fasta(accession: str) -> str | None:
    resp = requests.get(UNIPROT_FASTA_URL.format(accession=accession), timeout=30)
    if not resp.ok:
        return None
    return resp.text


def parse_fasta_sequence(fasta_text: str) -> str:
    lines = [line.strip() for line in fasta_text.splitlines() if line.strip()]
    return "".join(line for line in lines if not line.startswith(">"))


def apply_mutation(wt_fasta_text: str, protein_change: str):
    wt_seq = parse_fasta_sequence(wt_fasta_text)
    match = re.search(r"([A-Z])(\d+)([A-Z])", protein_change or "")
    if not match:
        return None, None, None

    wt_aa = match.group(1)
    pos = int(match.group(2))
    mut_aa = match.group(3)
    idx = pos - 1

    if idx < 0 or idx >= len(wt_seq):
        return None, None, None
    if wt_seq[idx] != wt_aa:
        return None, None, None

    seq_list = list(wt_seq)
    seq_list[idx] = mut_aa
    return wt_seq, "".join(seq_list), pos


def main():
    variants = json.loads(TOP_VARIANTS_PATH.read_text(encoding="utf-8"))
    results = []

    for var in variants:
        gene = var["hugo_symbol"]
        change = var["protein_change"]
        print(f"Processing {gene} {change}...")

        accession = get_uniprot_id(gene)
        if not accession:
            print(f"  Could not find UniProt accession for {gene}")
            continue

        wt_fasta = fetch_fasta(accession)
        if not wt_fasta:
            print(f"  Could not fetch FASTA for {accession}")
            continue

        wt_seq, mut_seq, pos = apply_mutation(wt_fasta, change)
        if not wt_seq or not mut_seq:
            print(f"  Could not apply mutation {change} to {gene}")
            continue

        wt_path = OUT_DIR / f"{gene}_WT.fasta"
        mut_path = OUT_DIR / f"{gene}_{change}.fasta"

        if not wt_path.exists():
            wt_path.write_text(f">{gene}_WT | {accession}\n{wt_seq}\n", encoding="utf-8")
        mut_path.write_text(f">{gene}_{change} | {accession}\n{mut_seq}\n", encoding="utf-8")

        results.append({
            "gene": gene,
            "protein_change": change,
            "uniprot_id": accession,
            "mutation_pos": pos,
            "wt_fasta": str(wt_path.relative_to(ROOT)).replace('\\', '/'),
            "mut_fasta": str(mut_path.relative_to(ROOT)).replace('\\', '/'),
        })
        print(f"  Saved WT and MUT FASTA for {gene} {change}")

    OUT_INDEX.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nGenerated {len(results)} FASTA pairs")
    print(f"Index written to {OUT_INDEX}")


if __name__ == "__main__":
    main()
