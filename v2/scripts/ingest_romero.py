#!/usr/bin/env python3
"""
GEMAMX v2 - Ingestion pipeline for Romero cohort.

Source files:
- data/romero/41467_2021_22478_MOESM8_ESM.xlsx  -> sheet 'HM maf'
- data/romero/41467_2021_22478_MOESM9_ESM.xlsx  -> sheet 'Driver mutation in HM'

Outputs in v2/data:
- raw_mutations.json
- top_variants.json
- mutation_stats.json
"""

from __future__ import annotations

import json
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
ROMERO_DIR = ROOT / "data" / "romero"
OUT_DIR = ROOT / "v2" / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MOESM8 = ROMERO_DIR / "41467_2021_22478_MOESM8_ESM.xlsx"
MOESM9 = ROMERO_DIR / "41467_2021_22478_MOESM9_ESM.xlsx"

KEY_GENES = [
    "PIK3CA", "TP53", "ERBB2", "BRCA1", "BRCA2",
    "CDH1", "GATA3", "MAP3K1", "KMT2C", "PTEN",
    "RB1", "ARID1A", "SF3B1", "AKT1", "CBFB"
]

NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkg": "http://schemas.openxmlformats.org/package/2006/relationships",
}


def load_xlsx_sheet_rows(path: Path, sheet_name: str):
    with zipfile.ZipFile(path) as zf:
        shared_strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("main:si", NS):
                text = "".join(t.text or "" for t in si.findall(".//main:t", NS))
                shared_strings.append(text)

        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels.findall("pkg:Relationship", NS)}

        target = None
        for sheet in workbook.find("main:sheets", NS):
            if sheet.attrib.get("name") == sheet_name:
                rid = sheet.attrib.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
                target = "xl/" + rel_map[rid]
                break
        if not target:
            raise ValueError(f"Sheet not found: {sheet_name}")

        sheet_xml = ET.fromstring(zf.read(target))
        rows = []
        for row in sheet_xml.findall(".//main:sheetData/main:row", NS):
            values = []
            for cell in row.findall("main:c", NS):
                cell_type = cell.attrib.get("t")
                v = cell.find("main:v", NS)
                value = "" if v is None or v.text is None else v.text
                if cell_type == "s" and value.isdigit():
                    idx = int(value)
                    value = shared_strings[idx] if idx < len(shared_strings) else value
                values.append(value)
            rows.append(values)
        return rows


def sheet_to_dicts(path: Path, sheet_name: str, header_row_idx: int = 1):
    rows = load_xlsx_sheet_rows(path, sheet_name)
    header = rows[header_row_idx]
    data_rows = rows[header_row_idx + 1 :]
    width = len(header)
    cleaned = []
    for row in data_rows:
        if not any(str(v).strip() for v in row):
            continue
        padded = row + [""] * (width - len(row))
        cleaned.append({str(header[i]).strip(): padded[i] for i in range(width)})
    return cleaned


def parse_int(value, default=0):
    try:
        return int(float(str(value).strip()))
    except Exception:
        return default


def parse_gdna(gdna: str):
    # Example: chr3:g.178952085A>G
    m = re.match(r"chr(?P<chr>[^:]+):g\.(?P<pos>\d+)(?P<ref>[A-Z\-]+)>(?P<alt>[A-Z\-]+)", gdna or "")
    if not m:
        return None
    return {
        "chromosome": m.group("chr"),
        "start_position": int(m.group("pos")),
        "reference_allele": m.group("ref"),
        "tumor_seq_allele2": m.group("alt"),
    }


def classify_variant(cohort_pct: float) -> str:
    if cohort_pct >= 5.0:
        return "Patogénica"
    if cohort_pct >= 1.0:
        return "Probable Patogénica"
    if cohort_pct >= 0.5:
        return "VUS"
    return "Baja Frecuencia"


def build_driver_lookup(driver_rows):
    lookup = {}
    for row in driver_rows:
        gdna = row.get("gdna", "")
        parsed = parse_gdna(gdna)
        if not parsed:
            continue
        gene = str(row.get("Gene", "")).strip()
        key = (
            gene,
            parsed["chromosome"],
            parsed["start_position"],
            parsed["reference_allele"],
            parsed["tumor_seq_allele2"],
        )
        lookup[key] = {
            "protein": str(row.get("Protein", "")).strip(),
            "protein_change": normalize_protein_change(row.get("Protein", ""), row.get("Protein change", "")),
            "driver": str(row.get("Driver", "")).strip(),
            "consequence": str(row.get("Consequence", "")).strip(),
            "gene_role_mutation": str(row.get("Gene role mutation", "")).strip(),
            "known_variant_reference": str(row.get("Known variant reference", "")).strip(),
        }
    return lookup


def normalize_protein_change(protein_value: str, short_change: str) -> str | None:
    protein_value = str(protein_value or "").strip()
    short_change = str(short_change or "").strip()
    m = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", protein_value)
    aa3_to_1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
        "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
        "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    }
    if m:
        wt = aa3_to_1.get(m.group(1), "")
        pos = m.group(2)
        mt = aa3_to_1.get(m.group(3), "")
        if wt and mt:
            return f"{wt}{pos}{mt}"
    if short_change and re.match(r"^[A-Z]/[A-Z]$", short_change):
        pos_match = re.search(r"(\d+)", protein_value)
        if pos_match:
            wt, mt = short_change.split("/")
            return f"{wt}{pos_match.group(1)}{mt}"
    return None


def main():
    hm_maf_rows = sheet_to_dicts(MOESM8, "HM maf", header_row_idx=1)
    driver_rows = sheet_to_dicts(MOESM9, "Driver mutation in HM", header_row_idx=1)

    driver_lookup = build_driver_lookup(driver_rows)

    raw_mutations = []
    for row in hm_maf_rows:
        gene = str(row.get("Hugo Symbol", "")).strip()
        chrom = str(row.get("Chromosome", "")).strip()
        start = parse_int(row.get("Start Position", 0))
        ref = str(row.get("Reference Allele", "")).strip()
        alt = str(row.get("Tumor Seq Allele2", "")).strip() or str(row.get("Tumor Seq Allele1", "")).strip()
        lookup_key = (gene, chrom, start, ref, alt)
        driver = driver_lookup.get(lookup_key, {})

        raw_mutations.append({
            "hugo_symbol": gene,
            "variant_classification": str(row.get("Variant Classification", "")).strip(),
            "variant_type": str(row.get("Variant Type", "")).strip(),
            "protein_change": driver.get("protein_change"),
            "protein_hgvsp": driver.get("protein"),
            "sample_id": str(row.get("Tumor Sample Barcode", "")).strip(),
            "chromosome": chrom,
            "start_position": start,
            "end_position": parse_int(row.get("End Position", 0)),
            "strand": str(row.get("Strand", "")).strip(),
            "reference_allele": ref,
            "tumor_seq_allele1": str(row.get("Tumor Seq Allele1", "")).strip(),
            "tumor_seq_allele2": alt,
            "source_center": str(row.get("Center", "")).strip(),
            "driver_status": driver.get("driver"),
            "consequence": driver.get("consequence"),
            "gene_role_mutation": driver.get("gene_role_mutation"),
            "known_variant_reference": driver.get("known_variant_reference"),
        })

    total_samples = len({m["sample_id"] for m in raw_mutations if m.get("sample_id")})

    filtered = [
        m for m in raw_mutations
        if m.get("hugo_symbol") in KEY_GENES
        and (m.get("variant_classification") or "").lower().find("missense") >= 0
        and m.get("protein_change")
    ]

    counts = defaultdict(lambda: {"frequency": 0, "chromosome": "?", "start_position": 0})
    for m in filtered:
        key = (m["hugo_symbol"], m["protein_change"])
        counts[key]["frequency"] += 1
        counts[key]["chromosome"] = m.get("chromosome") or counts[key]["chromosome"]
        counts[key]["start_position"] = m.get("start_position") or counts[key]["start_position"]

    top_variants = []
    for rank, ((gene, protein_change), meta) in enumerate(sorted(counts.items(), key=lambda kv: kv[1]["frequency"], reverse=True), start=1):
        cohort_pct = (meta["frequency"] / total_samples * 100.0) if total_samples else 0.0
        top_variants.append({
            "hugo_symbol": gene,
            "protein_change": protein_change,
            "chromosome": meta["chromosome"],
            "start_position": meta["start_position"],
            "frequency": meta["frequency"],
            "cohort_percentage": cohort_pct,
            "classification": classify_variant(cohort_pct),
            "rank": rank,
            "dataset": "Romero-Cordoba 2021 HM",
        })

    gene_dist = dict(Counter(m["hugo_symbol"] for m in raw_mutations if m.get("hugo_symbol")).most_common(20))
    type_dist = dict(Counter(m["variant_classification"] for m in raw_mutations if m.get("variant_classification")).most_common())
    stats = {
        "gene_distribution": gene_dist,
        "type_distribution": type_dist,
        "total_samples": total_samples,
        "total_mutations": len(raw_mutations),
        "dataset": "Romero-Cordoba 2021 HM",
    }

    (OUT_DIR / "raw_mutations.json").write_text(json.dumps(raw_mutations, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUT_DIR / "top_variants.json").write_text(json.dumps(top_variants, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUT_DIR / "mutation_stats.json").write_text(json.dumps(stats, indent=2, ensure_ascii=False), encoding="utf-8")

    print("GEMAMX v2 Romero ingestion complete")
    print(f"  raw_mutations: {len(raw_mutations)}")
    print(f"  filtered top variants: {len(top_variants)}")
    print(f"  total samples: {total_samples}")
    if top_variants:
        print("  top 10 variants:")
        for item in top_variants[:10]:
            print(f"    #{item['rank']} {item['hugo_symbol']} {item['protein_change']} -> {item['frequency']} cases ({item['cohort_percentage']:.2f}%)")


if __name__ == "__main__":
    main()
