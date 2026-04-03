#!/usr/bin/env python3
"""
GEMAMX v2 - Scan Boltz/ESMFold structures and refresh manifest.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STRUCTURES_DIR = ROOT / "v2" / "data" / "structures"
MANIFEST_PATH = ROOT / "v2" / "data" / "manifest_structures.json"


def main():
    manifest = {}
    if MANIFEST_PATH.exists():
        try:
            manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        except Exception:
            manifest = {}

    found = {}
    for ext in ("*.cif", "*.pdb"):
        for path in sorted(STRUCTURES_DIR.rglob(ext)):
            if path.parent.name == "meta":
                continue
            if path.stat().st_size < 1000:
                continue
            key = path.stem
            found[key] = str(path.relative_to(ROOT)).replace('\\', '/')

    manifest.update(found)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    print("=" * 55)
    print("SCAN DE ESTRUCTURAS v2")
    print("=" * 55)
    for key, rel in sorted(found.items()):
        size_kb = round((ROOT / rel).stat().st_size / 1024, 1)
        print(f"  OK {key:35s} {size_kb:7.1f} KB  {rel}")
    print(f"\nTotal estructuras encontradas: {len(found)}")
    print(f"Manifest actualizado: {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
