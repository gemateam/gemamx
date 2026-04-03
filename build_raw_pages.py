"""
build_raw_pages.py
==================
Pagina raw_mutations.json en chunks de 500 registros
y los escribe en GEMAMXST/data/raw_pages/ para Netlify.

Uso:
    cd D:\ProyectosIA\GEMAMX
    python build_raw_pages.py
"""

import json, os, shutil

SRC  = "data/raw_mutations.json"
DEST = "GEMAMXST/data/raw_pages"
CHUNK = 500

print(f"Leyendo {SRC}...")
with open(SRC, "r", encoding="utf-8") as f:
    data = json.load(f)
print(f"  {len(data):,} registros cargados")

# Limpiar y recrear destino
if os.path.exists(DEST):
    shutil.rmtree(DEST)
os.makedirs(DEST)

chunks = [data[i:i+CHUNK] for i in range(0, len(data), CHUNK)]

# Escribir index
index = {"total": len(data), "chunk_size": CHUNK, "chunks": len(chunks)}
with open(f"{DEST}/index.json", "w", encoding="utf-8") as f:
    json.dump(index, f)

# Escribir páginas
for i, chunk in enumerate(chunks):
    with open(f"{DEST}/page_{i}.json", "w", encoding="utf-8") as f:
        json.dump(chunk, f, ensure_ascii=False)
    if i % 20 == 0:
        print(f"  Página {i}/{len(chunks)}...")

total_kb = sum(
    os.path.getsize(f"{DEST}/{fn}")
    for fn in os.listdir(DEST)
) / 1024

print(f"\n✓ {len(chunks)} páginas escritas en {DEST}/")
print(f"✓ Tamaño total: {total_kb/1024:.1f} MB")
print(f"✓ vs raw_mutations.json original: {os.path.getsize(SRC)/1024/1024:.1f} MB")
print("\nBuild listo. Sube la carpeta GEMAMXST/ a Netlify.")
