"""
rebuild_stats.py
Regenera mutation_stats.json desde raw_mutations.json (el dataset completo).
Reemplaza el archivo incorrecto generado por Colab.
"""
import json, collections

print("Leyendo raw_mutations.json...")
with open("data/raw_mutations.json", encoding="utf-8") as f:
    mutations = json.load(f)

print(f"Total registros: {len(mutations):,}")

# Contar por gen
gene_counts = collections.Counter(m["hugo_symbol"] for m in mutations if m.get("hugo_symbol"))
top_genes   = dict(gene_counts.most_common(20))

# Contar por tipo de variante
type_counts = collections.Counter(m["variant_classification"] for m in mutations if m.get("variant_classification"))
type_dist   = dict(sorted(type_counts.items(), key=lambda x: x[1], reverse=True))

# Contar muestras unicas
samples = {m["sample_id"] for m in mutations if m.get("sample_id")}

stats = {
    "gene_distribution":  top_genes,
    "type_distribution":  type_dist,
    "total_samples":      len(samples),
    "total_mutations":    len(mutations)
}

with open("data/mutation_stats.json", "w", encoding="utf-8") as f:
    json.dump(stats, f, indent=4, ensure_ascii=False)

print(f"\n✓ mutation_stats.json reconstruido:")
print(f"  Muestras unicas:  {len(samples):,}")
print(f"  Total mutaciones: {len(mutations):,}")
print(f"\nTop 5 genes:")
for g, c in list(top_genes.items())[:5]:
    print(f"  {g:15s}: {c:,}")
print(f"\nTipos de variante ({len(type_dist)} tipos):")
for t, c in list(type_dist.items())[:5]:
    print(f"  {t:30s}: {c:,}")
