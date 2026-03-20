"""
generate_missing.py
Genera solo las estructuras que faltan y son alcanzables con ESMFold.
Corre rapido — omite todo lo que ya existe.
"""
import requests, os, re

ESMATLAS_URL   = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMATLAS_LIMIT = 390

def mutate(wt_seq, change):
    """Aplica mutacion puntual a la secuencia WT. Retorna secuencia mutada o None."""
    m = re.match(r"([A-Z])(\d+)([A-Z])", change)
    if not m:
        return None
    wt_aa, pos, mut_aa = m.group(1), int(m.group(2)), m.group(3)
    idx = pos - 1
    if idx >= len(wt_seq) or wt_seq[idx] != wt_aa:
        print(f"  ✗ Pos {pos}: esperado '{wt_aa}', encontrado '{wt_seq[idx] if idx < len(wt_seq) else 'OOB'}'")
        return None
    seq_list = list(wt_seq)
    seq_list[idx] = mut_aa
    return "".join(seq_list)

def read_seq(fasta_path):
    lines = open(fasta_path, encoding="utf-8").read().strip().split("\n")
    return "".join(l for l in lines if not l.startswith(">"))

def fragment(seq, mut_pos, window=180):
    start = max(0, mut_pos - window - 1)
    end   = min(len(seq), mut_pos + window)
    return seq[start:end], start + 1

def predict(variant_id, seq, out_pdb):
    print(f"  ESMFold ({len(seq)} aa)...")
    r = requests.post(ESMATLAS_URL,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data=seq, timeout=180)
    if r.ok and "ATOM" in r.text:
        open(out_pdb, "w").write(r.text)
        print(f"  ✓ {out_pdb} ({os.path.getsize(out_pdb)/1024:.1f} KB)")
        return True
    print(f"  ✗ {r.status_code}: {r.text[:150]}")
    return False

os.makedirs("data/structures", exist_ok=True)
os.makedirs("data/fasta", exist_ok=True)

# ── Variantes faltantes alcanzables con ESMFold ───────────────────────────
# Formato: variant_id -> (wt_fasta, change_or_None, mut_pos)
# change=None significa usar la secuencia del wt_fasta directamente (ya es mutante)

PENDIENTES = {
    # TP53 variantes — 393 aa, secuencia completa
    "TP53_C145F": ("data/fasta/TP53_WT.fasta", "C145F", 145),
    "TP53_H20Y":  ("data/fasta/TP53_WT.fasta", "H20Y",   20),
    "TP53_I193S": ("data/fasta/TP53_WT.fasta", "I193S", 193),
    "TP53_P146S": ("data/fasta/TP53_WT.fasta", "P146S", 146),
}

results = {}

for vid, (wt_path, change, mut_pos) in PENDIENTES.items():
    print(f"\n=== {vid} ===")
    out_pdb = f"data/structures/{vid}.pdb"

    if os.path.exists(out_pdb) and os.path.getsize(out_pdb) > 1000:
        print(f"  Ya existe — omitiendo")
        results[vid] = out_pdb
        continue

    if not os.path.exists(wt_path):
        print(f"  ✗ FASTA no encontrado: {wt_path}")
        continue

    wt_seq = read_seq(wt_path)

    # Aplicar mutacion
    seq = mutate(wt_seq, change) if change else wt_seq
    if seq is None:
        continue

    # Guardar FASTA mutante para referencia
    mut_fasta = f"data/fasta/{vid}.fasta"
    if not os.path.exists(mut_fasta):
        open(mut_fasta, "w").write(f">{vid}\n{seq}\n")

    # Secuencia completa o fragmento
    if len(seq) <= ESMATLAS_LIMIT:
        ok = predict(vid, seq, out_pdb)
    else:
        frag, frag_start = fragment(seq, mut_pos)
        print(f"  Fragmento aa {frag_start}-{frag_start+len(frag)-1} ({len(frag)} aa)")
        ok = predict(vid, frag, out_pdb)

    if ok:
        results[vid] = out_pdb

print("\n" + "="*50)
print(f"Generadas: {len(results)}/{len(PENDIENTES)}")
for vid, path in results.items():
    print(f'  "{vid}": "{path.replace(chr(92), "/")}",')
print("\nEjecuta scan_boltz.py para actualizar el manifest.")
