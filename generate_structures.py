"""
generate_structures.py
Genera estructuras PDB para variantes mutantes usando ESMFold API.
- Proteinas cortas (<400aa): secuencia completa
- Proteinas largas (>400aa): fragmento de dominio centrado en la mutacion
- Genera FASTA faltantes desde UniProt si no existen en disco

Variantes configuradas:
  Batch 1 (originales): AARS1_S645C, PIK3CA_H1047R, ERBB2_V777L, BRCA1_C61G
  Batch 2 (nuevas de Colab): TP53 variantes, PTEN_N48I, PIK3CA variantes adicionales,
                              SF3B1_E802Q, SF3B1_WT
"""

import requests, os, re, io

try:
    from Bio import SeqIO
    from Bio.Seq import MutableSeq
    HAS_BIO = True
except ImportError:
    HAS_BIO = False
    print("AVISO: biopython no disponible")

ESMATLAS_URL   = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMATLAS_LIMIT = 390

def get_domain_fragment(seq, mut_pos, window=180):
    start = max(0, mut_pos - window - 1)
    end   = min(len(seq), mut_pos + window)
    return seq[start:end], start + 1

def get_uniprot_fasta(uniprot_id):
    r = requests.get(f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta", timeout=30)
    return r.text if r.ok else None

def apply_mutation_and_save(gene, uniprot_id, change, out_fasta):
    if not HAS_BIO:
        print(f"  ✗ biopython requerido"); return False
    wt_fasta = get_uniprot_fasta(uniprot_id)
    if not wt_fasta:
        print(f"  ✗ No se pudo obtener FASTA de {uniprot_id}"); return False
    record = SeqIO.read(io.StringIO(wt_fasta), "fasta")
    seq = MutableSeq(str(record.seq))
    m = re.match(r"([A-Z])(\d+)([A-Z])", change)
    if m:
        pos = int(m.group(2))
        if pos <= len(seq) and seq[pos-1] == m.group(1):
            seq[pos-1] = m.group(3)
        else:
            print(f"  ✗ Posicion {pos} no coincide"); return False
    with open(out_fasta, "w") as f:
        f.write(f">{gene}_{change} | {uniprot_id}\n{str(seq)}\n")
    print(f"  ✓ FASTA generado: {out_fasta} ({len(seq)} aa)")
    return True

def predict_structure(variant_id, seq, out_pdb):
    print(f"  Prediciendo con ESMFold ({len(seq)} aa)...")
    r = requests.post(
        ESMATLAS_URL,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data=seq, timeout=180
    )
    if r.ok and "ATOM" in r.text:
        open(out_pdb, "w").write(r.text)
        size_kb = os.path.getsize(out_pdb) / 1024
        print(f"  ✓ Guardado: {out_pdb} ({size_kb:.1f} KB)")
        return True
    else:
        print(f"  ✗ ESMFold error {r.status_code}: {r.text[:200]}")
        return False

# ─────────────────────────────────────────────────────────────────────────
# Configuracion de variantes
#
# Para variantes con fasta_path existente: lee la secuencia del disco
# Para variantes sin fasta: descarga de UniProt y aplica mutacion
#
# Proteinas grandes (>390 aa): usa fragmento de dominio centrado en mut_pos
# Proteinas pequeñas (<=390 aa): secuencia completa
# ─────────────────────────────────────────────────────────────────────────

VARIANTES = {
    # ── Batch 1: originales (ya generadas, se omitiran) ───────────────────
    "AARS1_S645C": {
        "fasta_path": "data/fasta/AARS1_S645C.fasta",
        "mut_pos": 645,
    },
    "PIK3CA_H1047R": {
        "fasta_path": "data/fasta/PIK3CA_H1047R.fasta",
        "mut_pos": 1047,
    },
    "ERBB2_V777L": {
        "fasta_path": "data/fasta/ERBB2_V777L.fasta",
        "uniprot_id": "P04626",
        "mut_change": "V777L",
        "mut_pos": 777,
    },
    "BRCA1_C61G": {
        "fasta_path": "data/fasta/BRCA1_C61G.fasta",
        "uniprot_id": "P38398",
        "mut_change": "C61G",
        "mut_pos": 61,
    },

    # ── Batch 2: variantes TP53 (393 aa — secuencia completa) ────────────
    # TP53_R248Q ya existe, las siguientes son nuevas
    "TP53_C145F": {
        "fasta_path": "data/fasta/TP53_WT.fasta",   # usar WT como base y aplicar mutacion
        "mut_change": "C145F",
        "mut_pos": 145,
        "use_wt_fasta": True,   # flag: aplicar mutacion sobre el WT existente
    },
    "TP53_H20Y": {
        "fasta_path": "data/fasta/TP53_WT.fasta",
        "mut_change": "H20Y",
        "mut_pos": 20,
        "use_wt_fasta": True,
    },
    "TP53_I193S": {
        "fasta_path": "data/fasta/TP53_WT.fasta",
        "mut_change": "I193S",
        "mut_pos": 193,
        "use_wt_fasta": True,
    },
    "TP53_P146S": {
        "fasta_path": "data/fasta/TP53_WT.fasta",
        "mut_change": "P146S",
        "mut_pos": 146,
        "use_wt_fasta": True,
    },

    # ── PTEN_N48I (403 aa — fragmento de dominio) ─────────────────────────
    "PTEN_N48I": {
        "fasta_path": "data/fasta/PTEN_N48I.fasta",
        "uniprot_id": "P60484",
        "mut_change": "N48I",
        "mut_pos": 48,
    },

    # ── PIK3CA variantes adicionales (1068 aa — fragmento de dominio) ─────
    "PIK3CA_D350N": {
        "fasta_path": "data/fasta/PIK3CA_WT.fasta",
        "mut_change": "D350N",
        "mut_pos": 350,
        "use_wt_fasta": True,
    },
    "PIK3CA_R108H": {
        "fasta_path": "data/fasta/PIK3CA_WT.fasta",
        "mut_change": "R108H",
        "mut_pos": 108,
        "use_wt_fasta": True,
    },

    # ── SF3B1 (1304 aa — fragmento de dominio) ────────────────────────────
    "SF3B1_WT": {
        "fasta_path": "data/fasta/SF3B1_WT.fasta",
        "uniprot_id": "O75533",
        "mut_change": None,   # None = solo WT, no mutar
        "mut_pos": 802,
    },
    "SF3B1_E802Q": {
        "fasta_path": "data/fasta/SF3B1_E802Q.fasta",
        "uniprot_id": "O75533",
        "mut_change": "E802Q",
        "mut_pos": 802,
    },
}

os.makedirs("data/structures", exist_ok=True)
os.makedirs("data/fasta", exist_ok=True)

results = {}

for variant_id, cfg in VARIANTES.items():
    print(f"\n=== {variant_id} ===")
    out_pdb = f"data/structures/{variant_id}.pdb"

    # Omitir si ya existe y es valido
    if os.path.exists(out_pdb) and os.path.getsize(out_pdb) > 1000:
        print(f"  Ya existe y parece valido — omitiendo")
        results[variant_id] = out_pdb
        continue

    fasta_path = cfg["fasta_path"]
    mut_change = cfg.get("mut_change")
    mut_pos    = cfg.get("mut_pos", 1)
    use_wt     = cfg.get("use_wt_fasta", False)

    # ── Caso 1: usar FASTA WT existente + aplicar mutacion en memoria ─────
    if use_wt and os.path.exists(fasta_path) and mut_change:
        lines = open(fasta_path, encoding="utf-8").read().strip().split("\n")
        wt_seq = "".join(l for l in lines if not l.startswith(">"))
        m = re.match(r"([A-Z])(\d+)([A-Z])", mut_change)
        if m:
            pos = int(m.group(2))
            if pos <= len(wt_seq) and wt_seq[pos-1] == m.group(1):
                seq_list = list(wt_seq)
                seq_list[pos-1] = m.group(3)
                seq = "".join(seq_list)
                # Guardar el FASTA mutante
                mut_fasta = f"data/fasta/{variant_id}.fasta"
                if not os.path.exists(mut_fasta):
                    open(mut_fasta, "w").write(f">{variant_id}\n{seq}\n")
                    print(f"  ✓ FASTA mutante generado: {mut_fasta}")
            else:
                print(f"  ✗ Posicion {pos}: esperado {m.group(1)}, encontrado {wt_seq[pos-1] if pos <= len(wt_seq) else 'OOB'}")
                continue
        else:
            seq = wt_seq  # No mutar (para WT de proteina grande)

    # ── Caso 2: FASTA no existe — generarlo desde UniProt ─────────────────
    elif not os.path.exists(fasta_path):
        print(f"  FASTA no encontrado: {fasta_path}")
        if "uniprot_id" in cfg and mut_change:
            print(f"  Generando desde UniProt ({cfg['uniprot_id']})...")
            gene = variant_id.split("_")[0]
            ok = apply_mutation_and_save(gene, cfg["uniprot_id"], mut_change, fasta_path)
            if not ok: continue
        elif "uniprot_id" in cfg and not mut_change:
            # Solo WT — descargar directamente
            wt_text = get_uniprot_fasta(cfg["uniprot_id"])
            if not wt_text: print("  ✗ No se pudo descargar"); continue
            if HAS_BIO:
                record = SeqIO.read(io.StringIO(wt_text), "fasta")
                seq_str = str(record.seq)
            else:
                lines = wt_text.strip().split("\n")
                seq_str = "".join(l for l in lines if not l.startswith(">"))
            open(fasta_path, "w").write(f">{variant_id}\n{seq_str}\n")
            print(f"  ✓ FASTA WT guardado: {fasta_path} ({len(seq_str)} aa)")
        else:
            print("  ✗ Sin uniprot_id ni fasta — saltando")
            continue

    # ── Leer secuencia del FASTA ──────────────────────────────────────────
    if not use_wt or not mut_change:
        lines = open(fasta_path, encoding="utf-8").read().strip().split("\n")
        seq = "".join(l for l in lines if not l.startswith(">"))
    # (si use_wt con mutacion, seq ya fue generada arriba)

    # ── Decidir secuencia completa o fragmento de dominio ─────────────────
    if len(seq) <= ESMATLAS_LIMIT:
        print(f"  Secuencia completa ({len(seq)} aa)")
        ok = predict_structure(variant_id, seq, out_pdb)
    else:
        fragment, frag_start = get_domain_fragment(seq, mut_pos, window=180)
        print(f"  Proteina larga ({len(seq)} aa) — fragmento: aa {frag_start}–{frag_start+len(fragment)-1} ({len(fragment)} aa)")
        ok = predict_structure(variant_id, fragment, out_pdb)

    if ok:
        results[variant_id] = out_pdb

# ── Resumen ───────────────────────────────────────────────────────────────
print("\n" + "="*50)
print("RESUMEN — pegar en manifest_structures.json:")
print("="*50)
for variant_id, pdb_path in results.items():
    print(f'  "{variant_id}": "{pdb_path.replace(chr(92), "/")}",')
if not results:
    print("  Ninguna estructura generada.")

print("\nEjecuta scan_boltz.py para actualizar el manifest completo.")
