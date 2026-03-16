import requests
import os
from pymongo import MongoClient
from Bio import SeqIO
from Bio.Seq import MutableSeq
import io

# MongoDB Setup
try:
    client = MongoClient("mongodb://localhost:27017/", serverSelectionTimeoutMS=5000)
    db = client["gen_mama_db"]
    ranking_col = db["top_variants"]
    fasta_col = db["fasta_sequences"]
    print("Connected to MongoDB on localhost:27017")
except Exception as e:
    print(f"Error connecting to MongoDB: {e}")
    exit(1)

# Ensure output directory exists
os.makedirs("data/fasta", exist_ok=True)

def get_uniprot_id(gene_symbol):
    """Search UniProt for the best accession ID for a human gene symbol."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f'gene:{gene_symbol} AND organism_id:9606 AND reviewed:true',
        "fields": "accession",
        "format": "json",
        "size": 1
    }
    res = requests.get(url, params=params)
    if res.status_code == 200:
        data = res.json()
        if data.get("results"):
            return data["results"][0]["primaryAccession"]
    return None

def fetch_fasta(accession):
    """Fetch the FASTA sequence from UniProt."""
    url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
    res = requests.get(url)
    if res.status_code == 200:
        return res.text
    return None

def apply_mutation(wt_fasta, protein_change):
    """
    Apply a missense mutation (e.g., H1047R) to a FASTA sequence.
    H1047R means: Pos 1047, H (WT) -> R (Mutant)
    """
    try:
        # Parse FASTA
        record = SeqIO.read(io.StringIO(wt_fasta), "fasta")
        wt_seq = record.seq
        
        # Parse mutation string (e.g., H1047R)
        # Use regex or simple split
        import re
        match = re.match(r"([A-Z])(\d+)([A-Z])", protein_change)
        if not match:
            # Handle cases like p.H1047R
            match = re.search(r"([A-Z])(\d+)([A-Z])", protein_change)
            
        if match:
            wt_aa = match.group(1)
            pos = int(match.group(2))
            mut_aa = match.group(3)
            
            # Position is 1-indexed in biological notation
            idx = pos - 1
            
            if idx < len(wt_seq) and wt_seq[idx] == wt_aa:
                mut_seq = MutableSeq(wt_seq)
                mut_seq[idx] = mut_aa
                return str(wt_seq), str(mut_seq), pos
            else:
                print(f"Warning: Position {pos} mismatch or out of bounds for {protein_change}. WT char at index: {wt_seq[idx] if idx < len(wt_seq) else 'OOB'}")
        
    except Exception as e:
        print(f"Error applying mutation {protein_change}: {e}")
    
    return None, None, None

def run():
    variants = list(ranking_col.find().sort("rank", 1))
    print(f"Processing {len(variants)} variants from ranking...")
    
    results = []
    for var in variants:
        gene = var["hugo_symbol"]
        change = var["protein_change"]
        print(f"Processing {gene} {change}...")
        
        acc = get_uniprot_id(gene)
        if not acc:
            print(f"Could not find UniProt ID for {gene}")
            continue
            
        wt_fasta = fetch_fasta(acc)
        if not wt_fasta:
            print(f"Could not fetch FASTA for {acc}")
            continue
            
        wt_str, mut_str, pos = apply_mutation(wt_fasta, change)
        
        if wt_str and mut_str:
            # Save to MongoDB
            doc = {
                "gene": gene,
                "protein_change": change,
                "uniprot_id": acc,
                "wt_sequence": wt_str,
                "mut_sequence": mut_str,
                "mutation_pos": pos
            }
            results.append(doc)
            
            # Save to disk for AlphaFold/Boltz
            filename_wt = f"data/fasta/{gene}_WT.fasta"
            filename_mut = f"data/fasta/{gene}_{change}.fasta"
            
            with open(filename_wt, "w") as f:
                f.write(f">{gene}_WT | {acc}\n{wt_str}\n")
            with open(filename_mut, "w") as f:
                f.write(f">{gene}_{change} | {acc}\n{mut_str}\n")
                
            print(f"Saved FASTA for {gene} {change}")

    if results:
        fasta_col.delete_many({})
        fasta_col.insert_many(results)
        print(f"Total {len(results)} sequences stored in MongoDB.")

if __name__ == "__main__":
    run()
