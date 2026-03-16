import requests
import json
import pandas as pd
from pymongo import MongoClient
from tqdm import tqdm
import time

# MongoDB Setup
try:
    client = MongoClient("mongodb://localhost:27017/", serverSelectionTimeoutMS=5000)
    db = client["gen_mama_db"]
    mutations_col = db["mutations"]
    ranking_col = db["top_variants"]
    print("Connected to MongoDB on localhost:27017")
except Exception as e:
    print(f"Error connecting to MongoDB: {e}")
    exit(1)

# Configuration
STUDY_ID = "brca_tcga"
MOLECULAR_PROFILE_ID = "brca_tcga_mutations"
KEY_GENES = ["PIK3CA", "TP53", "ERBB2", "BRCA1", "BRCA2", "CDH1", "GATA3", "MAP3K1"]
BASE_URL = "https://www.cbioportal.org/api"

def get_gene_mapping():
    print("Fetching gene mapping (Symbol -> EntrezID)...")
    genes_map = {}
    for gene in KEY_GENES:
        url = f"{BASE_URL}/genes/{gene}"
        res = requests.get(url)
        if res.status_code == 200:
            data = res.json()
            entrez_id = data.get("entrezGeneId")
            if entrez_id:
                genes_map[entrez_id] = gene
    print(f"Mapped {len(genes_map)} genes.")
    return genes_map

def fetch_mutations():
    print(f"Fetching mutations for study: {STUDY_ID}...")
    url = f"{BASE_URL}/molecular-profiles/{MOLECULAR_PROFILE_ID}/mutations/fetch"
    
    print("Getting sample list...")
    sample_list_url = f"{BASE_URL}/studies/{STUDY_ID}/sample-lists"
    response = requests.get(sample_list_url)
    sample_lists = response.json()
    
    all_samples_id = next((sl["sampleListId"] for sl in sample_lists if "all" in sl["sampleListId"]), None)
    if not all_samples_id:
        all_samples_id = sample_lists[0]["sampleListId"]
    
    print(f"Using sample list: {all_samples_id}")
    
    samples_url = f"{BASE_URL}/sample-lists/{all_samples_id}/sample-ids"
    sample_ids = requests.get(samples_url).json()
    total_samples = len(sample_ids)
    print(f"Found {total_samples} samples.")

    all_mutations = []
    chunk_size = 500
    for i in tqdm(range(0, len(sample_ids), chunk_size), desc="Downloading Mutations"):
        chunk = sample_ids[i:i + chunk_size]
        payload = { "sampleIds": chunk }
        res = requests.post(url, json=payload)
        if res.status_code == 200:
            all_mutations.extend(res.json())
        else:
            print(f"Error fetching chunk: {res.status_code}")
        time.sleep(0.5)

    return all_mutations, total_samples

def process_and_store(data, gene_map, total_samples):
    if not data:
        print("No data received.")
        return

    df = pd.DataFrame(data)
    print(f"Total mutations fetched: {len(df)}")

    # Map entrezGeneId to Symbol
    df['hugo_symbol'] = df['entrezGeneId'].map(gene_map)

    # US-01 Columns mapping
    mapping = {
        'mutationType': 'variant_classification',
        'proteinChange': 'protein_change',
        'sampleId': 'sample_id',
        'chr': 'chromosome',
        'startPosition': 'start_position'
    }
    
    df = df.rename(columns=mapping)

    # US-02: Filter Missense in Key Genes
    print("Filtering key genes and missense variants...")
    filtered_df = df[
        (df['hugo_symbol'].notnull()) & 
        (df['variant_classification'].str.contains('Missense', case=False, na=False))
    ].copy()
    
    print(f"Filtered mutations: {len(filtered_df)}")

    # Store in MongoDB
    print("Storing filtered mutations in MongoDB...")
    mutations_col.delete_many({})
    if not filtered_df.empty:
        mutations_col.insert_many(filtered_df.to_dict('records'))

    # US-03: Ranking Top-10
    print("Calculating Top-10 ranking...")
    if not filtered_df.empty:
        # Group and keep coordinates
        ranking = filtered_df.groupby(['hugo_symbol', 'protein_change']).agg({
            'chromosome': 'first',
            'start_position': 'first'
        }).reset_index()
        
        # Add frequency
        counts = filtered_df.groupby(['hugo_symbol', 'protein_change']).size().reset_index(name='frequency')
        ranking = ranking.merge(counts, on=['hugo_symbol', 'protein_change'])
        
        ranking = ranking.sort_values(by='frequency', ascending=False).head(10)
        
        # Calculate Percentage and Metadata
        ranking['cohort_percentage'] = (ranking['frequency'] / total_samples) * 100
        
        # Heuristic for classification (simple for now: Frequent Hotspots = Pathogenic)
        def classify(freq):
            if freq > 50: return "Patogénica"
            if freq > 10: return "Probable Patogénica"
            return "VUS"
            
        ranking['classification'] = ranking['frequency'].apply(classify)
        ranking['rank'] = range(1, len(ranking) + 1)
        
        print("Top 10 Variants with Metrics:")
        print(ranking)

        ranking_col.delete_many({})
        ranking_col.insert_many(ranking.to_dict('records'))
        print("Ranking stored in MongoDB.")
    else:
        print("No variants to rank.")

if __name__ == "__main__":
    gene_map = get_gene_mapping()
    raw_data, total_samples = fetch_mutations()
    process_and_store(raw_data, gene_map, total_samples)
