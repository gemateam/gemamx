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
KEY_GENES = [
    "PIK3CA", "TP53", "ERBB2", "BRCA1", "BRCA2",   # Tier 1: accionables
    "CDH1", "GATA3", "MAP3K1", "KMT2C", "PTEN",     # Tier 2: drivers frecuentes
    "RB1", "ARID1A", "SF3B1", "AKT1", "CBFB"        # Tier 3: relevantes TCGA-BRCA
]
BASE_URL = "https://www.cbioportal.org/api"

def get_gene_mapping(entrez_ids):
    print(f"Fetching mapping for {len(entrez_ids)} unique Gene IDs...")
    genes_map = {}
    url = f"{BASE_URL}/genes/fetch"
    # cBioPortal supports fetching multiple genes via POST
    chunk_size = 500
    for i in range(0, len(entrez_ids), chunk_size):
        chunk = [int(x) for x in entrez_ids[i:i + chunk_size] if pd.notna(x)]
        if not chunk: continue
        res = requests.post(url, json=chunk)
        if res.status_code == 200:
            genes_data = res.json()
            for g in genes_data:
                genes_map[g['entrezGeneId']] = g['hugoGeneSymbol']
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

    # Map entrezGeneId to Symbol for ALL genes
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

    # US-02: Filter Missense in Key Genes for the RANKING and MongoDB
    print("Filtering key genes and missense variants for Ranking...")
    filtered_df = df[
        (df['hugo_symbol'].isin(KEY_GENES)) & 
        (df['variant_classification'].str.contains('Missense', case=False, na=False))
    ].copy()
    
    print(f"Filtered mutations (Key Genes + Missense): {len(filtered_df)}")

    # Store in MongoDB (Only filtered ones for Epic 1 focus)
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
        
        ranking = ranking.sort_values(by='frequency', ascending=False).head(20)
        
        # Calculate Percentage and Metadata
        ranking['cohort_percentage'] = (ranking['frequency'] / total_samples) * 100
        
        # Improved classification based on cohort percentage
        def classify(freq, cohort_pct):
            if cohort_pct >= 5.0:  return "Patogénica"           # ≥5% cohorte
            if cohort_pct >= 1.0:  return "Probable Patogénica"  # 1-5%
            if cohort_pct >= 0.5:  return "VUS"                  # 0.5-1%
            return "Baja Frecuencia"                              # <0.5%
            
        ranking['classification'] = ranking.apply(
            lambda r: classify(r['frequency'], r['cohort_percentage']), axis=1
        )
        ranking['rank'] = range(1, len(ranking) + 1)
        
        print(f"Top {len(ranking)} Variants with Metrics:")
        print(ranking)

        # Export for Web UI Ranking (TCGA Exploration page)
        ranking_export = ranking[['hugo_symbol', 'protein_change', 'chromosome', 
                                   'start_position', 'frequency', 'cohort_percentage', 
                                   'classification', 'rank']].copy()
        ranking_export.to_json('data/top_variants.json', orient='records', indent=4)
        print(f"Top variants exported to data/top_variants.json ({len(ranking_export)} records)")

        ranking_col.delete_many({})
        ranking_col.insert_many(ranking.to_dict('records'))
        print("Ranking stored in MongoDB.")

    # US-06 & US-01: Exporting stats and raw data for Web UI
    print("Exporting data for Web UI...")
    
    # Export for US-01: Full Consultation page (Showing all mapped mutations in the cohort)
    raw_cols = ['hugo_symbol', 'variant_classification', 'protein_change', 'sample_id', 'chromosome', 'start_position']
    available_cols = [c for c in raw_cols if c in df.columns]
    raw_export = df[available_cols]
    raw_export.to_json('data/raw_mutations.json', orient='records', indent=4)
    print(f"Raw mutation data ({len(raw_export)} records) exported to data/raw_mutations.json")

    # Stats for US-06: Plotly Dashboard
    # For stats, we also use the full mapped df to show global distribution
    gene_dist = df.groupby('hugo_symbol').size().sort_values(ascending=False).head(20).to_dict()
    type_dist = df.groupby('variant_classification').size().to_dict()
    
    stats = {
        "gene_distribution": gene_dist,
        "type_distribution": type_dist,
        "total_samples": total_samples,
        "total_mutations": len(data)
    }
    
    with open('data/mutation_stats.json', 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=4)
    print("Analytical stats exported to data/mutation_stats.json")

def main():
    raw_data, total_samples = fetch_mutations()
    if raw_data:
        # Get all unique Gene IDs to fetch symbols once
        temp_df = pd.DataFrame(raw_data)
        unique_ids = temp_df['entrezGeneId'].unique().tolist()
        gene_map = get_gene_mapping(unique_ids)
        process_and_store(raw_data, gene_map, total_samples)

if __name__ == "__main__":
    main()
