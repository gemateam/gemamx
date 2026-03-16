import os
import json
import subprocess

# --- CONFIGURATION ---
# Check if running in Google Colab
IS_COLAB = os.path.exists("/content")
BASE_PATH = "/content" if IS_COLAB else "."

DATA_DIR = os.path.join(BASE_PATH, "data")
FASTA_DIR = os.path.join(DATA_DIR, "fasta")
STRUCTURES_DIR = os.path.join(DATA_DIR, "structures")
TOP_VARIANTS_FILE = os.path.join(DATA_DIR, "top_variants.json")

def setup_boltz():
    """Installs Boltz-1 in the environment."""
    print("Installing Boltz-1...")
    subprocess.run(["pip", "install", "boltz"], check=True)
    print("Boltz-1 installed.")

def predict_variant(hugo_symbol, protein_change):
    """Runs Boltz prediction for WT and MUT sequences."""
    wt_fasta = os.path.join(FASTA_DIR, f"{hugo_symbol}_WT.fasta")
    mut_fasta = os.path.join(FASTA_DIR, f"{hugo_symbol}_{protein_change}.fasta")
    
    for fasta_path, tag in [(wt_fasta, "WT"), (mut_fasta, "MUT")]:
        if not os.path.exists(fasta_path):
            print(f"Skipping {tag}: {fasta_path} not found.")
            continue
            
        output_name = f"{hugo_symbol}_{protein_change}_{tag}"
        print(f"Predicting structure for {output_name}...")
        
        # Boltz-1 CLI command (simplified for demonstration)
        # Note: Actual boltz command may vary based on version and requirements
        try:
            subprocess.run([
                "boltz", "predict", 
                fasta_path,
                "--out_dir", STRUCTURES_DIR,
                "--name", output_name
            ], check=True)
            print(f"Successfully predicted {output_name}")
        except Exception as e:
            print(f"Error predicting {output_name}: {e}")

def main():
    if not os.path.exists(STRUCTURES_DIR):
        os.makedirs(STRUCTURES_DIR)
        
    if not os.path.exists(TOP_VARIANTS_FILE):
        print(f"Error: {TOP_VARIANTS_FILE} not found.")
        return

    with open(TOP_VARIANTS_FILE, 'r') as f:
        variants = json.load(f)

    # setup_boltz() # Uncomment in Colab if running directly
    
    # Check if boltz is in path, if not install it
    try:
        subprocess.run(["boltz", "--version"], check=True, capture_output=True)
    except:
        setup_boltz()

    for v in variants:
        gene = v.get('hugo_symbol', v.get('Gene')) # Handle both formats
        change = v.get('protein_change', v.get('AA_Change'))
        if gene and change:
            predict_variant(gene, change)

if __name__ == "__main__":
    main()
