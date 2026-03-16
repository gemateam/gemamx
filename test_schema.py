import requests
import json

BASE_URL = "https://www.cbioportal.org/api"
MOLECULAR_PROFILE_ID = "brca_tcga_mutations"

def check_schema():
    url = f"{BASE_URL}/molecular-profiles/{MOLECULAR_PROFILE_ID}/mutations/fetch"
    # Just fetch for one sample
    payload = { "sampleIds": ["TCGA-A1-A0SB-01"] }
    res = requests.post(url, json=payload)
    if res.status_code == 200:
        data = res.json()
        if data:
            print("Keys in one mutation record:", data[0].keys())
            print("Sample record:", json.dumps(data[0], indent=2))
        else:
            print("No data returned for sample.")
    else:
        print(f"Error: {res.status_code}")

if __name__ == "__main__":
    check_schema()
