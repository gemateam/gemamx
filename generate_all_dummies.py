import json
import os

base_dir = r"d:\ProyectosIA\GEMAMX\data\structures"
if not os.path.exists(base_dir):
    os.makedirs(base_dir)

# More robust dummy PDB (Small peptide fragment)
pdb_content = """HEADER    DUMMY PROTEIN
MODEL        1
ATOM      1  N   MET A   1      -0.015   1.465   0.033  1.00 20.00           N
ATOM      2  CA  MET A   1       0.142   0.010   0.021  1.00 20.00           C
ATOM      3  C   MET A   1       1.597  -0.344  -0.038  1.00 20.00           C
ATOM      4  O   MET A   1       2.002  -1.488  -0.198  1.00 20.00           O
ATOM      5  CB  MET A   1      -0.623  -0.638  -1.127  1.00 20.00           C
ATOM      6  CG  MET A   1      -0.573  -2.152  -1.139  1.00 20.00           C
ATOM      7  SD  MET A   1      -1.493  -2.894  -2.484  1.00 20.00           S
ATOM      8  CE  MET A   1      -0.297  -2.618  -3.791  1.00 20.00           C
TER       9      MET A   1
ENDMDL
END
"""

with open(r"d:\ProyectosIA\GEMAMX\data\top_variants.json", "r") as f:
    variants = json.load(f)

for v in variants:
    hugo = v["hugo_symbol"]
    change = v["protein_change"]
    
    # WT
    with open(os.path.join(base_dir, f"{hugo}_WT.pdb"), "w") as out:
        out.write(pdb_content)
        
    # MUT
    with open(os.path.join(base_dir, f"{hugo}_{change}_MUT.pdb"), "w") as out:
        out.write(pdb_content)

# Global placeholder
with open(os.path.join(base_dir, "placeholder.pdb"), "w") as out:
    out.write(pdb_content)

print(f"Comprehensive PDB mocks generated.")
