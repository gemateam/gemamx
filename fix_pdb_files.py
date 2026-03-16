import os
import shutil

base_dir = r"d:\ProyectosIA\GEMAMX\data\structures"
placeholder = os.path.join(base_dir, "placeholder.pdb")

# Standard PDB format placeholder (Glycine)
pdb_content = """HEADER    DUMMY PROTEIN
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLY A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   GLY A   1       1.233   2.356   0.000  1.00  0.00           O
TER       5      GLY A   1
END
"""

files = [
    "PIK3CA_WT.pdb",
    "PIK3CA_H1047R_MUT.pdb",
    "TP53_WT.pdb",
    "TP53_R248Q_MUT.pdb",
    "ERBB2_WT.pdb",
    "ERBB2_V777L_MUT.pdb",
    "BRCA1_WT.pdb",
    "BRCA1_C61G_MUT.pdb"
]

if not os.path.exists(base_dir):
    os.makedirs(base_dir)

# Create placeholder
with open(placeholder, "w") as f:
    f.write(pdb_content)

for f in files:
    target = os.path.join(base_dir, f)
    with open(target, "w") as out:
        out.write(pdb_content)
    print(f"Created {target}")
