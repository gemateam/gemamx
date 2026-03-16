import json

path = r"d:\ProyectosIA\GEMAMX\01_tcga_exploration.ipynb"
with open(path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

new_cells = [
    {
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 4. Predicción de Estructuras (Boltz-1)\n",
            "\n",
            "Este bloque ejecuta la predicción 3D para la proteína Wild-Type y la Mutada de forma automática."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Ejecutar script de predicción automatizada\n",
            "!python scripts/predict_structures.py"
        ]
    }
]

nb['cells'].extend(new_cells)

with open(path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print("Notebook updated successfully.")
