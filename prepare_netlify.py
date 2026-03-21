import os
import json
import shutil

def prepare_netlify():
    source_dir = os.getcwd()
    target_dir = os.path.join(source_dir, "GEMAMXST")
    
    # 1. Limpiar/Crear directorio destino
    if os.path.exists(target_dir):
        print(f"Limpiando directorio existente: {target_dir}")
        shutil.rmtree(target_dir)
    os.makedirs(target_dir)
    
    # 2. Archivos raíz indispensables
    core_files = [
        "index.html", "styles.css", 
        "app.js", "dashboard.js", "analysis.js", 
        "structural.js", "docking.js", "report.js",
        "banner.png"
    ]
    
    for f in core_files:
        src = os.path.join(source_dir, f)
        if os.path.exists(src):
            shutil.copy2(src, target_dir)
            print(f"Copiado: {f}")
        else:
            print(f"ADVERTENCIA: No se encontró {f}")

    # 3. Carpeta de datos (JSONs)
    data_src = os.path.join(source_dir, "data")
    data_target = os.path.join(target_dir, "data")
    os.makedirs(data_target)
    
    essential_jsons = [
        "top_variants.json", "raw_mutations.json", 
        "structural_metrics.json", "manifest_structures.json", 
        "mutation_stats.json"
    ]
    
    for j in essential_jsons:
        src = os.path.join(data_src, j)
        if os.path.exists(src):
            shutil.copy2(src, data_target)
            print(f"Copiado: data/{j}")

    # 3.1 Copiar FASTA (secuencias)
    fasta_src = os.path.join(data_src, "fasta")
    if os.path.exists(fasta_src):
        fasta_target = os.path.join(data_target, "fasta")
        if os.path.exists(fasta_target): shutil.rmtree(fasta_target)
        shutil.copytree(fasta_src, fasta_target)
        print("Copiada: data/fasta (secuencias)")

    # 4. Estructuras filtradas
    manifest_path = os.path.join(data_src, "manifest_structures.json")
    if os.path.exists(manifest_path):
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
        
        # Copiar solo archivos referenciados en el manifest
        struct_count = 0
        for rel_path in manifest.values():
            if rel_path and isinstance(rel_path, str):
                # rel_path suele ser "data/structures/..."
                src_path = os.path.join(source_dir, rel_path)
                dest_path = os.path.join(target_dir, rel_path)
                
                if os.path.exists(src_path):
                    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
                    shutil.copy2(src_path, dest_path)
                    struct_count += 1
                else:
                    print(f"ADVERTENCIA: Estructura no encontrada: {rel_path}")
        
        print(f"Total de estructuras copiadas: {struct_count}")

    # 5. Resumen
    total_size = sum(os.path.getsize(os.path.join(dp, f)) for dp, dn, filenames in os.walk(target_dir) for f in filenames)
    print(f"\n--- Resumen GEMAMXST ---")
    print(f"Directorio creado en: {target_dir}")
    print(f"Tamaño total: {total_size/1024/1024:.2f} MB")
    print("Listo para subir a Netlify.")

if __name__ == "__main__":
    prepare_netlify()
