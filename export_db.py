import pymongo
import json
import gzip
import os

def export_db():
    try:
        client = pymongo.MongoClient("mongodb://localhost:27017/", serverSelectionTimeoutMS=5000)
        db = client["gen_mama_db"]
        
        # Obtener lista de colecciones existentes
        existing_collections = db.list_collection_names()
        print(f"Colecciones encontradas: {existing_collections}")
        
        os.makedirs("export", exist_ok=True)
        
        for coll_name in existing_collections:
            print(f"Exportando colección: {coll_name}...")
            # Excluimos _id para que sea puramente portátil
            data = list(db[coll_name].find({}, {"_id": 0}))
            
            output_path = f"export/{coll_name}.json.gz"
            with gzip.open(output_path, "wt", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            
            size = os.path.getsize(output_path) / 1024
            print(f"Finalizado: {output_path} ({size:.2f} KB)")
            
        print("\nRespaldo de base de datos completado exitosamente en la carpeta 'export/'.")
        
    except Exception as e:
        print(f"Error durante el export: {e}")

if __name__ == "__main__":
    export_db()
