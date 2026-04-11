# GEMAMX GROMACS

Modulo para extender el pipeline de docking con una capa posterior de
dinamica molecular en Google Colab.

## Objetivo

Usar los resultados de AutoDock Vina como filtro rapido y ejecutar GROMACS
solo sobre los mejores pares variante-farmaco. La salida queda normalizada
para la app web de `http://localhost:8000/`.

## Flujo operativo

1. Generar lote desde Vina:
   - `python gromacs/build_gromacs_batch.py`
   - salida: `gromacs/results/gromacs_batch_manifest.json`
2. Ejecutar en Colab:
   - `gromacs/01_prepare_gromacs_batch_from_vina.ipynb`
   - `gromacs/02_run_gromacs_md_colab.ipynb`
   - o piloto corto: `gromacs/03_run_gromacs_pilot_colab.ipynb`
   - o flujo consolidado por par: `gromacs/04_run_gromacs_pair_definitive_colab.ipynb`
3. Guardar resultados en Google Drive:
   - `gromacs/results/gromacs_metrics.json`
4. Ingerir resultados en local:
   - `python gromacs/ingest_gromacs_results.py`
5. Visualizar:
   - abrir `http://localhost:8000/`
   - navegar a `Dinamica Molecular`

## Archivos clave

- `build_gromacs_batch.py`
  Genera el manifiesto canonicamente desde `data/docking/vina_affinities.json`.
- `ingest_gromacs_results.py`
  Normaliza el JSON de Colab, actualiza `data/structural_metrics.json` y deja
  una copia local consumible por la app.
- `gromacs.js`
  Renderiza la nueva pagina web.
- `results/gromacs_batch_manifest.json`
  Lista exacta de pares y rutas de entrada para Colab.
- `results/gromacs_metrics.json`
  Contrato de salida esperado desde Colab.

## Contrato de salida

Cada par debe incluir como minimo:

- `tag`
- `variant_id`
- `drug_name`
- `vina_score_kcalmol`
- `status`

Campos recomendados:

- `mmgbsa_delta_g`
- `interaction_energy`
- `rmsd_mean_nm`
- `rmsd_std_nm`
- `rmsf_mean_nm`
- `hbonds_mean`
- `hbonds_max`
- `run_dir`
- `trajectory_path`
- `plot_paths`
- `notes`

## Notas metodologicas

- La recomendacion base queda en `top_n_per_variant = 5`.
- Para una primera validacion practica se recomienda empezar con el piloto corto:
  - `PIK3CA_H1047R + Alpelisib`
  - `PIK3CA_H1047R + Pictilisib`
  - `ERBB2_V777L + Pictilisib`
  - `ERBB2_V777L + Alpelisib`
  - `TP53_R248Q + Abemaciclib`
  - `TP53_R248Q + Pictilisib`
- El modulo esta pensado para Linux/Colab, no para ejecutar GROMACS pesado en
  Windows.
- La app web no bloquea si aun no hay simulaciones: muestra el manifiesto y el
  estado del pipeline.
