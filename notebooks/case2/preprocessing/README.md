# Case study 2 – preprocessing from raw reads

Attempt to rebuild the expression matrix of case study 2 from the raw sequencing reads (He *et al.*, GEO GSE147424),
to compare it with the CellChat `data_humanSkin` dataset used by the authors. Independent from notebooks 03–07
(which start from CellChat's matrix). The notebooks are in Python (scanpy) and are shared without outputs.

## Pipeline

| Step | File | What it does |
| ---- | ---- | ------------ |
| 1 | `scripts/case2/preprocessing/download_SRR11396173.sh` | Download FASTQ from ENA (SRR11396173 = GSM4430473 = S15_NL, the only small sample in the CellChat subset) |
| 2 | `scripts/case2/preprocessing/run_kb.sh` | kallisto\|bustools quantification (10x v2), run in WSL |
| 3 | `scripts/case2/preprocessing/export_cellchat_S15.R` | Export sample S15 of the CellChat `.rda` as plain files readable from Python |
| 4 | `01_qc_clustering_SRR11396171.ipynb` | QC, normalisation, clustering of one skin sample |
| 5 | `02_qc_clustering_combined.ipynb` | Two samples (SRR11396171 + SRR11396174), Harmony integration |
| 6 | `03_compare_with_cellchat_SRR11396173.ipynb` | Our S15_NL vs. CellChat's cells: barcode overlap, expression, labels, genes of the mosaic |

Run all scripts from the root of the repository.

## Where the files are

Everything big goes to `data/case2/preprocessing/` (not tracked by Git): `fastq/`, `kb_out/`, `write/` (`.h5ad`),
`cache/`. The notebooks start by moving to that folder. The environment is described in `envs/environment_preprocessing.yml`.
