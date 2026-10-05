# Data

Data files are **not tracked by Git** (size, licensing). This file records where every dataset comes from and how to
obtain it, so that `data/` can be rebuilt from scratch.

- `data/raw/` – files exactly as downloaded, never modified.
- `data/processed/` – files produced by the notebooks and scripts.

| Dataset | Used for | Source | How to download | Notes |
|---|---|---|---|---|
| CausalCCC demo (Seurat object) | Getting started | https://miic.curie.fr/tutorial_causalCCC.php | `wget -c -P data/raw https://miic.curie.fr/datasets/causalCCC_demo_Seurat.RData` | 274 MB |
| Medaglia *et al.* 2017, lymph node T cells | Case study 1 | TBD | TBD | |
| He *et al.* 2020, atopic dermatitis skin | Case study 2 | TBD | TBD | |
