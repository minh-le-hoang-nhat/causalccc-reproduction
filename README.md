# Reproducing and extending CausalCCC

Computational project of the course **Computational Systems Biology of Cancer** (Institut Curie / ENS-PSL), 5–16 October 2026.

[CausalCCC](https://miic.curie.fr/causalCCC.php) reconstructs intracellular causal pathways upstream and downstream of
ligand–receptor interactions between cell types, by combining cell–cell communication (CCC) methods with the causal
network reconstruction method MIIC (Dupuis *et al.*, *Nucleic Acids Research* 2025, [doi:10.1093/nar/gkaf404](https://doi.org/10.1093/nar/gkaf404)).

## Goals

- [ ] **Case study 1** – CD4⁺ (sender) → CD8⁺ (receiver) T cells in inguinal lymph nodes, with NicheNet ligand–receptor pairs
      (data: Medaglia *et al.*, *Science* 2017). Notebooks 01 (authors' files, paper's claims), 02 (MIIC rerun), 03 (CausalCCC
      tutorial, default parameters) and 04 (tutorial, updated parameters).
- [ ] **Case study 2** – fibroblasts, dendritic cells and T cells in atopic dermatitis, with CellChat ligand–receptor pairs
      (data: He *et al.*, *J Allergy Clin Immunol* 2020). Notebooks 03 (authors' files, paper's claims), 04 (MIIC rerun), 05 (rebuilding the input), 06 (MIIC on the rebuilt input), 07 (why our results differ), 08 (reproducibility limits) and 09 (effect of the duplicated genes).
- [ ] Assess the preprocessing of the input data used in the case studies.
- [ ] Extensions of the published analyses (to be defined with the teaching team).
- [ ] Analysis of a dataset shared with the BoNesis team, and comparison of the two approaches.

## Repository structure

```
├── README.md
├── .Rprofile            # pins the CRAN snapshot used by this project
├── docs/INSTALL.md      # how to build the computing environment
├── envs/                # environment definitions and recorded package versions
├── data/                # not tracked by Git – see data/README.md; subfolders case1/, case2/
├── notebooks/           # one folder per case study, notebooks numbered in execution order
│   ├── case1/           # case study 1: 01_* to 04_*
│   └── case2/           # case study 2: 03_* to 07_*, and preprocessing/ (raw reads -> matrix, Python)
├── scripts/
│   ├── helpers.R, run_miic.R   # shared by all case studies
│   ├── case1/           # download_case1.sh, case1_helpers.R (functions of notebooks 03 and 04)
│   └── case2/           # download_case2.sh, MIIC job, preprocessing/ (kallisto|bustools)
└── results/
    ├── case1/           # one folder per notebook, same name
    ├── case2/           # one folder per notebook, same name
    └── figures/
```

## Getting started

1. Build the environment: see [docs/INSTALL.md](docs/INSTALL.md).
2. Download the data: see [data/README.md](data/README.md).
3. Run the notebooks in numerical order.

## Team

Minh Le and Jeanne Leclerc (ENS-PSL).

## Acknowledgements

We thank the teaching team (Loredana Martignetti, Denis Thieffry, Alexandre Lanaud) and the developers of CausalCCC and MIIC
at Institut Curie (Hervé Isambert, Franck Simon) for their guidance.
