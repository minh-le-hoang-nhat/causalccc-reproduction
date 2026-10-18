# Reproducing and extending CausalCCC

Computational project of the course **Computational Systems Biology of Cancer** (Institut Curie / ENS-PSL), 5–16 October 2026.

[CausalCCC](https://miic.curie.fr/causalCCC.php) reconstructs intracellular causal pathways upstream and downstream of
ligand–receptor interactions between cell types, by combining cell–cell communication (CCC) methods with the causal
network reconstruction method MIIC (Dupuis *et al.*, *Nucleic Acids Research* 2025, [doi:10.1093/nar/gkaf404](https://doi.org/10.1093/nar/gkaf404)).

## Where to start

Every R notebook is rendered as a self-contained `.html` file next to its `.qmd` source; Python notebooks (`.ipynb`)
are stored with their outputs and display directly on GitHub. GitHub shows the *source* of `.html` files: open them
from the rendered site, <https://minh-le-hoang-nhat.github.io/causalccc-reproduction/>, or download them.

| Notebook | Question | Answer in one line |
| -------- | -------- | ------------------ |
| [Case 1, 02](notebooks/case1/02_rerun_miic_authors_input.html) | Same input, same settings: same network? | Yes: the same 364 links, orientations and 13 genuine causations as the web server. |
| [Case 1, 04](notebooks/case1/04_tutorial_updated_params.html) | Does the published code rebuild Figure 2? | With three settings inferred from the authors' input: 347 of their 364 links, and all 17 links described in the paper. |
| [Case 2, 07](notebooks/case2/07_differences.html) | Why do our results differ from the authors'? | Their gene selection is reproduced exactly with non-default settings; 14 duplicated columns make some orientations contradictory. |
| [Case 2, 08](notebooks/case2/08_reproducibility_limits.html) | Where is the pipeline not reproducible? | Twelve points, from missing raw counts to undocumented settings and no recorded random seed. |
| [Case 2, 10](notebooks/case2/10_performance.html) | Can CausalCCC be faster without changing its results? | Three patches to `miic` make 8 threads usable: network reconstruction about 3× faster, same decisions. |
| [Tutorial](notebooks/tutorial/causalccc_tutorial.html) | How to run CausalCCC on one's own data? | Step by step on case study 2, with the pitfalls and the measured effect of each setting. |
| [Thymus, 09](notebooks/thymus_early/09_runB_interpretation.html) | Do the early thymic niche claims of Park *et al.* hold? | The ligand (expression) side of every claim is reproduced; the receptor → target side is not resolved. |
| [Thymus XCL1, 05](notebooks/thymus_case2_xcl1/05_interpretation.html) | Does the XCL1 → XCR1 relay to DC1 hold? | Partly: NK cells are the main XCL1 source, and DC1 activation is not linked to XCR1 in the network. |
| [Thymus XCL1 by age, 05](notebooks/thymus_case3_xcl1_age/05_interpretation.html) | Does the CD8αα(I) source of XCL1 change with age? | Yes: XCL1 drops at birth, and at both ages it goes with an activation program (TNFRSF9 / 4-1BB) of CD8αα(I) cells. |

## Goals

- [x] **Case study 1** – CD4⁺ (sender) → CD8⁺ (receiver) T cells in inguinal lymph nodes, with NicheNet ligand–receptor pairs
      (data: Medaglia *et al.*, *Science* 2017). `notebooks/case1/`: 01 (authors' files, paper's claims), 02 (MIIC
      rerun), 03 (CausalCCC tutorial, default parameters), 04 (tutorial, updated parameters) and 05 (patched `miic`
      against the web server).
- [x] **Case study 2** – fibroblasts, dendritic cells and T cells in atopic dermatitis, with CellChat ligand–receptor pairs
      (data: He *et al.*, *J Allergy Clin Immunol* 2020). `notebooks/case2/`: 03 (authors' files, paper's claims), 04
      (MIIC rerun), 05 (rebuilding the input), 06 (MIIC on the rebuilt input), 07 (why our results differ), 08
      (reproducibility limits), 09 (effect of the duplicated genes) and 10 (making CausalCCC faster without changing its
      results: patched `miic`, `envs/patches/miic/`, benchmarks in `scripts/benchmarks/` and `results/benchmarks/`).
- [x] **Tutorial** – `notebooks/tutorial/causalccc_tutorial.qmd`: CausalCCC step by step on case study 2, from raw counts to the
      network, with the pitfalls of each step and the measured effect of the settings of each function and of MIIC.
- [ ] Assess the preprocessing of the input data used in the case studies (case study 2: `notebooks/case2/preprocessing/`,
      from the raw reads; notebook 08, section 11).
- [ ] Extensions of the published analyses, on the human thymus atlas of Park *et al.*, *Science* 2020
      ([doi:10.1126/science.aay3224](https://doi.org/10.1126/science.aay3224)), with the CellPhoneDB ligand–receptor
      pairs of the authors.
  - [x] **CellPhoneDB re-run** on the atlas, overall and by stage, and reproduction of fig. S11
        (`scripts/thymus/01_*.py` to `05_*.py`, `results/thymus/`).
  - [x] **Early thymic niche** – stromal cells (cTEC, mTEC(II), Fb_1, endothelium) and early lymphoid cells (ETP, DN,
        ILC3, NK, γδT) of the prenatal human thymus at 7–11 post-conception weeks, with the CellPhoneDB ligand–receptor
        pairs of the authors (data: Park *et al.*, *Science* 2020, [doi:10.1126/science.aay3224](https://doi.org/10.1126/science.aay3224)).
        `notebooks/thymus_early/`: 01 (data inventory), 02 (cell populations), 03 (candidate ligand–receptor pairs),
        04–06 (run A, stroma → lymphoid: input, CausalCCC, interpretation of the Notch claim), 07–09 (run B,
        lymphoid → stroma: links, CausalCCC, interpretation of the LTB, RANKL and CD40LG claims, comparison of both runs),
        10 (a Notch activity score per cell from the CollecTRI RBPJ regulon with AUCell: it fails its positive control,
        ETP and DN(early) scoring lowest of the lymphoid populations, so it is not used in a network).
        Notebooks 01–04 and 07 are Python (`.ipynb`), 05, 06, 08, 09 and 10 R (Quarto `.qmd`).
  - [ ] CellPhoneDB re-run on the thymus atlas, overall and by stage (`scripts/thymus/01_*.py` to `05_*.py`,
        `results/thymus/`).
  - [ ] Medullary niche of the same atlas (Fig. 4 of Park *et al.*).

  **Limits of the thymus analysis.** Sender and receiver cells are unpaired: the ligand → receptor links come from
  CellPhoneDB and are imposed in the network, and the causal inference is within each group of cells. In the stroma,
  cell population and donor are confounded (cTEC and mTEC(II) come almost only from two donors), and sequencing depth
  still structures part of the networks (thymus_early notebooks 06 and 09).
- [ ] Analysis of a dataset shared with the BoNesis team, and comparison of the two approaches.

## Repository structure

```
├── README.md
├── .Rprofile            # pins the CRAN snapshot used by this project
├── docs/                # INSTALL.md (how to build the computing environment), text of the paper's case study 1
├── envs/                # environment definitions, recorded package versions, patches/miic/ (faster miic)
├── data/                # not tracked by Git – see data/README.md; subfolders case1/, case2/, thymus/
├── notebooks/           # one folder per analysis, notebooks numbered in execution order
│   ├── case1/           # case study 1: 01_* to 05_*
│   ├── case2/           # case study 2: 03_* to 10_*, and preprocessing/ (raw reads -> matrix, Python)
│   ├── tutorial/        # CausalCCC step by step: preprocessing, pitfalls, choice of settings
│   └── thymus_early/    # early thymic niche: 01_* to 10_* (Python .ipynb and R Quarto .qmd)
├── scripts/
│   ├── helpers.R, run_miic.R   # shared by all case studies
│   ├── case1/           # download_case1.sh, case1_helpers.R (functions of notebooks 03 and 04)
│   ├── case2/           # download_case2.sh, MIIC job, preprocessing/ (kallisto|bustools)
│   ├── tutorial/        # tutorial_helpers.R (functions of the tutorial)
│   ├── benchmarks/      # building the miic variants and timing them (case1 notebook 05, case2 notebook 10)
│   └── thymus/          # download_thymus.sh; CellPhoneDB re-run (01_*.py to 05_*.py)
└── results/
    ├── case1/, case2/, tutorial/, thymus_early/, thymus_case2_xcl1/, thymus_case3_xcl1_age/
    │                    # one folder per notebook, same name (exported count matrices *.mtx.gz not tracked by Git);
    │                    # tutorial/server_<code>/: local MIIC rerun of a CausalCCC web-server job (tutorial, section 8b)
    ├── benchmarks/      # timings of the miic variants, see results/benchmarks/README.md
    ├── thymus/          # CellPhoneDB re-run, fig. S11, first_results.html (summary page of these results)
    └── figures/
```

Python notebooks (`.ipynb`) prepare the cells and the ligand–receptor links; R notebooks (Quarto `.qmd`) run CausalCCC
and MIIC and read the networks.

## Getting started

1. Build the environment: see [docs/INSTALL.md](docs/INSTALL.md).
2. Download the data: see [data/README.md](data/README.md).
3. Run the notebooks in numerical order, from the root of the repository.

## Team

Minh Le and Jeanne Leclerc (ENS-PSL). Who did what: [docs/CONTRIBUTIONS.md](docs/CONTRIBUTIONS.md).

## Acknowledgements

We thank the teaching team (Loredana Martignetti, Denis Thieffry, Alexandre Lanaud) and the developers of CausalCCC and MIIC
at Institut Curie (Hervé Isambert, Franck Simon) for their guidance.
