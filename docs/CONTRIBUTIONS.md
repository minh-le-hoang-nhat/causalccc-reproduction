# Contributions

Who did what in this project. Choices were discussed between the two of us, and notebooks list both names as
authors; this file records who carried out each part. Details and results are in the notebooks themselves (each ends
with a summary).

| Part                                             | Minh Le | Jeanne Leclerc | Notes                                                                                        |
| ------------------------------------------------ | :-----: | :------------: | -------------------------------------------------------------------------------------------- |
| Repository, computing environment, data download |    ●    |       ○        | Jeanne: case study 2 download, preprocessing from raw reads                                  |
| Case study 1 – lymph node, NicheNet              |    ●    |       ○        | Jeanne: performance compared with the web server                                             |
| Case study 2 – atopic dermatitis, CellChat       |    ○    |       ●        | Follows the approach of case study 1 (authors' files → MIIC rerun → rebuilt input)           |
| Running MIIC on several threads                  |    ○    |       ●        | Minh: stack size limit, used in the thymus runs; Jeanne: `miic` patches, benchmarks          |
| CausalCCC tutorial                               |    ○    |       ●        | One pipeline combining the workflows of both case studies                                    |
| Thymus – early niche                             |    ●    |                |                                                                                              |
| Thymus – other niches, CellPhoneDB re-run        |    ○    |       ●        | XCL1 niche built from the early-niche notebooks                                              |

● led · ○ contributed or built on

*Status (2026-10-18): Minh Le's part is complete. Jeanne Leclerc's part is in progress; her section below is a draft
to be completed by her.*

## Minh Le

### Project structure and computing environment

- Repository layout, conventions (one folder per case study, notebooks numbered in execution order, one results
  folder per notebook) and `README.md`.
- Reproducible environment: pinned conda environment (`envs/environment.yml`), CRAN snapshot of 2026-10-03
  (`.Rprofile`), R install script with `miic` (branch `causalccc`), `causalCCC`, LIANA and AUCell at fixed versions
  (`envs/install_R_packages.R`), recorded package versions (`envs/*_versions.txt`), `docs/INSTALL.md` (sections 0–3).
- Data download scripts and data sources: `scripts/case1/download_case1.sh`, `scripts/thymus/download_thymus.sh`,
  `data/README.md`.
- Integration of both parts of the work into one repository (folder reorganisation, merges).

### Case study 1 – CD4⁺ → CD8⁺ T cells, lymph node (`notebooks/case1/01`–`04`)

- **01** – The files the authors uploaded to the CausalCCC server, and the paper's statements about Figure 2.
- **02** – MIIC rerun on the authors' input with the server's settings: the network is reproduced exactly (same
  decision for all 3,798 gene pairs, same 364 links, same orientations and 13 genuine causations).
- **03** – The CausalCCC tutorial with default parameters: 61 of the authors' 144 variables, 107 of their 364 links,
  12 of the 17 links described in the paper. The differences come from the gene selection.
- **04** – Three settings explain the authors' gene selection (ribosomal protein genes removed, 60 genes per side,
  sender genes added to receiver genes): all 144 variables, 347 of 364 links, all 17 paper links with the authors'
  orientation.
- Points raised for the authors: `proba` column order in the `miic` documentation, `NA` counted as a value when
  `causalCCC.state_order()` decides whether a variable is discrete, undocumented L–R link strengths in the edited edge
  list, duplicated L–R rows in the `causalccc` branch of `miic`.

### Extension – early thymic niche (`notebooks/thymus_early/01`–`10`)

CausalCCC on the prenatal human thymus at 7–11 post-conception weeks (Park *et al.*, *Science* 2020), to test the
paper's ligand–receptor claims in both directions.

- **01–03** – Data inventory (which files hold the counts, labels and metadata, and which are unreliable), choice of
  cell populations, orientation of the authors' CellPhoneDB pairs into ligand → receptor.
- **04–06** – Run A, stroma → early lymphoid cells (Notch ligands JAG2 and DLL4 → NOTCH1). Sequencing depth found to
  drive part of the network and added as a variable.
- **07–09** – Run B, early lymphoid cells → stroma (LTB, RANKL, CD40LG), and comparison of both runs.
- **10** – A per-cell Notch activity score (CollecTRI RBPJ regulon, AUCell). It fails its positive control, so it is not
  used in a network.
- Running MIIC on 8 threads (with an unlimited stack size) gives the one-thread network on case study 1 in a quarter of
  the time (notebook 05).
- Conclusion: the expression side of every claim is reproduced in both directions; the signalling side is not
  resolved, because receptor mRNA does not reflect receptor activity, and sender and receiver cells are unpaired.
  Sequencing depth and the confounding of population with donor in the stroma limit the interpretation.
- Fix in the shared script `scripts/run_miic.R` (state order levels read as text).

## Jeanne Leclerc

*Draft from the repository history, to be completed and corrected by Jeanne.*

- **Case study 2** – fibroblasts, dendritic cells and T cells in atopic dermatitis: `notebooks/case2/03`–`10`,
  preprocessing from raw reads (`notebooks/case2/preprocessing/`, `scripts/case2/`), network plotting functions
  (`scripts/helpers.R`).
- **Case study 1** – `notebooks/case1/05` (performance compared with the web server).
- **Faster `miic`** – patches that make multithreaded runs crash-free and deterministic (`envs/patches/miic/`),
  benchmarks (`scripts/benchmarks/`, `results/benchmarks/`), `docs/INSTALL.md` section 4; shared script
  `scripts/run_miic.R`.
- **Tutorial** – `notebooks/tutorial/causalccc_tutorial.qmd`: CausalCCC step by step on case study 2.
- **Thymus** – CellPhoneDB re-run on the atlas, overall and by stage (`scripts/thymus/01`–`05`, `results/thymus/`);
  XCL1 niche (`notebooks/thymus_case2_xcl1/`).
- *To add: …*
