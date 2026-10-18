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
| Thymus – XCL1 niche, CellPhoneDB re-run          |    ○    |       ●        | XCL1 niche (all ages, then fetal vs postnatal) built from the early-niche notebooks          |

● led · ○ contributed or built on

*Status (2026-10-18): both parts are complete.*

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

### Case study 2 – fibroblasts, dendritic cells and T cells, atopic dermatitis (`notebooks/case2/03`–`10`)

- **03** – The files the authors uploaded to the CausalCCC server, and the paper's statements about Figure 3C (all
  supported within the neighbourhood the figure draws). 14 genes appear twice, as identical columns under two names.
- **04** – MIIC rerun on the authors' input with the server's settings: the network is reproduced. As in case study 1,
  the MIIC step is reproducible from the authors' files.
- **05** – Input rebuilt from CellChat's `data_humanSkin`: cells and counts reproduced exactly (normalisation undone).
  `causalCCC.wrapper()` does not run as published (Seurat 5, data frame input, one sender and one receiver only).
- **06** – MIIC on the rebuilt input: the skeleton of the network (Jaccard 0.85–0.92 per cell type) and the Figure 3C
  neighbourhoods are reproduced; orientations mostly (96 % of the shared links without a duplicated gene), genuine
  causations less so (8 in common out of 14 and 12).
- **07** – The causes of the differences. The authors' gene selection is reproduced exactly with undocumented settings
  (25 genes per run, no cap of 8 ligands or receptors, both metadata). Removing the duplicated columns changes no link,
  but the 111 links of the duplicated genes have contradictory orientations in the published network.
- **08** – The 12 places where the pipeline is not reproducible as published (summary table).
- **09** – Adding the 14 duplicated columns to the rebuilt input gives the authors' network exactly (same links,
  orientations and genuine causations). The duplicates also set the strength of all 43 L–R links (11,579 instead of
  3,050).
- **10** – Making CausalCCC faster without changing its results (see below).
- Preprocessing from raw reads (`notebooks/case2/preprocessing/`, `scripts/case2/`): kallisto|bustools, QC and
  clustering in scanpy, comparison with CellChat's cells. Agreement is only partial, because the preprocessing of the
  CellChat dataset is not described.
- Network plotting functions (`scripts/helpers.R`), used by the later notebooks.

### Faster MIIC (`notebooks/case2/10`, `notebooks/case1/05`)

- Patches to `miic` (`envs/patches/miic/`): no console output from worker threads (stock `miic` crashed with several
  threads), early decision of the confidence cut, exact ties broken by position (results independent of thread timing).
- On 8 threads, the network reconstruction is 3 times faster with the same links, orientations and genuine causations:
  22 instead of 70 minutes on the rebuilt case 2 input, 30 instead of 91 on the authors' input. Only the `confidence`
  values move slightly (up to 0.02).
- `compute_mi_table()` for the gene selection: identical MI values, 21 % faster, a third of the memory.
- Case study 1: every local run gives the server's network; 288 s on 8 threads against 935 s on the server. The patches
  are not faster at equal threads; the gain comes only from using 8 threads.
- Benchmark scripts and results (`scripts/benchmarks/`, `results/benchmarks/`), `docs/INSTALL.md` section 4, shared
  script `scripts/run_miic.R`.

### CausalCCC tutorial (`notebooks/tutorial/causalccc_tutorial.qmd`)

CausalCCC step by step on case study 2, from a raw count matrix to the network. It combines what both case studies
taught: preprocessing, pitfalls of each function, how to choose the settings, how to compare a web server job with a
local run, and a checklist of 22 points.

### Extension – thymus, XCL1 relay to DC1 (`scripts/thymus/`, `notebooks/thymus_case2_xcl1/`, `notebooks/thymus_case3_xcl1_age/`)

- **CellPhoneDB re-run** on the atlas, overall and by developmental stage, following the paper's protocol, and
  reproduction of fig. S11 (`scripts/thymus/01`–`05`, `results/thymus/first_results.html`).
- **XCL1 niche, all ages** (`thymus_case2_xcl1/01`–`05`) – claim: CD8αα(I) and NK cells express XCL1 and recruit XCR1⁺
  DC1, which become activated and up-regulate CCR7. Built on the early-niche pipeline, with sequencing depth as a
  variable (runs with and without depth compared). The claim is partly supported: NK cells are the main XCL1 source;
  XCR1 (resting DC1) and CCR7 (activated DC1) are linked only through the population label; costimulation pairs reach
  the activation program of DC1 more directly than XCL1 → XCR1.
- **XCL1 niche, fetal vs postnatal** (`thymus_case3_xcl1_age/01`–`05`) – CD8αα(I) senders only, one run per age. XCL1
  in CD8αα(I) cells drops at birth (29 % → 12 % of cells) and goes with an activation program (TNFRSF9 / 4-1BB) at both
  ages; the receiver side is the same at both ages; costimulation links to activated DC1 are a fetal feature.
