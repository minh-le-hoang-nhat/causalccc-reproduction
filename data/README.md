# Data

Data files are **not tracked by Git** (size, licensing); only this README is. It records where every dataset comes from
and how to obtain it, so that `data/` can be rebuilt from scratch.

**Never modify a file in this folder.** Notebooks read from `data/` and write to `results/`.

## Rebuilding this folder

| Dataset      | Command                                                     | Size                                         |
| ------------ | ----------------------------------------------------------- | -------------------------------------------- |
| Case study 1 | `bash scripts/download_case1.sh`                            | about 770 MB                                 |
|              | `bash scripts/download_case1.sh --with-geo`                 | also GEO raw data, about 860 MB in total     |
| Case study 2 | `bash scripts/case2/download_case2.sh`                            | about 5 MB                                   |

Run the scripts from the root of the repository. They can be run again at any time: complete files are skipped,
interrupted downloads resume, and files that were modified since the last run are reported. Each download is recorded
in `download_log.tsv` (date, file, URL) and all checksums in `checksums.md5` (case study 1) or `checksums_case2.md5`
(case study 2). Check by hand with `cd data && md5sum -c checksums.md5`.

## Layout

```
data/
├── README.md                                    (this file, tracked by Git)
├── download_log.tsv                             when and where each file was downloaded
├── checksums.md5                                md5 checksums of the case study 1 files
├── checksums_case2.md5                          md5 checksums of the case study 2 files
├── case1/
│   ├── causalccc_server_eaP9Br0llvzxN07zZZEC/   files of the authors' web server job
│   ├── nichenet_seurat/seuratObj.rds            processed Medaglia data (Seurat object)
│   └── geo_GSE104054/                           raw data (optional, --with-geo)
├── case2/
│   ├── causalccc_server_w0788tM1wzXGHLzNaKKa/   files of the authors' web server job
│   ├── cellchat/data_humanSkin_CellChat.rda     CellChat example dataset (copied by hand, see below)
│   └── preprocessing/                           raw reads (fastq/), kallisto|bustools output (kb_out/), write/, cache/ – see notebooks/case2/preprocessing/README.md
└── nichenet_v2/                                 NicheNet v2 mouse prior-knowledge networks
```

## Sources – case study 1

CD4⁺ (sender) → CD8⁺ (receiver) T cells, T-cell area of inguinal lymph nodes, steady state and 72 h after LCMV infection
(Medaglia *et al.*, *Science* 2017, [doi:10.1126/science.aao4277](https://doi.org/10.1126/science.aao4277)).
L-R pairs inferred with NicheNet. Figure 2 of Dupuis *et al.* 2025.

| Folder                                          | Source                                                                                                                                          | Content                                                                                                                                                                         | Size   |
| ----------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ |
| `case1/causalccc_server_eaP9Br0llvzxN07zZZEC/`  | CausalCCC web server job, linked as "NicheNet CausalCCC input files" in the [tutorial](https://miic.curie.fr/tutorial_causalCCC.php#datasets) | Authors' input mosaic table, L-R pairs, state order, MIIC settings, run log and result network. Job run on 2025-01-21.                                                          | 4.8 MB |
| `case1/nichenet_seurat/`                        | Zenodo [10.5281/zenodo.3531889](https://doi.org/10.5281/zenodo.3531889) (NicheNet demo data)                                                    | Seurat object of the Medaglia data after Seurat alignment, steady state + LCMV                                                                                                  | 512 MB |
| `nichenet_v2/`                                  | Zenodo [10.5281/zenodo.7074291](https://doi.org/10.5281/zenodo.7074291) (NicheNet v2)                                                           | Mouse L-R network, ligand–target matrix, weighted networks                                                                                                                      | 257 MB |
| `case1/geo_GSE104054/` (optional)               | GEO [GSE104054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104054)                                                                   | `GSE104054_RAW.tar`: raw MARS-seq UMI tables (149 amplification batches). `GSE104054_metadata.txt.gz`: per-well metadata (57,216 wells; lymph node, spleen and tumour; naive and LCMV). Superset of case study 1. | 89 MB  |

Checksums of the Zenodo files are verified against the values published by Zenodo. GEO and the CausalCCC server publish
no checksums: ours are recorded at first download.

### Main files of the authors' job

| File                              | What it is                                                                                                                                 |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `input.tsv` (from `input.zip`)    | Mosaic table: 4,207 cells (2,562 CD4⁺ senders, 1,645 CD8⁺ receivers) × 144 genes, raw UMI counts, `NA` outside each cell type's own genes |
| `clean_interact_edges.tsv`        | The 18 NicheNet L-R pairs (10 ligands, 8 receptors), no header, no score                                                                   |
| `clean_state_order.tsv`           | Type (continuous/discrete) and group (ligand, receptor, sender gene, receiver gene) of each variable                                       |
| `script.sh`, `whole_run.log`      | Exact MIIC settings used by the server, and the log of the run                                                                             |
| `edgesList.miic.summary.txt`      | Result network (364 retained edges, `type = P`). Edited on 2025-04-24.                                                                     |
| `edgesList.miic.summary.txt_save` | Result network as produced by the 2025-01-21 run. Same edges; only the values of the 18 L-R links differ from the edited version.         |

### Notes and open questions

- The tutorial's "NicheNet dataset" link points to the NicheNet v1 model files (Zenodo 3260758, human), not to the
  expression data. That `seuratObj.rds` is the origin of the authors' input is our hypothesis, checked in notebook 03.
- The mosaic has no condition column, although the paper mentions steady state and infection: both are probably pooled.
- The NicheNet version and parameters used to obtain the 18 L-R pairs are not stated.
- The CausalCCC server is not an archive: these files could disappear, hence the local copy and its checksums.

## Sources – case study 2

Fibroblasts, dendritic cells and T cells in atopic dermatitis skin (He *et al.*, *J Allergy Clin Immunol* 2020,
[doi:10.1016/j.jaci.2020.01.042](https://doi.org/10.1016/j.jaci.2020.01.042)). L-R pairs inferred with CellChat
(Jin *et al.*, *Nat Commun* 2021). Figure 3 and Supplementary Figure S1 of Dupuis *et al.* 2025.

| Folder                                         | Source                                                                                         | Content                                                                                                      | Size   |
| ---------------------------------------------- | ---------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ | ------ |
| `case2/causalccc_server_w0788tM1wzXGHLzNaKKa/` | CausalCCC web server job, linked as "CellChat CausalCCC input files" in the [tutorial](https://miic.curie.fr/tutorial_causalCCC.php#datasets) | Authors' input mosaic table, L-R pairs, state order, MIIC settings, run log and result network. Job run on 2025-01-21. | 5 MB   |

### Main files of the authors' job

| File                           | What it is                                                                                                                       |
| ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| `input.tsv` (from `input.zip`) | Mosaic table: 7,563 cells (4,734 fibroblasts, 2,200 T cells, 629 dendritic cells) × 238 variables (232 genes, raw counts, and 6 metadata) |
| `clean_interact_edges.tsv`     | The 43 CellChat L-R pairs, in all six directions between the three cell types, no header, no score                              |
| `clean_state_order.tsv`        | Type and group (ligand, receptor, genes of each cell type, metadata) of each variable                                            |
| `script.sh`, `whole_run.log`   | Exact MIIC settings used by the server (as in case study 1, but `--n_threads 1`), and the log of the run (about 3 hours)        |
| `edgesList.miic.summary.txt`   | Result network (715 retained links out of 8,458 pairs examined). Unlike case study 1, there is no `_save` version.               |

### Notes and open questions

- `input.tsv` is identical (same md5) to the copy found online at the start of the project (`../data/case_study2/input.tsv`
  outside this repository).
- The 7,563 cells are exactly the cells of CellChat's example dataset `data_humanSkin` (`data_humanSkin_CellChat.rda`,
  figshare), which contains these three cell types from He *et al.* The figshare link refused scripted downloads on
  2026-10-05, so it is not downloaded by the script: download it by hand from the
  [CellChat tutorial](https://github.com/jinworks/CellChat) and put it in `data/case2/cellchat/`
  (md5 of our copy: `9d31f4c88be64912b937b5d6a2cb6c2d`). Its matrix is log-normalised; notebook 05 recovers the raw
  counts exactly and shows that they are the authors' counts.
- 14 genes appear twice in the same cell type as identical columns (notebook 03, section 3).

## Other datasets (optional)

Not downloaded by the scripts.

| Dataset                        | Used for        | Source                                                                    | How to download                                                                            | Size   |
| ------------------------------ | --------------- | ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------ |
| CausalCCC demo (Seurat object) | Getting started | [CausalCCC tutorial](https://miic.curie.fr/tutorial_causalCCC.php)        | `wget -c -P data/demo https://miic.curie.fr/datasets/causalCCC_demo_Seurat.RData`          | 274 MB |
