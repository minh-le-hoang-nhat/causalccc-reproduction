# Data

Data files are **not tracked by Git** (size, licensing); only this README is. It records where every dataset comes from
and how to obtain it, so that `data/` can be rebuilt from scratch.

**Never modify a file in this folder.** Notebooks read from `data/` and write to `results/`.

## Rebuilding this folder

| Dataset      | Command                                                     | Size                                         |
| ------------ | ----------------------------------------------------------- | -------------------------------------------- |
| Case study 1 | `bash scripts/download_case1.sh`                            | about 770 MB                                 |
|              | `bash scripts/download_case1.sh --with-geo`                 | also GEO raw data, about 860 MB in total     |
| Case study 2 | not scripted yet (see [below](#sources--case-study-2))      |                                              |

Run the scripts from the root of the repository. They can be run again at any time: complete files are skipped,
interrupted downloads resume, and files that were modified since the last run are reported. Each download is recorded
in `download_log.tsv` (date, file, URL) and all checksums in `checksums.md5` (check by hand with
`cd data && md5sum -c checksums.md5`).

## Layout

```
data/
├── README.md                                    (this file, tracked by Git)
├── download_log.tsv                             when and where each file was downloaded
├── checksums.md5                                md5 checksums of all files
├── case1/
│   ├── causalccc_server_eaP9Br0llvzxN07zZZEC/   files of the authors' web server job
│   ├── nichenet_seurat/seuratObj.rds            processed Medaglia data (Seurat object)
│   └── geo_GSE104054/                           raw data (optional, --with-geo)
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

Fibroblasts, dendritic cells and T cells in atopic dermatitis skin (He *et al.*, *J Allergy Clin Immunol* 2020).
L-R pairs inferred with CellChat. **Not downloaded yet.**

- The authors' CausalCCC input files are on the web server (job `w0788tM1wzXGHLzNaKKa`).
- The CellChat dataset link on figshare refused scripted downloads on 2026-10-05.

## Other datasets (optional)

Not downloaded by the scripts.

| Dataset                        | Used for        | Source                                                                    | How to download                                                                            | Size   |
| ------------------------------ | --------------- | ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | ------ |
| CausalCCC demo (Seurat object) | Getting started | [CausalCCC tutorial](https://miic.curie.fr/tutorial_causalCCC.php)        | `wget -c -P data/demo https://miic.curie.fr/datasets/causalCCC_demo_Seurat.RData`          | 274 MB |
