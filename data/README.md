# Data

Data files are **not tracked by Git** (size, licensing); only this README is. It records where every dataset comes from
and how to obtain it, so that `data/` can be rebuilt from scratch.

**Never modify a file in this folder.** Notebooks read from `data/` and write to `results/`.

## Rebuilding this folder

| Dataset      | Command                                           | Size                                     |
| ------------ | ------------------------------------------------- | ---------------------------------------- |
| Case study 1 | `bash scripts/case1/download_case1.sh`            | about 770 MB                             |
|              | `bash scripts/case1/download_case1.sh --with-geo` | also GEO raw data, about 860 MB in total |
| Case study 2 | `bash scripts/case2/download_case2.sh`            | about 5 MB                               |
| Thymus       | `bash scripts/thymus/download_thymus.sh`          | about 4.2 GB, plus the unzipped archives |
|              | `bash scripts/thymus/download_thymus.sh --skip-zenodo` | table S5 check, CellPhoneDB and CollecTRI files only, about 12 MB |

Run the scripts from the root of the repository. They can be run again at any time: complete files are skipped,
interrupted downloads resume, and files that were modified since the last run are reported. Each download is recorded
in `download_log.tsv` (date, file, URL) and all checksums in `checksums.md5` (case study 1), `checksums_case2.md5`
(case study 2) or `checksums_thymus.md5` (thymus). Check by hand with `cd data && md5sum -c checksums.md5`.

## Layout

```
data/
├── README.md                                    (this file, tracked by Git)
├── download_log.tsv                             when and where each file was downloaded
├── checksums.md5                                md5 checksums of the case study 1 files
├── checksums_case2.md5                          md5 checksums of the case study 2 files
├── checksums_thymus.md5                         md5 checksums of the thymus files
├── case1/
│   ├── causalccc_server_eaP9Br0llvzxN07zZZEC/   files of the authors' web server job
│   ├── causalccc_server_ours_*_params/          our own server jobs (not needed, see below)
│   ├── nichenet_seurat/seuratObj.rds            processed Medaglia data (Seurat object)
│   └── geo_GSE104054/                           raw data (optional, --with-geo)
├── case2/
│   ├── causalccc_server_w0788tM1wzXGHLzNaKKa/   files of the authors' web server job
│   ├── causalccc_server_<code>/                 other server jobs, downloaded by the tutorial (section 8b)
│   ├── cellchat/data_humanSkin_CellChat.rda     CellChat example dataset (copied by hand, see below)
│   └── preprocessing/                           raw reads (fastq/), kallisto|bustools output (kb_out/), write/, cache/ – see notebooks/case2/preprocessing/README.md
├── nichenet_v2/                                 NicheNet v2 mouse prior-knowledge networks
└── thymus/
    ├── zenodo_5500511/                          atlas (raw counts, annotated matrices), sample metadata, authors' code
    ├── park2020_supplement/                     supplementary tables of Park et al. (copied by hand, see below)
    ├── cellphonedb_v2_db/db_2019-09/, db_2019-07/   CellPhoneDB v2 database files (two versions, see below)
    └── collectri/collectri_human.tsv            CollecTRI transcription factor regulons, human (see below)
```

## Sources – case study 1

CD4⁺ (sender) → CD8⁺ (receiver) T cells, T-cell area of inguinal lymph nodes, steady state and 72 h after LCMV infection
(Medaglia *et al.*, *Science* 2017, [doi:10.1126/science.aao4277](https://doi.org/10.1126/science.aao4277)).
L-R pairs inferred with NicheNet. Figure 2 of Dupuis *et al.* 2025.

| Folder                                         | Source                                                                                                                                        | Content                                                                                                                                                                                                           | Size   |
| ---------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ |
| `case1/causalccc_server_eaP9Br0llvzxN07zZZEC/` | CausalCCC web server job, linked as "NicheNet CausalCCC input files" in the [tutorial](https://miic.curie.fr/tutorial_causalCCC.php#datasets) | Authors' input mosaic table, L-R pairs, state order, MIIC settings, run log and result network. Job run on 2025-01-21.                                                                                            | 4.8 MB |
| `case1/nichenet_seurat/`                       | Zenodo [10.5281/zenodo.3531889](https://doi.org/10.5281/zenodo.3531889) (NicheNet demo data)                                                  | Seurat object of the Medaglia data after Seurat alignment, steady state + LCMV                                                                                                                                    | 512 MB |
| `nichenet_v2/`                                 | Zenodo [10.5281/zenodo.7074291](https://doi.org/10.5281/zenodo.7074291) (NicheNet v2)                                                         | Mouse L-R network, ligand–target matrix, weighted networks                                                                                                                                                        | 257 MB |
| `case1/geo_GSE104054/` (optional)              | GEO [GSE104054](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104054)                                                                 | `GSE104054_RAW.tar`: raw MARS-seq UMI tables (149 amplification batches). `GSE104054_metadata.txt.gz`: per-well metadata (57,216 wells; lymph node, spleen and tumour; naive and LCMV). Superset of case study 1. | 89 MB  |

Checksums of the Zenodo files are verified against the values published by Zenodo. GEO and the CausalCCC server publish
no checksums: ours are recorded at first download.

### Main files of the authors' job

| File                              | What it is                                                                                                                                |
| --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `input.tsv` (from `input.zip`)    | Mosaic table: 4,207 cells (2,562 CD4⁺ senders, 1,645 CD8⁺ receivers) × 144 genes, raw UMI counts, `NA` outside each cell type's own genes |
| `clean_interact_edges.tsv`        | The 18 NicheNet L-R pairs (10 ligands, 8 receptors), no header, no score                                                                  |
| `clean_state_order.tsv`           | Type (continuous/discrete) and group (ligand, receptor, sender gene, receiver gene) of each variable                                      |
| `script.sh`, `whole_run.log`      | Exact MIIC settings used by the server, and the log of the run                                                                            |
| `edgesList.miic.summary.txt`      | Result network (364 retained edges, `type = P`). Edited on 2025-04-24.                                                                    |
| `edgesList.miic.summary.txt_save` | Result network as produced by the 2025-01-21 run. Same edges; only the values of the 18 L-R links differ from the edited version.         |

### Notes and open questions

- The tutorial's "NicheNet dataset" link points to the NicheNet v1 model files (Zenodo 3260758, human), not to the
  expression data. `seuratObj.rds` is the origin of the authors' input: their counts are identical to ours for all 144
  variables, and four genes carry a newer name in their table (Mir703, Rack1, Ifi209, Tpt1-ps3; notebook 03, section 1).
- The mosaic has no condition column, although the paper mentions steady state and infection: both are probably pooled.
- The NicheNet version and parameters used to obtain the 18 L-R pairs are not stated, though, it is likely we won't need NicheNet v2 networks for now.
- The CausalCCC server is not an archive: these files could disappear, hence the local copy and its checksums.
- `case1/causalccc_server_ours_default_params/` and `…_updated_params/`: networks returned by the web server on
  2026-10-08 for our inputs of notebooks 03 and 04, uploaded by hand (no download script). The server ran with its
  default settings, not those of the authors' job (no confidence filter, `confidence` is `NA`), so the notebooks run MIIC
  locally with the authors' settings instead. Not needed to run the notebooks; kept for reference.

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

## Sources – thymus

Human thymus cell atlas, 7 post-conception weeks (PCW) to adult (Park *et al.*, *Science* 2020,
[doi:10.1126/science.aay3224](https://doi.org/10.1126/science.aay3224), open text
[PMC7611066](https://pmc.ncbi.nlm.nih.gov/articles/PMC7611066/)). Used for the early thymic niche
(`notebooks/thymus_early/`: stromal cells ↔ early lymphoid cells at 7–11 PCW, CausalCCC in both directions), with the
ligand–receptor pairs of the authors' CellPhoneDB analysis (table S5).

| Folder                              | Source                                                                                     | Content                                                                                                                       | Size    |
| ----------------------------------- | ------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------- | ------- |
| `thymus/zenodo_5500511/`            | Zenodo [10.5281/zenodo.5500511](https://doi.org/10.5281/zenodo.5500511) (version 1.0.2, 2021-09-11) | `HTA07.A01.v02.entire_data_raw_count.h5ad` (raw counts, 272,554 cells), `thymus_annotated_matrix_files.zip` (annotated matrices; `HTA08.v01.A05.Science_human_fig1.h5ad` is the Fig. 1 atlas, 255,901 cells, log-normalised), `sample_metadata_fix.xlsx` (table S1: per-sample QC, sort, sex, chemistry, enzyme; no age), `thymus_code_package.zip` (authors' code). The archives are unzipped into folders of the same name. `thymus_vdj.zip` (TCR data, 0.5 GB) is not downloaded. | 4.2 GB  |
| `thymus/park2020_supplement/`       | Supplementary material of Park *et al.* (Science website)                                 | `aay3224_table-s5.csv`: CellPhoneDB means, 710 interactions × 53 × 53 cell-type pairs (no p-values). Copied by hand (the Science supplement refuses scripted downloads); the script checks its md5 (`f2b5c69fa8468414a9db1df93e3aad5f`). Tables S1–S4 and data S1 are copied by hand too and not used by the notebooks. | 17 MB   |
| `thymus/cellphonedb_v2_db/db_2019-09/` | GitHub [ventolab/cellphonedb-data](https://github.com/ventolab/cellphonedb-data), tag `v2.0.0`, commit `05d95ad` (2019-09-05) | `protein_input.csv` (receptor flags), `complex_input.csv` (subunits), `gene_input.csv` (UniProt → gene), `interaction_input.csv`. Same files as in CellPhoneDB code releases v2.1.0 to v2.1.7. | < 1 MB  |
| `thymus/cellphonedb_v2_db/db_2019-07/` | GitHub [ventolab/CellphoneDB](https://github.com/ventolab/CellphoneDB), tag `v2.0.6`, commit `095d4ae` (2019-07-23), `cellphonedb/src/core/data/` | Same four files, as bundled with code releases v2.0.0 to v2.0.6. | < 1 MB  |
| `thymus/collectri/`                 | CollecTRI (Müller-Dott *et al.*, *Nucleic Acids Res* 2023, [doi:10.1093/nar/gkad841](https://doi.org/10.1093/nar/gkad841)), from the [OmniPath](https://omnipathdb.org) web service, downloaded 2026-10-17 | `collectri_human.tsv`: 64,516 rows, 62,411 distinct transcription factor → target gene pairs, including complexes such as `NFKB1_RELB` (columns `source_genesymbol`, `target_genesymbol`, `is_stimulation`, `is_inhibition`, literature references). Gene sets of the activity scores (notebook 10). md5 `2f626e289712311596bc5a6687e44384`. | 12 MB   |

Checksums of the Zenodo files are verified against the values published by Zenodo. GitHub and OmniPath publish no
checksums: ours are recorded at first download.

### Notes and open questions

- The expression values come from the raw-count file, the annotations and metadata from the `.obs` of the Fig. 1 atlas:
  every atlas cell is in the raw-count file, with the same sample, donor, sort and chemistry (`notebooks/thymus_early/`,
  notebook 01). The `.csv` next to each annotated `.h5ad` is a copy of its `.obs`.
- The chemistry column (`Method`) of `sample_metadata_fix.xlsx` agrees with the atlas for only about half of the
  samples; the atlas column `method` agrees with the sample names, so it is the one used (notebook 01).
- Park *et al.* cite CellPhoneDB v2 but give no database version. The two v2 databases above give the same receptor
  and integrin flags for all 519 partners of table S5, and `db_2019-09` is used. Table S5 was nevertheless computed with
  neither: only 259 of its 710 interaction IDs are in them, and its 53 cell types differ from the published annotation
  (e.g. `TEC_early`, `Fb_early`; notebook 01).
- Only the CellPhoneDB database files are used by `notebooks/thymus_early/`; the `cellphonedb` package is not installed
  in `envs/`. The CellPhoneDB re-run of `scripts/thymus/01_*.py` to `05_*.py` uses CellPhoneDB 5.0.1 with its own
  database, and reads the atlas and table S5 from other paths (`data/thymus/HTA08_v01_A05_Science_human_fig1.h5ad`,
  `data/thymus/aay3224_table-s5.csv`); the two layouts are still to be harmonised.
- CollecTRI has no versioned release: the OmniPath service returns the current version, so a later download can differ.
  The script keeps the first download; compare a new one with the md5 above. The downloaded version is larger than the
  one described in the paper (43,175 interactions) and cites `CollecTRI2` references, so it is a later update.
  `OmnipathR::collectri()` (OmnipathR 4.1.0, installed with LIANA) is not used: it fails with
  `Join columns in y must be present in the data … ncbi_tax_id`.
- In CollecTRI, the NF-κB family (RELA, RELB, REL, NFKB1, NFKB2) share 803 targets through their complexes, so the
  RELB and NFKB2 regulons do not separate the non-canonical from the canonical NF-κB pathway.

## Server jobs downloaded by the tutorial

Section 8b of the tutorial (`notebooks/tutorial/causalccc_tutorial.qmd`, rendered with `-P server_job:<link or code>`)
compares a CausalCCC web server job with a local run. It downloads six files of the job (`script.sh`,
`edgesList.miic.summary.txt`, `input.zip` unzipped into `input.tsv`, `clean_state_order.tsv`,
`clean_interact_edges.tsv`, `whole_run.log`) from `https://miic.curie.fr/server/server_outputs/<code>/` into
`case1/causalccc_server_<code>/` for the authors' case 1 job and `case2/causalccc_server_<code>/` otherwise, skips the
files already there (so the authors' folders are never overwritten) and logs each download in `download_log.tsv`.
These jobs are not downloaded by the scripts; render the tutorial with the same link to get them back while the server
keeps them.

## Other datasets (optional)

Not downloaded by the scripts.

| Dataset                        | Used for        | Source                                                             | How to download                                                                   | Size   |
| ------------------------------ | --------------- | ------------------------------------------------------------------ | --------------------------------------------------------------------------------- | ------ |
| CausalCCC demo (Seurat object) | Getting started | [CausalCCC tutorial](https://miic.curie.fr/tutorial_causalCCC.php) | `wget -c -P data/demo https://miic.curie.fr/datasets/causalCCC_demo_Seurat.RData` | 274 MB |
