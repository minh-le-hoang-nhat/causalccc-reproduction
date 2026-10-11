# Installation

Tested on Ubuntu 24.04 (native or WSL2). macOS should work with the adjustments noted below.

## Requirements

- **git**, and a GitHub account with access to this repository.
- **conda** with the conda-forge channel only ([Miniforge](https://github.com/conda-forge/miniforge) recommended).
- **R 4.6** ([rig](https://github.com/r-lib/rig) recommended: it configures pre-compiled package binaries and installs `pak`).
- **C/C++ compilers** to build `miic` from source (Ubuntu: `sudo apt install build-essential`; macOS: `xcode-select --install`).

## 0. Get the code

```bash
git clone https://github.com/minh-le-hoang-nhat/causalccc-reproduction.git
cd causalccc-reproduction
```

All commands below are run from this folder.

## 1. Python environment

```bash
conda env create -f envs/environment.yml
conda activate csbc-causalccc
```

Key package versions are pinned in `envs/environment.yml`. Never `conda install` directly into this environment:
add the package (with its version) to `envs/environment.yml`, then run

```bash
conda env update -n csbc-causalccc -f envs/environment.yml
conda list -n csbc-causalccc > envs/python_versions.txt
```

and commit both files.

## 2. R packages

`.Rprofile` points R to the Posit Package Manager snapshot of CRAN from **2026-10-03**, so every installation resolves
to the same package versions. R reads it automatically when it starts in this folder. On Linux it selects pre-compiled
binaries for your distribution; on macOS and Windows it uses the standard snapshot.

```bash
Rscript envs/install_R_packages.R
```

This installs Seurat, tidyverse and rjson from the CRAN snapshot, and `miic` (branch `causalccc`), `causalCCC` and
`liana` (ligand-receptor inference used by `causalCCC.links()`) from GitHub at fixed commits. It stops with an error if
the installed commits do not match. LIANA also installs several Bioconductor packages (scran, scater, OmnipathR, …);
their versions follow the Bioconductor release that matches R 4.6.

## 3. Check

```bash
python -c "import scanpy; print(scanpy.__version__)"     # 1.12.4 (with csbc-causalccc active)
Rscript -e 'library(Seurat); library(miic); library(causalCCC); packageVersion("Seurat")'
```

Compare your versions with `envs/python_versions.txt` and `envs/R_versions.txt`.

## Known issues

- The `causalccc` branch of `miic` uses `rjson` and `tidyverse` without declaring them; the install script installs them first.
- The Python helper functions shown in the CausalCCC online tutorial are not compatible with pandas ≥ 2
  (`DataFrame.append` was removed). This project uses the R package `causalCCC` for the CausalCCC step.
- **Positron users:** set `python.useBundledIpykernel` to `false` and restart Positron. Otherwise the bundled ipykernel
  shadows the environment's `typing_extensions` and `import scanpy` fails with
  `ImportError: cannot import name 'sentinel'` ([posit-dev/positron#15574](https://github.com/posit-dev/positron/issues/15574)).
