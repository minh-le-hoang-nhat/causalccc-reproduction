#!/usr/bin/env bash
# =============================================================================
# Download the data for case study 1 of CausalCCC (Dupuis et al., NAR 2025):
# CD4+ (sender) -> CD8+ (receiver) T cells, inguinal lymph nodes,
# data from Medaglia et al., Science 2017, L-R pairs from NicheNet.
#
# Usage (from the root of the repository):
#   bash scripts/download_case1.sh              # everything except GEO raw data
#   bash scripts/download_case1.sh --with-geo   # also GEO raw data (+ 89 MB)
#
# Safe to run again: complete files are skipped, interrupted downloads resume.
# Sources are described in data/README.md.
# =============================================================================

set -euo pipefail   # stop at the first error, unset variable or failed pipe


# --- 1. Settings -------------------------------------------------------------

DATA_DIR="data"
LOG="${DATA_DIR}/download_log.tsv"

WITH_GEO=false
if [[ "${1:-}" == "--with-geo" ]]; then WITH_GEO=true; fi

# CausalCCC web server job used by the authors for case study 1 (Figure 2)
JOB_ID="eaP9Br0llvzxN07zZZEC"
JOB_URL="https://miic.curie.fr/server/server_outputs/${JOB_ID}"
JOB_DIR="${DATA_DIR}/case1/causalccc_server_${JOB_ID}"

# Zenodo records (official checksums are published by Zenodo)
SEURAT_URL="https://zenodo.org/records/3531889/files"   # NicheNet demo Seurat object
NICHENET_URL="https://zenodo.org/records/7074291/files" # NicheNet v2 networks

GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104054/suppl"


# --- 2. Checks before starting -----------------------------------------------

if [[ ! -f README.md || ! -d scripts ]]; then
  echo "ERROR: run this script from the root of the repository." >&2
  exit 1
fi

for tool in wget md5sum unzip; do
  if ! command -v "$tool" > /dev/null; then
    echo "ERROR: '$tool' is not installed (try: sudo apt install $tool)." >&2
    exit 1
  fi
done


# --- 3. Helper function: download one file -----------------------------------
# Arguments: URL, destination path, expected md5 checksum (optional).
# The file is downloaded as <name>.part and renamed only when complete,
# so a file without .part is always a finished download.

download() {
  local url="$1"
  local dest="$2"
  local expected_md5="${3:-}"

  if [[ -f "$dest" ]]; then
    echo "  already here: $dest"
    return 0
  fi

  echo "  downloading:  $dest"
  mkdir -p "$(dirname "$dest")"
  wget --continue --quiet --show-progress --tries=5 --output-document="${dest}.part" "$url"
  mv "${dest}.part" "$dest"

  if [[ -n "$expected_md5" ]]; then
    local actual_md5
    actual_md5=$(md5sum "$dest" | cut -d " " -f 1)
    if [[ "$actual_md5" != "$expected_md5" ]]; then
      mv "$dest" "${dest}.corrupt"
      echo "ERROR: wrong checksum for $dest (renamed to .corrupt). Delete it and run again." >&2
      exit 1
    fi
    echo "  checksum OK"
  fi

  # Record what was downloaded, from where, and when
  if [[ ! -f "$LOG" ]]; then
    printf "date\tfile\turl\n" > "$LOG"
  fi
  printf "%s\t%s\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$dest" "$url" >> "$LOG"
}


# --- 4. Files used by the CausalCCC authors (web server job) -----------------
# Small files (about 3 MB in total). The server is not an archive: keep them.

echo "[1/4] CausalCCC authors' files (job ${JOB_ID})"

AUTHOR_FILES=(
  input.zip                                    # mosaic table (input.tsv)
  clean_interact_edges.tsv                     # the 18 NicheNet L-R pairs
  clean_state_order.tsv                        # variable types and groups
  script.sh                                    # MIIC settings used on the server
  whole_run.log                                # log of the server run
  edgesList.miic.summary.txt                   # result network (April 2025 version)
  edgesList.miic.summary.txt_save              # result network (January 2025 run)
  NicheNet_edgesList.miic.summary_relativeLR.txt
  NicheNet_edgesList.miic.summary_relativeLR_hi.txt
  adjacencyMatrix.miic.txt
  adjacencyMatrix.miic.orientProba.txt
  orientations.probaArrowhead.txt
  cutpoints.json
  layout.json
  miic_graph.json
  groupsNew.txt
  template_category_order.tsv
)

for f in "${AUTHOR_FILES[@]}"; do
  download "${JOB_URL}/${f}" "${JOB_DIR}/${f}"
done

if [[ ! -f "${JOB_DIR}/input.tsv" ]]; then
  unzip -q "${JOB_DIR}/input.zip" -d "$JOB_DIR"
  echo "  unzipped:     ${JOB_DIR}/input.tsv"
fi


# --- 5. Processed expression data: NicheNet's Seurat object (536 MB) ---------

echo "[2/4] Seurat object from Medaglia et al. (Zenodo 3531889)"

download "${SEURAT_URL}/seuratObj.rds?download=1" \
         "${DATA_DIR}/case1/nichenet_seurat/seuratObj.rds" \
         "c5f94621f4644d17ae6c5afdc0ab3bf1"


# --- 6. NicheNet v2 prior knowledge, mouse (about 270 MB) ---------------------

echo "[3/4] NicheNet v2 mouse networks (Zenodo 7074291)"

download "${NICHENET_URL}/lr_network_mouse_21122021.rds?download=1" \
         "${DATA_DIR}/nichenet_v2/lr_network_mouse_21122021.rds" \
         "cf33ee8b6bf84bdf2d11cab9c8f94b9e"

download "${NICHENET_URL}/ligand_target_matrix_nsga2r_final_mouse.rds?download=1" \
         "${DATA_DIR}/nichenet_v2/ligand_target_matrix_nsga2r_final_mouse.rds" \
         "ac80d846fe0bfc4879a5b52ca85ffeb9"

download "${NICHENET_URL}/weighted_networks_nsga2r_final_mouse.rds?download=1" \
         "${DATA_DIR}/nichenet_v2/weighted_networks_nsga2r_final_mouse.rds" \
         "6620b465a566a9e556ee0f49f5f19cfa"


# --- 7. Optional: raw data from GEO (89 MB) -----------------------------------

if [[ "$WITH_GEO" == true ]]; then
  echo "[4/4] GEO GSE104054 raw data"
  download "${GEO_URL}/GSE104054_RAW.tar"         "${DATA_DIR}/case1/geo_GSE104054/GSE104054_RAW.tar"
  download "${GEO_URL}/GSE104054_metadata.txt.gz" "${DATA_DIR}/case1/geo_GSE104054/GSE104054_metadata.txt.gz"
else
  echo "[4/4] GEO raw data skipped (add --with-geo to download it)"
fi


# --- 8. Check, then record, the checksums of everything we have ---------------
# Data files must never be modified. If a previous run recorded checksums,
# first make sure the files are unchanged, then record the new list.
# You can also check by hand at any time:  cd data && md5sum -c checksums.md5

if [[ -f "${DATA_DIR}/checksums.md5" ]]; then
  echo "Checking that previously downloaded files are unchanged..."
  if ! (cd "$DATA_DIR" && md5sum --check --quiet --ignore-missing checksums.md5); then
    echo "ERROR: the files listed above changed since the last run." >&2
    echo "       Delete them and run this script again to get clean copies." >&2
    exit 1
  fi
fi

echo "Computing checksums (a few seconds)..."
(
  cd "$DATA_DIR"
  find case1 nichenet_v2 -type f ! -name "*.part" ! -name "*.corrupt" | sort | xargs md5sum
) > "${DATA_DIR}/checksums.md5"

echo "Done. Files are in ${DATA_DIR}/, log in ${LOG}, checksums in ${DATA_DIR}/checksums.md5"
