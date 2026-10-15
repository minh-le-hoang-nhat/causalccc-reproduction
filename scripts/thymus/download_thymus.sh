#!/usr/bin/env bash
# =============================================================================
# Download the data for the thymus analyses (early and medullary niches):
# human thymus cell atlas of Park et al., Science 2020 (doi:10.1126/science.aay3224),
# and the CellPhoneDB v2 database files needed to read their table S5.
#
# Usage (from the root of the repository):
#   bash scripts/thymus/download_thymus.sh                 # everything (about 4.2 GB)
#   bash scripts/thymus/download_thymus.sh --skip-zenodo   # table S5 check + CellPhoneDB only (< 1 MB)
#
# Safe to run again: complete files are skipped, interrupted downloads resume.
# Sources are described in data/README.md.
# =============================================================================

set -euo pipefail   # stop at the first error, unset variable or failed pipe


# --- 1. Settings -------------------------------------------------------------

DATA_DIR="data"
LOG="${DATA_DIR}/download_log.tsv"
CHECKSUMS="checksums_thymus.md5"   # separate from the case studies' checksum files
THYMUS_DIR="${DATA_DIR}/thymus"

SKIP_ZENODO=false
if [[ "${1:-}" == "--skip-zenodo" ]]; then SKIP_ZENODO=true; fi

# Table S5 of Park et al. (CellPhoneDB means). The Science supplement cannot be
# downloaded by a script: it is copied by hand, then checked here.
S5_FILE="${THYMUS_DIR}/park2020_supplement/aay3224_table-s5.csv"
S5_MD5="f2b5c69fa8468414a9db1df93e3aad5f"

# CellPhoneDB v2 database files, pinned to commits (a tag can move, a commit cannot).
# Two different v2 databases exist; they give the same receptor flags for the S5 partners, and
# notebooks/thymus_early/ uses db_2019-09 (notebook 03).
#   db_2019-09: data release v2.0.0 (2019-09-05), the database of code releases v2.1.0 to v2.1.7
#   db_2019-07: database bundled with code release v2.0.6 (2019-07-23), same as v2.0.0 to v2.0.6
CPDB_DIR="${THYMUS_DIR}/cellphonedb_v2_db"
CPDB_2019_09_URL="https://raw.githubusercontent.com/ventolab/cellphonedb-data/05d95adc6ebb3546f302c8ec422d1b86594cd8df/data"
CPDB_2019_07_URL="https://raw.githubusercontent.com/ventolab/CellphoneDB/095d4ae851e14355677b4f44f8fd37f40ba78fad/cellphonedb/src/core/data"
CPDB_FILES=(
  protein_input.csv       # proteins, with the receptor flag
  complex_input.csv       # complexes and their subunits (UniProt IDs)
  gene_input.csv          # UniProt ID -> gene name
  interaction_input.csv   # interactions, with partner order (partner_a, partner_b)
)

# Zenodo record 5500511 (version 1.0.2 of the authors' data and code; official checksums published by Zenodo)
ZENODO_URL="https://zenodo.org/records/5500511/files"
ZENODO_DIR="${THYMUS_DIR}/zenodo_5500511"


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
# Same as in download_case1.sh. Arguments: URL, destination path, expected md5 checksum (optional).
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

# Unzip an archive into its own folder (named after the archive), once.
unzip_once() {
  local archive="$1"
  local dest="${archive%.zip}"
  if [[ -d "$dest" ]]; then
    echo "  already unzipped: $dest/"
    return 0
  fi
  rm -rf "${dest}.part"   # left over by an interrupted run of this script
  unzip -q "$archive" -d "${dest}.part"
  mv "${dest}.part" "$dest"
  echo "  unzipped:     $dest/"
}


# --- 4. Table S5 (copied by hand) --------------------------------------------

echo "[1/3] Table S5 of Park et al."

if [[ ! -f "$S5_FILE" ]]; then
  echo "ERROR: $S5_FILE is missing. Download table S5 from the supplementary materials of" >&2
  echo "       doi:10.1126/science.aay3224 and put it there (see data/README.md)." >&2
  exit 1
fi
if [[ "$(md5sum "$S5_FILE" | cut -d " " -f 1)" != "$S5_MD5" ]]; then
  echo "ERROR: wrong checksum for $S5_FILE (expected $S5_MD5)." >&2
  exit 1
fi
echo "  checksum OK:  $S5_FILE"


# --- 5. CellPhoneDB v2 database files (less than 1 MB) ------------------------
# GitHub publishes no md5: ours are recorded at first download (section 7).

echo "[2/3] CellPhoneDB v2 database files"

for f in "${CPDB_FILES[@]}"; do
  download "${CPDB_2019_09_URL}/${f}" "${CPDB_DIR}/db_2019-09/${f}"
  download "${CPDB_2019_07_URL}/${f}" "${CPDB_DIR}/db_2019-07/${f}"
done


# --- 6. Park et al. atlas and code (Zenodo 5500511, about 4.2 GB) -------------

if [[ "$SKIP_ZENODO" == true ]]; then
  echo "[3/3] Zenodo files skipped (run without --skip-zenodo to download them)"
else
  echo "[3/3] Thymus atlas, metadata and code (Zenodo 5500511)"
  download "${ZENODO_URL}/sample_metadata_fix.xlsx?download=1" \
           "${ZENODO_DIR}/sample_metadata_fix.xlsx" \
           "16fcb558c27016665f0c946d29c2b9aa"
  download "${ZENODO_URL}/thymus_code_package.zip?download=1" \
           "${ZENODO_DIR}/thymus_code_package.zip" \
           "57d189b2c87501a63d2771776df4f2c8"
  download "${ZENODO_URL}/thymus_annotated_matrix_files.zip?download=1" \
           "${ZENODO_DIR}/thymus_annotated_matrix_files.zip" \
           "5eb2cb0ccc7f87dcddb68d785517de9d"
  download "${ZENODO_URL}/HTA07.A01.v02.entire_data_raw_count.h5ad?download=1" \
           "${ZENODO_DIR}/HTA07.A01.v02.entire_data_raw_count.h5ad" \
           "496cfaeb6d227b2bb32d52abd2249085"

  unzip_once "${ZENODO_DIR}/thymus_code_package.zip"
  unzip_once "${ZENODO_DIR}/thymus_annotated_matrix_files.zip"
fi


# --- 7. Check, then record, the checksums of everything we have ---------------
# Data files must never be modified. If a previous run recorded checksums,
# first make sure the files are unchanged, then record the new list.
# You can also check by hand at any time:  cd data && md5sum -c checksums_thymus.md5

if [[ -f "${DATA_DIR}/${CHECKSUMS}" ]]; then
  echo "Checking that previously downloaded files are unchanged..."
  if ! (cd "$DATA_DIR" && md5sum --check --quiet --ignore-missing "$CHECKSUMS"); then
    echo "ERROR: the files listed above changed since the last run." >&2
    echo "       Delete them and run this script again to get clean copies." >&2
    exit 1
  fi
fi

echo "Computing checksums (up to a few minutes after the Zenodo download)..."
(
  cd "$DATA_DIR"
  find thymus -type f ! -name "*.part" ! -name "*.corrupt" ! -path "*.part/*" -print0 | sort -z | xargs -0 md5sum
) > "${DATA_DIR}/${CHECKSUMS}"

echo "Done. Files are in ${THYMUS_DIR}/, log in ${LOG}, checksums in ${DATA_DIR}/${CHECKSUMS}"
