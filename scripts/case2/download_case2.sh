#!/usr/bin/env bash
# =============================================================================
# Download the data for case study 2 of CausalCCC (Dupuis et al., NAR 2025):
# fibroblasts, dendritic cells and T cells in atopic dermatitis skin,
# data from He et al., J Allergy Clin Immunol 2020, L-R pairs from CellChat.
#
# Usage (from the root of the repository):
#   bash scripts/case2/download_case2.sh
#
# Safe to run again: complete files are skipped, interrupted downloads resume.
# Sources are described in data/README.md.
# =============================================================================

set -euo pipefail   # stop at the first error, unset variable or failed pipe


# --- 1. Settings -------------------------------------------------------------

DATA_DIR="data"
LOG="${DATA_DIR}/download_log.tsv"
CHECKSUMS="checksums_case2.md5"   # separate from case 1's checksums.md5

# CausalCCC web server job used by the authors for case study 2 (Figure 3)
JOB_ID="w0788tM1wzXGHLzNaKKa"
JOB_URL="https://miic.curie.fr/server/server_outputs/${JOB_ID}"
JOB_DIR="${DATA_DIR}/case2/causalccc_server_${JOB_ID}"


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
# Same as in download_case1.sh. Arguments: URL, destination path.
# The file is downloaded as <name>.part and renamed only when complete,
# so a file without .part is always a finished download.

download() {
  local url="$1"
  local dest="$2"

  if [[ -f "$dest" ]]; then
    echo "  already here: $dest"
    return 0
  fi

  echo "  downloading:  $dest"
  mkdir -p "$(dirname "$dest")"
  wget --continue --quiet --show-progress --tries=5 --output-document="${dest}.part" "$url"
  mv "${dest}.part" "$dest"

  # Record what was downloaded, from where, and when
  if [[ ! -f "$LOG" ]]; then
    printf "date\tfile\turl\n" > "$LOG"
  fi
  printf "%s\t%s\t%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$dest" "$url" >> "$LOG"
}


# --- 4. Files used by the CausalCCC authors (web server job) -----------------
# Small files (about 5 MB in total). The server is not an archive: keep them.
# Unlike case 1, this job has no edgesList.miic.summary.txt_save.

echo "[1/2] CausalCCC authors' files (job ${JOB_ID})"

AUTHOR_FILES=(
  input.zip                          # mosaic table (input.tsv)
  clean_interact_edges.tsv           # the 43 CellChat L-R pairs
  clean_state_order.tsv              # variable types and groups
  script.sh                          # MIIC settings used on the server
  whole_run.log                      # log of the server run
  log.txt                            # log of the server wrapper
  edgesList.miic.summary.txt         # result network
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


# --- 5. Check, then record, the checksums of everything we have ---------------
# Data files must never be modified. If a previous run recorded checksums,
# first make sure the files are unchanged, then record the new list.
# You can also check by hand at any time:  cd data && md5sum -c checksums_case2.md5

echo "[2/2] Checksums"

if [[ -f "${DATA_DIR}/${CHECKSUMS}" ]]; then
  echo "Checking that previously downloaded files are unchanged..."
  if ! (cd "$DATA_DIR" && md5sum --check --quiet --ignore-missing "$CHECKSUMS"); then
    echo "ERROR: the files listed above changed since the last run." >&2
    echo "       Delete them and run this script again to get clean copies." >&2
    exit 1
  fi
fi

(
  cd "$DATA_DIR"
  find case2 -type f ! -name "*.part" ! -name "*.corrupt" | sort | xargs md5sum
) > "${DATA_DIR}/${CHECKSUMS}"

echo "Done. Files are in ${JOB_DIR}/, log in ${LOG}, checksums in ${DATA_DIR}/${CHECKSUMS}"
