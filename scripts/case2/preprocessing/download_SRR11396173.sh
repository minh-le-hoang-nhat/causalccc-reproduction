#!/usr/bin/env bash
# Téléchargement FASTQ (ENA) de SRR11396173 = GSM4430473, S15_NL (He et al. 2020)
# Seul échantillon "petit" présent dans le subset CellChat du case study 2 (321 cellules S15_*).
# Usage (dans WSL) : bash scripts/case2/preprocessing/download_SRR11396173.sh   puis   bash scripts/case2/preprocessing/run_kb.sh SRR11396173
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"   # repository root
DATA="$PROJ/data/case2/preprocessing"
SRR=SRR11396173
BASE="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR113/073/$SRR"
declare -A SIZE=([1]=1983560996 [2]=6647337797)   # tailles annoncées par ENA (fastq_bytes)
mkdir -p "$DATA/fastq"
for R in 1 2; do
  F="$DATA/fastq/${SRR}_${R}.fastq.gz"
  wget -c -O "$F" "$BASE/${SRR}_${R}.fastq.gz"
  GOT=$(stat -c %s "$F")
  if [ "$GOT" != "${SIZE[$R]}" ]; then
    echo "ERREUR : $F fait $GOT octets, attendu ${SIZE[$R]}" >&2; exit 1
  fi
  echo "OK : $F ($GOT octets)"
done
