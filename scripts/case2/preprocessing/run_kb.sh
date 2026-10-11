#!/usr/bin/env bash
# Quantification FASTQ -> matrice de comptage (kallisto|bustools), 10x Chromium v2
# Usage (dans WSL) : bash scripts/case2/preprocessing/run_kb.sh SRR11396171 [SRR... ]
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"   # repository root
DATA="$PROJ/data/case2/preprocessing"
KB=~/kb-env/bin/kb
REF=~/kb-ref/human            # index sur le FS Linux (plus rapide)
mkdir -p "$REF"
cd ~   # dossiers tmp de kb sur le FS Linux (sinon le rename vers ~/kb-ref échoue)
if [ ! -f "$REF/index.idx" ]; then
  "$KB" ref --tmp ~/kb-tmp-ref -d human -i "$REF/index.idx" -g "$REF/t2g.txt"
fi
for SRR in "$@"; do
  OUT="$DATA/kb_out/$SRR"
  mkdir -p "$OUT"
  "$KB" count --tmp ~/kb-tmp-"$SRR" -i "$REF/index.idx" -g "$REF/t2g.txt" -x 10xv2 \
    -o "$OUT" -t 8 -m 5G --filter bustools \
    "$DATA/fastq/${SRR}_1.fastq.gz" "$DATA/fastq/${SRR}_2.fastq.gz"
done
