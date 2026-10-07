#!/bin/bash
# MIIC on our rebuilt case study 2 input (notebook 06), on the IFB Core cluster (Slurm).
# Prerequisites: R environment installed as in docs/INSTALL.md, notebook 05 and 06 (up to the input files) run once,
# so that results/case2/06_rebuilt_network/{input,state_order,interact_edges}.tsv exist.
# Submit from the repository root:  sbatch scripts/case2/run_miic_case2_rebuilt.sbatch
# Then copy back results/case2/06_rebuilt_network/miic_result.rds: notebook 06 reloads it.
#SBATCH --job-name=miic_case2_rebuilt
#SBATCH --account=YOUR_IFB_PROJECT
# One thread: a 4-thread run crashed on this miic branch (notebook 06)
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.log

set -euo pipefail
mkdir -p logs
# module load r/4.6   # or: conda activate <your R environment>

DIR=results/case2/06_rebuilt_network
Rscript scripts/run_miic.R "$DIR/input.tsv" "$DIR/state_order.tsv" "$DIR/interact_edges.tsv" "$DIR" "$SLURM_CPUS_PER_TASK"
