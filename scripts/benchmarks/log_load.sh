#!/usr/bin/env bash
# Logs the machine load every minute while benchmarks run, so that timings biased by other processes can be spotted
# (notebook 10). Stops when results/benchmarks/stop_load_log exists.
# Usage (from the root of the repository): bash scripts/benchmarks/log_load.sh &
# Each line of results/benchmarks/load.log: date, time, 1/5/15-minute load averages, number of R processes, and how
# many of them are not benchmark runs (e.g. a notebook being rendered).
log=results/benchmarks/load.log
stop=results/benchmarks/stop_load_log
rm -f "$stop"
while [[ ! -f "$stop" ]]; do
  r_all=$(pgrep -c -f 'exec/R')
  r_other=$(pgrep -a -f 'exec/R' | grep -v -c -E 'run_miic.R|bench_mi.R')
  echo "$(date '+%F %T') $(cut -d ' ' -f 1-3 /proc/loadavg) R_processes=$r_all other_R=$r_other" >> "$log"
  sleep 60
done
