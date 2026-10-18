# Benchmarks

Timings of MIIC and of the mutual information (MI) step, for the patched `miic` variants of `envs/patches/miic/`.
Read by notebooks `case2/10_performance.qmd` and `case1/05_performance_vs_server.qmd`; written by the scripts of
`scripts/benchmarks/` (run from the root of the repository, in WSL/Linux, one benchmark at a time).

| File or folder              | Content                                                                                                   | Written by                                                     |
| --------------------------- | --------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------- |
| `benchmarks.tsv`            | Summary of case study 2: one row per tested change, verdict (identical or not) and time                  | notebook `case2/10_performance.qmd`                            |
| `case1_server_vs_local.tsv` | Summary of case study 1: web server vs. local runs, stock vs. patched `miic`, 1 to 12 threads           | `compare_case1.R`                                              |
| `miic/<dataset>/<run>/`     | One MIIC run: `run.log`, `status.txt` (ok/failed), `time.txt` (GNU time: wall time, peak memory), `miic_path.txt` (the `miic` actually loaded) | `bench_miic.sh`, called by `run_benchmarks.sh` and `run_case1_benchmarks.sh` |
| `mi/`                       | One MI run per file (`<method>_<cell type>_b<block size>_r<repeat>`): `.log` and GNU time output        | `run_mi_benchmarks.sh`                                         |
| `queue.log`                 | Start and end of each run of the queue                                                                    | `run_*benchmarks.sh` (via `tee`)                               |
| `load.log`                  | Machine load every minute during the runs, to spot timings disturbed by other processes                  | `log_load.sh`                                                  |
| `stop_load_log`             | Empty flag file: stops `log_load.sh`                                                                      | by hand                                                        |

Run names: `<variant>_t<threads>`. Datasets: `case1` (authors' input of case study 1), `authors` (authors' input of
case study 2), `rebuilt` (our rebuilt input, notebook `case2/06`), `dedup` (authors' input without duplicated genes,
notebook `case2/07`). The MIIC result files (`.rds`) are not tracked by Git.
