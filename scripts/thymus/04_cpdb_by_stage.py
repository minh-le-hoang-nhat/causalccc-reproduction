"""CellPhoneDB with cell types split by developmental stage (as fig. S11 / Table S5).
Groups = Anno_level_fig1 x Anno_stage (early 7-11 pcw, middle 12 pcw-10 m, late >=15 m),
groups with < 50 cells dropped (fig. S11 rule), <= 1000 cells sampled per group (paper protocol).
One CellPhoneDB run per stage: only same-stage couples are interpreted, and a single run with all
99 groups (9,801 couples) exhausted the 7.8 GB of the WSL machine.
Run from the root of the repository: python scripts/thymus/04_cpdb_by_stage.py [stage ...]"""
from pathlib import Path
import os, sys
import anndata as ad, numpy as np, pandas as pd
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

ROOT = Path(__file__).resolve().parents[2]   # root of the repository
DB = os.path.expanduser('~/cpdb_db')
LEVEL, MIN_CELLS, MAX_CELLS = 'Anno_level_fig1', 50, 1000
STAGES = sys.argv[1:] or ['early', 'middle', 'late']

a = ad.read_h5ad(f'{ROOT}/data/thymus/HTA08_v01_A05_Science_human_fig1.h5ad', backed='r')
stage = a.obs['Anno_stage'].astype(str).str.rsplit('_', n=1).str[-1]
group = a.obs[LEVEL].astype(str) + '_' + stage
sizes = group.value_counts()
keep = sizes[sizes >= MIN_CELLS].index
print(f'{len(keep)} groups kept, {len(sizes) - len(keep)} dropped (<{MIN_CELLS} cells):',
      sorted(sizes[sizes < MIN_CELLS].index))
cpdb_genes = set(pd.read_csv(f'{DB}/gene_input.csv')['hgnc_symbol'].dropna())

for st in STAGES:
    rng = np.random.default_rng(0)
    idx = []
    for g in sorted(k for k in keep if k.rsplit('_', 1)[1] == st):
        pos = np.flatnonzero((group == g).values)
        idx += list(rng.choice(pos, min(MAX_CELLS, len(pos)), replace=False))
    idx = np.sort(idx)
    s = a[idx].to_memory()
    s = s[:, ~s.var_names.duplicated()]
    s = s[:, s.var_names.isin(cpdb_genes)].copy()
    s.obs['group'] = group.values[idx]
    print(st, 'input:', s.shape, s.obs['group'].nunique(), 'groups', flush=True)
    s.write_h5ad(f'{ROOT}/data/thymus/cpdb_{st}_input.h5ad')
    pd.DataFrame({'Cell': s.obs_names, 'cell_type': s.obs['group'].values}).to_csv(
        f'{ROOT}/data/thymus/cpdb_{st}_meta.tsv', sep='\t', index=False)
    del s

    out = f'{ROOT}/results/thymus/cpdb_by_stage/{st}'
    os.makedirs(out, exist_ok=True)
    cpdb_statistical_analysis_method.call(
        cpdb_file_path=f'{DB}/cellphonedb.zip',
        meta_file_path=f'{ROOT}/data/thymus/cpdb_{st}_meta.tsv',
        counts_file_path=f'{ROOT}/data/thymus/cpdb_{st}_input.h5ad',
        counts_data='hgnc_symbol', iterations=1000, threshold=0.1, threads=6,
        output_path=out, output_suffix=st, score_interactions=False, debug_seed=42)
    print(st, 'DONE', flush=True)
