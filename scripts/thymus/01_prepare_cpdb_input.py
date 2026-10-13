"""Subsample 1000 cells per cell type (as in Park et al. 2020) and write CellPhoneDB inputs.
Run in WSL: python 01_prepare_cpdb_input.py"""
import os
import anndata as ad, numpy as np, pandas as pd
ROOT = '/mnt/c/Users/Dell 5540/Documents/computational_biology_project/causalccc-reproduction'
LEVEL = 'Anno_level_fig1'
a = ad.read_h5ad(f'{ROOT}/data/thymus/HTA08_v01_A05_Science_human_fig1.h5ad', backed='r')
print('var_names head:', list(a.var_names[:8]), '| duplicated:', a.var_names.duplicated().sum())
rng = np.random.default_rng(0)
idx = []
for ct, g in a.obs.groupby(LEVEL, observed=True):
    pos = np.flatnonzero((a.obs[LEVEL] == ct).values)
    idx += list(rng.choice(pos, min(1000, len(pos)), replace=False))
idx = np.sort(idx)
s = a[idx].to_memory()
s = s[:, ~s.var_names.duplicated()].copy()
s.write_h5ad(f'{ROOT}/data/thymus/subsample_1000_allgenes.h5ad')   # full gene set, kept for the MIIC step
# CellPhoneDB densifies the matrix: keep only genes known to its database (result unchanged, memory x15 lower)
cpdb_genes = set(pd.read_csv(os.path.expanduser('~/cpdb_db/gene_input.csv'))['hgnc_symbol'].dropna())
s = s[:, s.var_names.isin(cpdb_genes)].copy()
print('genes kept for CellPhoneDB:', s.n_vars)
s.write_h5ad(f'{ROOT}/data/thymus/cpdb_input.h5ad')
pd.DataFrame({'Cell': s.obs_names, 'cell_type': s.obs[LEVEL].astype(str).values}).to_csv(
    f'{ROOT}/data/thymus/cpdb_meta.tsv', sep='\t', index=False)
print(s.shape); print(s.obs[LEVEL].value_counts().to_string())
