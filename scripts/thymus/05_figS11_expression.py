"""Reproduce fig. S11: expression of thymic-development signalling genes per cell type x stage.
Groups = Anno_level_fig1 x Anno_stage, groups with < 50 cells dropped (fig. S11 rule).
Statistics on all cells (no subsampling): mean log-normalised expression and fraction of cells > 0.
Run in WSL: python 05_figS11_expression.py"""
import os
import anndata as ad, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = '/mnt/c/Users/Dell 5540/Documents/computational_biology_project/causalccc-reproduction'
OUT = f'{ROOT}/results/thymus/figS11'
os.makedirs(OUT, exist_ok=True)
GENES = ['LTB', 'LTBR', 'TNFSF11', 'TNFRSF11A', 'FGF7', 'FGFR2', 'NOTCH1', 'JAG1', 'JAG2', 'DLL4']
STAGES = ['early', 'middle', 'late']

a = ad.read_h5ad(f'{ROOT}/data/thymus/HTA08_v01_A05_Science_human_fig1.h5ad', backed='r')
gidx = [int(np.flatnonzero(a.var_names == g)[0]) for g in GENES]
cell_type = a.obs['Anno_level_fig1'].astype(str).values
stage = a.obs['Anno_stage'].astype(str).str.rsplit('_', n=1).str[-1].values
group = pd.Series(cell_type) + '_' + pd.Series(stage)

# stream the matrix by row blocks, keep only the 10 genes
blocks = []
for i in range(0, a.n_obs, 20000):
    blocks.append(a.X[i:i + 20000][:, gidx].toarray())
X = pd.DataFrame(np.vstack(blocks), columns=GENES)
X['group'], X['cell_type'], X['stage'] = group.values, cell_type, stage

n = X.groupby('group').size()
mean = X.groupby('group')[GENES].mean()
frac = X[GENES].gt(0).groupby(X['group']).mean()
keep = n[n >= 50].index
long = (mean.loc[keep].stack().rename('mean').to_frame()
        .join(frac.loc[keep].stack().rename('frac')).reset_index()
        .rename(columns={'level_1': 'gene'}))
long['n_cells'] = long['group'].map(n)
long['cell_type'] = long['group'].str.rsplit('_', n=1).str[0]
long['stage'] = long['group'].str.rsplit('_', n=1).str[-1]
long.to_csv(f'{OUT}/expression_by_type_stage.tsv', sep='\t', index=False)
print(f'{len(keep)} groups kept (>=50 cells) of {len(n)}')

# dot plot: rows = cell type x stage, ordered by type then stage;
# colour = mean scaled to the gene's maximum (scanpy standard_scale='var'), so weak genes stay visible
order = sorted(keep, key=lambda g: (g.rsplit('_', 1)[0], STAGES.index(g.rsplit('_', 1)[1])))
long['scaled'] = long['mean'] / long.groupby('gene')['mean'].transform('max')
fig, ax = plt.subplots(figsize=(6.5, 0.2 * len(order) + 1.5))
for j, gname in enumerate(GENES):
    sub = long[long.gene == gname].set_index('group').loc[order]
    sc = ax.scatter(np.full(len(order), j), np.arange(len(order)), s=sub['frac'] * 120,
                    c=sub['scaled'], cmap='Reds', vmin=0, vmax=1, edgecolors='none')
ax.set_xticks(range(len(GENES)), GENES, rotation=90, style='italic')
ax.set_yticks(range(len(order)), order, fontsize=6.5)
ax.set_ylim(len(order) - 0.5, -0.5)
ax.set_xlim(-0.6, len(GENES) - 0.4)
fig.colorbar(sc, ax=ax, shrink=0.3, label='mean expression / gene max')
for f in (0.25, 0.5, 1.0):
    ax.scatter([], [], s=f * 120, c='grey', label=f'{int(f * 100)} %')
ax.legend(title='cells > 0', loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(f'{OUT}/figS11_dotplot.png', dpi=150)
print('saved', f'{OUT}/figS11_dotplot.png')

# "NOTCH1 is the main receptor in ETPs": compare NOTCH1-4 in ETP (and DN) per stage
NG = ['NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4']
rows = np.flatnonzero(np.isin(cell_type, ['ETP', 'DN']))
ngi = [int(np.flatnonzero(a.var_names == g)[0]) for g in NG]
N = pd.DataFrame(a.X[rows][:, ngi].toarray(), columns=NG)
N['group'] = group.values[rows]
notch = (N.groupby('group')[NG].mean().stack().rename('mean').to_frame()
         .join(N[NG].gt(0).groupby(N['group']).mean().stack().rename('frac')).reset_index()
         .rename(columns={'level_1': 'gene'}))
notch.to_csv(f'{OUT}/notch_receptors_ETP_DN.tsv', sep='\t', index=False)

# top groups per gene (used for the claim checks)
for gname in GENES:
    t = long[long.gene == gname].sort_values('mean', ascending=False).head(8)
    print(f'\n## {gname}')
    print(t[['group', 'n_cells', 'mean', 'frac']].round(3).to_string(index=False))
