"""CellPhoneDB v5 statistical analysis on the subsampled thymus data (Park et al. 2020 protocol).
Run from the root of the repository: python scripts/thymus/02_run_cellphonedb.py"""
from pathlib import Path
import os
from cellphonedb.utils import db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
ROOT = Path(__file__).resolve().parents[2]   # root of the repository
db_dir = os.path.expanduser('~/cpdb_db'); os.makedirs(db_dir, exist_ok=True)
db_utils.download_database(db_dir, 'v5.0.0')
out = f'{ROOT}/results/thymus/cpdb'; os.makedirs(out, exist_ok=True)
cpdb_statistical_analysis_method.call(
    cpdb_file_path=f'{db_dir}/cellphonedb.zip',
    meta_file_path=f'{ROOT}/data/thymus/cpdb_meta.tsv',
    counts_file_path=f'{ROOT}/data/thymus/cpdb_input.h5ad',
    counts_data='hgnc_symbol', iterations=1000, threshold=0.1, threads=10,
    output_path=out, output_suffix='thymus', score_interactions=False, debug_seed=42)
print('DONE')
