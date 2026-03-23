"""
analysis_utils.py
-----------------
Core analysis functions for the RNA-seq interactive explorer.
Handles filtering, normalisation, and differential expression.

Dataset: paired design (control vs treatment, same donors).
DE method: paired t-test on log-CPM values + Benjamini-Hochberg FDR.
This is a simplified educational demonstration — not a substitute for
DESeq2/edgeR in a real analysis.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


# ── Data loading ─────────────────────────────────────────────────────────────

def load_data(counts_path: str, metadata_path: str):
    """Load count matrix and metadata from CSV files.

    Returns
    -------
    counts : pd.DataFrame, shape (genes, samples)
    metadata : pd.DataFrame, shape (samples, features)
    """
    counts = pd.read_csv(counts_path, index_col=0)
    metadata = pd.read_csv(metadata_path, index_col=0)
    return counts, metadata


def get_paired_columns(metadata: pd.DataFrame):
    """Extract ordered control and treatment column lists from metadata.

    Preserves donor pairing so index i in ctrl_cols matches index i in treat_cols.
    """
    donors = metadata['donor'].unique()
    ctrl_cols  = [f"{d}_control"   for d in donors]
    treat_cols = [f"{d}_treatment" for d in donors]
    return donors, ctrl_cols, treat_cols


# ── Filtering ────────────────────────────────────────────────────────────────

def filter_low_expression(counts: pd.DataFrame,
                           min_count: int = 10,
                           min_samples: int = 5) -> pd.DataFrame:
    """Remove genes that fail to reach min_count in at least min_samples samples.

    Parameters
    ----------
    counts      : raw count matrix (genes × samples)
    min_count   : minimum count per sample to be considered 'expressed'
    min_samples : minimum number of samples that must pass the count threshold

    Returns
    -------
    Filtered DataFrame (subset of rows).
    """
    mask = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[mask]


# ── Normalisation ────────────────────────────────────────────────────────────

def normalise(counts: pd.DataFrame, method: str = "log-CPM") -> pd.DataFrame:
    """Normalise count matrix.

    Methods
    -------
    'Raw counts'  : no transformation (returns float copy)
    'CPM'         : counts per million (library-size normalised)
    'log-CPM'     : log2(CPM + 1)  — default; used for DE testing
    'Size-factor' : simple median-ratio normalisation (like DESeq2 concept)
    """
    lib_sizes = counts.sum(axis=0)

    if method == "Raw counts":
        return counts.astype(float)

    cpm = counts.divide(lib_sizes, axis=1) * 1e6

    if method == "CPM":
        return cpm

    if method == "log-CPM":
        return np.log2(cpm + 1)

    if method == "Size-factor":
        # Geometric mean per gene (ignoring zeros)
        counts_nz = counts.replace(0, np.nan)
        geo_means = np.exp(np.log(counts_nz).mean(axis=1))
        ratios = counts.divide(geo_means, axis=0)
        size_factors = ratios.median(axis=0)
        return counts.divide(size_factors, axis=1)

    raise ValueError(f"Unknown normalisation method: {method}")


# ── Differential expression ───────────────────────────────────────────────────

def run_paired_de(counts_filtered: pd.DataFrame,
                  ctrl_cols: list,
                  treat_cols: list,
                  fdr_cutoff: float = 0.05,
                  lfc_cutoff: float = 1.0) -> pd.DataFrame:
    """Run paired DE analysis using log-CPM + paired t-test + BH correction.

    The paired t-test accounts for donor-to-donor variability.
    log-CPM is always used for the test regardless of what the user
    chose for the display normalisation — this is the standard approach.

    Returns
    -------
    pd.DataFrame with columns:
        gene, log2FC, mean_expr, t_stat, pvalue, padj, significant
    """
    # Always test in log-CPM space
    lib_sizes = counts_filtered.sum(axis=0)
    cpm = counts_filtered.divide(lib_sizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    ctrl_mat  = log_cpm[ctrl_cols].values   # shape: (genes, donors)
    treat_mat = log_cpm[treat_cols].values

    # Mean log2 fold change (paired difference in log space = log ratio)
    diff = treat_mat - ctrl_mat
    log2fc = diff.mean(axis=1)

    # Mean expression (average of control and treatment, log-CPM)
    mean_expr = ((ctrl_mat + treat_mat) / 2).mean(axis=1)

    # Paired t-test per gene
    t_stats, p_values = stats.ttest_rel(treat_mat, ctrl_mat, axis=1)

    # Benjamini-Hochberg FDR correction
    _, padj, _, _ = multipletests(p_values, method='fdr_bh')

    results = pd.DataFrame({
        'gene':      counts_filtered.index,
        'log2FC':    log2fc,
        'mean_expr': mean_expr,
        't_stat':    t_stats,
        'pvalue':    p_values,
        'padj':      padj,
    })

    results['significant'] = (
        (results['padj'] < fdr_cutoff) &
        (results['log2FC'].abs() >= lfc_cutoff)
    )
    results['neg_log10_padj'] = -np.log10(results['padj'].clip(lower=1e-300))
    results['direction'] = 'Not significant'
    results.loc[
        results['significant'] & (results['log2FC'] > 0), 'direction'
    ] = 'Up in treatment'
    results.loc[
        results['significant'] & (results['log2FC'] < 0), 'direction'
    ] = 'Down in treatment'

    return results.set_index('gene')
