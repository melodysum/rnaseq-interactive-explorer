"""
app.py
------
RNA-seq Interactive Explorer — Streamlit app.

Loads real count data (counts.csv / metadata.csv) and lets the user
interactively explore how low-expression filtering and normalisation
choices affect differential expression results.

Run:
    streamlit run app.py
"""

import os
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from analysis_utils import (
    load_data,
    get_paired_columns,
    filter_low_expression,
    normalise,
    run_paired_de,
)

# ── Page setup ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="RNA-seq Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .main-header {
        font-size: 2rem;
        font-weight: 700;
        color: #1a1a2e;
        margin-bottom: 0.2rem;
    }
    .sub-caption {
        color: #6c757d;
        font-size: 0.9rem;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background: #f8f9fa;
        border-radius: 8px;
        padding: 0.8rem 1rem;
        border-left: 4px solid #4361ee;
    }
    .warning-box {
        background: #fff3cd;
        border: 1px solid #ffc107;
        border-radius: 6px;
        padding: 0.7rem 1rem;
        font-size: 0.85rem;
    }
    div[data-testid="stTabs"] button {
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

# ── Load data ─────────────────────────────────────────────────────────────────
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

@st.cache_data(show_spinner="Loading count matrix…")
def cached_load():
    counts, metadata = load_data(
        os.path.join(DATA_DIR, "counts.csv"),
        os.path.join(DATA_DIR, "metadata.csv"),
    )
    donors, ctrl_cols, treat_cols = get_paired_columns(metadata)
    return counts, metadata, donors, ctrl_cols, treat_cols

counts_raw, metadata, donors, ctrl_cols, treat_cols = cached_load()
N_GENES_TOTAL = len(counts_raw)
N_DONORS      = len(donors)

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown('<div class="main-header">🧬 RNA-seq Interactive Explorer</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="sub-caption">'
    f'Dataset: <b>{N_GENES_TOTAL:,} genes</b> × <b>{N_DONORS} donors</b> '
    f'(paired control / treatment design) — 40 samples total'
    '</div>',
    unsafe_allow_html=True,
)

# ── Sidebar controls ──────────────────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Analysis Parameters")

    st.subheader("1 · Low-expression filtering")
    min_count = st.slider(
        "Minimum count per sample",
        min_value=0, max_value=50, value=10, step=1,
        help="A gene must reach this count in at least N samples to be kept.",
    )
    min_samples = st.slider(
        "Minimum samples passing threshold",
        min_value=1, max_value=N_DONORS, value=5, step=1,
        help="How many samples must exceed the count threshold.",
    )

    st.subheader("2 · Normalisation method")
    norm_method = st.selectbox(
        "Method",
        ["log-CPM", "CPM", "Raw counts", "Size-factor"],
        help="Affects distribution plots. DE testing always uses log-CPM internally.",
    )

    st.subheader("3 · DE significance thresholds")
    fdr_cutoff = st.slider(
        "FDR cutoff (adjusted p-value)",
        min_value=0.01, max_value=0.20, value=0.05, step=0.01,
    )
    lfc_cutoff = st.slider(
        "Minimum |log₂FC|",
        min_value=0.0, max_value=3.0, value=1.0, step=0.1,
    )

    st.divider()
    st.caption(
        "⚠️ **Educational demo only.** DE uses a paired t-test on log-CPM "
        "values with BH correction. A production analysis would use "
        "DESeq2 or edgeR."
    )

# ── Run analysis ──────────────────────────────────────────────────────────────
@st.cache_data(show_spinner="Running analysis…")
def run_analysis(min_count, min_samples, norm_method, fdr_cutoff, lfc_cutoff):
    filtered  = filter_low_expression(counts_raw, min_count, min_samples)
    normed    = normalise(filtered, norm_method)
    de_result = run_paired_de(filtered, ctrl_cols, treat_cols, fdr_cutoff, lfc_cutoff)
    return filtered, normed, de_result

counts_filtered, counts_norm, de_results = run_analysis(
    min_count, min_samples, norm_method, fdr_cutoff, lfc_cutoff
)

n_after = len(counts_filtered)
n_removed = N_GENES_TOTAL - n_after
n_up   = (de_results['direction'] == 'Up in treatment').sum()
n_down = (de_results['direction'] == 'Down in treatment').sum()
n_sig  = n_up + n_down

# ── Summary metrics ───────────────────────────────────────────────────────────
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Total genes", f"{N_GENES_TOTAL:,}")
c2.metric("After filtering", f"{n_after:,}", delta=f"−{n_removed:,} removed")
c3.metric("Significant DE", f"{n_sig:,}")
c4.metric("⬆ Up in treatment", f"{n_up:,}")
c5.metric("⬇ Down in treatment", f"{n_down:,}")

st.divider()

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab_dist, tab_volcano, tab_ma, tab_table, tab_info = st.tabs([
    "📊 Distributions",
    "🌋 Volcano Plot",
    "📈 MA Plot",
    "📋 Top DE Genes",
    "ℹ️ Explanations",
])

# ─── Tab 1 · Distributions ───────────────────────────────────────────────────
with tab_dist:
    col_left, col_right = st.columns(2)

    with col_left:
        st.subheader("Expression distribution: before vs after filtering")
        # Sample genes for speed (plotting all ~10k is slow in browser)
        rng = np.random.default_rng(42)
        idx_before = rng.choice(len(counts_raw),    size=min(3000, len(counts_raw)),    replace=False)
        idx_after  = rng.choice(len(counts_filtered), size=min(3000, len(counts_filtered)), replace=False)

        vals_before = np.log2(counts_raw.iloc[idx_before].values.flatten() + 1)
        vals_after  = np.log2(counts_filtered.iloc[idx_after].values.flatten() + 1)

        fig_hist = go.Figure()
        fig_hist.add_trace(go.Histogram(
            x=vals_before, name="Before filtering",
            opacity=0.55, nbinsx=80, marker_color="#4361ee",
        ))
        fig_hist.add_trace(go.Histogram(
            x=vals_after, name="After filtering",
            opacity=0.55, nbinsx=80, marker_color="#e63946",
        ))
        fig_hist.update_layout(
            barmode="overlay",
            xaxis_title="log₂(count + 1)",
            yaxis_title="Frequency",
            height=370,
            legend=dict(x=0.65, y=0.95),
            margin=dict(t=20),
        )
        st.plotly_chart(fig_hist, use_container_width=True)
        st.caption(
            "Low-expression filtering removes the left-hand spike of near-zero counts. "
            "These genes add noise without contributing reliable signal."
        )

    with col_right:
        st.subheader(f"Per-sample distributions ({norm_method})")

        # Show first 5 pairs for clarity
        sample_cols = ctrl_cols[:5] + treat_cols[:5]
        display = counts_norm[sample_cols]

        # Always plot in log space for readability
        if norm_method in ("Raw counts", "CPM", "Size-factor"):
            plot_vals = np.log2(display.replace(0, np.nan) + 1)
            y_label   = f"log₂({norm_method} + 1)"
        else:
            plot_vals = display
            y_label   = "log₂(CPM + 1)"

        colors = ["#4361ee"] * 5 + ["#e63946"] * 5
        fig_box = go.Figure()
        for i, col in enumerate(sample_cols):
            group = "Control" if "control" in col else "Treatment"
            fig_box.add_trace(go.Violin(
                y=plot_vals[col].dropna(),
                name=col.replace("_", " "),
                box_visible=True,
                meanline_visible=True,
                marker_color=colors[i],
                legendgroup=group,
                showlegend=(i in (0, 5)),
                legendgrouptitle_text=group if i in (0, 5) else None,
            ))
        fig_box.update_layout(
            yaxis_title=y_label,
            height=370,
            margin=dict(t=20),
            violinmode="overlay",
        )
        st.plotly_chart(fig_box, use_container_width=True)
        st.caption(
            "After normalisation, sample distributions should be roughly aligned. "
            "Shift normalisation method in the sidebar to see the effect."
        )

# ─── Tab 2 · Volcano plot ─────────────────────────────────────────────────────
with tab_volcano:
    st.subheader("Volcano Plot — Treatment vs Control (paired t-test, BH FDR)")

    de_plot = de_results.reset_index()
    color_map = {
        "Not significant":    "#adb5bd",
        "Up in treatment":    "#e63946",
        "Down in treatment":  "#4361ee",
    }
    # Cap y-axis for display
    y_cap = min(de_plot["neg_log10_padj"].quantile(0.999) * 1.1, 50)

    fig_vol = px.scatter(
        de_plot,
        x="log2FC",
        y="neg_log10_padj",
        color="direction",
        color_discrete_map=color_map,
        hover_name="gene",
        hover_data={
            "log2FC":          ":.3f",
            "pvalue":          ":.2e",
            "padj":            ":.2e",
            "direction":       False,
            "neg_log10_padj":  False,
        },
        labels={
            "log2FC":         "log₂ Fold Change (Treatment / Control)",
            "neg_log10_padj": "−log₁₀(adjusted p-value)",
        },
        category_orders={"direction": ["Up in treatment", "Down in treatment", "Not significant"]},
        height=520,
        opacity=0.65,
    )
    fig_vol.update_traces(marker=dict(size=4))
    fig_vol.add_vline(x= lfc_cutoff, line_dash="dot", line_color="#6c757d", opacity=0.7)
    fig_vol.add_vline(x=-lfc_cutoff, line_dash="dot", line_color="#6c757d", opacity=0.7)
    fig_vol.add_hline(y=-np.log10(fdr_cutoff), line_dash="dot", line_color="#6c757d", opacity=0.7)
    fig_vol.update_layout(
        yaxis_range=[0, y_cap],
        legend_title="Direction",
        margin=dict(t=20),
    )
    st.plotly_chart(fig_vol, use_container_width=True)
    st.info(
        f"Dashed lines mark FDR < {fdr_cutoff} (horizontal) and |log₂FC| > {lfc_cutoff} (vertical). "
        f"Adjust thresholds in the sidebar to see how the significant set changes."
    )

# ─── Tab 3 · MA plot ──────────────────────────────────────────────────────────
with tab_ma:
    st.subheader("MA Plot — Mean expression vs fold change")
    st.caption(
        "An MA plot reveals whether fold changes are biased by expression level. "
        "Ideally the cloud should be centered on log₂FC = 0 across all expression levels."
    )

    de_ma = de_results.reset_index()
    fig_ma = px.scatter(
        de_ma,
        x="mean_expr",
        y="log2FC",
        color="direction",
        color_discrete_map=color_map,
        hover_name="gene",
        hover_data={"log2FC": ":.3f", "padj": ":.2e", "direction": False},
        labels={
            "mean_expr": "Mean log₂(CPM + 1)",
            "log2FC":    "log₂ Fold Change",
        },
        category_orders={"direction": ["Up in treatment", "Down in treatment", "Not significant"]},
        height=480,
        opacity=0.55,
    )
    fig_ma.update_traces(marker=dict(size=3))
    fig_ma.add_hline(y=0, line_dash="solid", line_color="black", opacity=0.3)
    fig_ma.add_hline(y= lfc_cutoff, line_dash="dot", line_color="#6c757d", opacity=0.6)
    fig_ma.add_hline(y=-lfc_cutoff, line_dash="dot", line_color="#6c757d", opacity=0.6)
    fig_ma.update_layout(margin=dict(t=20))
    st.plotly_chart(fig_ma, use_container_width=True)

# ─── Tab 4 · Top DE genes ─────────────────────────────────────────────────────
with tab_table:
    st.subheader("Top DE genes (sorted by adjusted p-value)")

    sig_genes = de_results[de_results["significant"]].sort_values("padj")

    if len(sig_genes) == 0:
        st.warning(
            "No significant DE genes found with the current thresholds. "
            "Try lowering the FDR cutoff or the |log₂FC| threshold in the sidebar."
        )
    else:
        top_n = st.slider("Show top N genes", 10, min(200, len(sig_genes)), 30, 10)
        display_df = sig_genes.head(top_n)[["log2FC", "mean_expr", "t_stat", "pvalue", "padj", "direction"]]
        display_df = display_df.rename(columns={
            "log2FC":    "log₂FC",
            "mean_expr": "Mean log-CPM",
            "t_stat":    "t statistic",
            "pvalue":    "p-value",
            "padj":      "adj. p-value",
            "direction": "Direction",
        })
        st.dataframe(
            display_df.style.format({
                "log₂FC":        "{:.3f}",
                "Mean log-CPM":  "{:.2f}",
                "t statistic":   "{:.2f}",
                "p-value":       "{:.2e}",
                "adj. p-value":  "{:.2e}",
            }),
            use_container_width=True,
            height=420,
        )

        st.download_button(
            "⬇️ Download full DE results (CSV)",
            data=de_results.reset_index().to_csv(index=False),
            file_name="de_results.csv",
            mime="text/csv",
        )

# ─── Tab 5 · Explanations ─────────────────────────────────────────────────────
with tab_info:
    st.subheader("📚 Conceptual guide")

    with st.expander("🔍 Low-expression filtering — what and why", expanded=True):
        st.markdown(f"""
**What it does:**
Removes genes whose count never (or rarely) exceeds a minimum threshold across samples.
With the current settings, a gene must have ≥ **{min_count} counts** in at least
**{min_samples} samples** to be retained.

**Why it matters:**
- A count of 0 vs 1 vs 2 is dominated by sampling noise, not biology.
- These genes inflate the total number of tests performed. After FDR correction,
  every extra test dilutes your statistical power.
- They also distort fold-change estimates: a true count of 0 in one condition and 2 in another
  gives an enormous but meaningless fold change.

**Right now:** {N_GENES_TOTAL:,} → {n_after:,} genes kept ({n_removed:,} removed).
        """)

    with st.expander("⚖️ Normalisation — making samples comparable"):
        st.markdown("""
**The problem:**
RNA-seq libraries have different total read counts (library sizes). A gene with 100 counts
in a 5M-read library is *more* expressed than the same count in a 20M-read library.

**Methods available here:**

| Method | Formula | Notes |
|---|---|---|
| Raw counts | none | Only useful to see the raw numbers |
| CPM | count / lib_size × 10⁶ | Removes library-size bias |
| log-CPM | log₂(CPM + 1) | Also stabilises variance; used for DE testing |
| Size-factor | count / median ratio | Concept from DESeq2; robust to outlier genes |

**Important:** Regardless of which display normalisation you choose in the sidebar,
the DE analysis always uses **log-CPM** internally. This is the standard approach
for t-test-based methods.
        """)

    with st.expander("📐 Paired design & the t-test"):
        st.markdown(f"""
**Why paired?**
Each of the {N_DONORS} donors contributes one control sample and one treatment sample.
Donor-to-donor expression differences (genetics, age, sex, environment) are often much
larger than the treatment effect. A paired test subtracts this donor noise before
asking "does treatment change expression?"

**How it works here:**
1. Compute log₂CPM for each sample.
2. For each gene, calculate the *difference* (treatment − control) per donor.
3. Run a one-sample t-test on these {N_DONORS} differences (equivalent to a paired t-test).
4. Apply Benjamini-Hochberg FDR correction across all genes.

This is a simplified demonstration. Real paired RNA-seq analysis uses DESeq2 or limma-voom,
which model count-specific variance and handle batch effects more rigorously.
        """)

    with st.expander("🎯 FDR — why p-values alone are not enough"):
        st.markdown(f"""
**The multiple testing problem:**
With ~{n_after:,} genes tested, even a p-value threshold of 0.05 would give
~{int(n_after * 0.05):,} false positives by chance alone.

**Benjamini-Hochberg FDR:**
Controls the *expected proportion* of false positives among significant results.
An FDR of {fdr_cutoff} means that, on average, {int(fdr_cutoff*100)}% of your
significant genes are false discoveries.

**The interaction with filtering:**
Removing low-expression genes *before* the test reduces the number of hypotheses tested.
This makes the FDR correction less severe and increases power to detect real DE genes.
Try raising the filtering thresholds and watch what happens to the significant gene count.
        """)

    with st.expander("⚠️ Limitations of this demo"):
        st.markdown("""
- **Not DESeq2/edgeR.** Real tools use negative binomial models that properly
  handle the count-specific variance of RNA-seq data.
- **No batch correction.** The metadata includes batch (batch1/batch2).
  A real analysis would model or correct for batch.
- **No outlier/QC checks.** Sample-level QC (PCA, library size checks) is skipped.
- **t-test assumptions.** A parametric test on log-CPM is a reasonable approximation
  but can be anti-conservative at very low counts.
- **Synthetic-looking genes.** The count values here follow patterns typical of real
  RNA-seq but the gene-to-function interpretation still requires proper context.
        """)
