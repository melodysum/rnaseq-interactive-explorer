# RNA-seq Interactive Explorer

An educational Streamlit app for exploring how **low-expression filtering** and
**normalisation** choices affect differential expression results in a real paired
RNA-seq dataset.

---

## Dataset

| File | Description |
|---|---|
| `data/counts.csv` | 9,936 genes × 40 samples (raw counts) |
| `data/metadata.csv` | Sample metadata: donor, group, batch, sex, age |

**Design:** 20 donors, each with one control and one treatment sample (paired).

---

## Setup

```bash
# 1. Create a conda environment (recommended)
conda create -n rnaseq-app python=3.11
conda activate rnaseq-app

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run the app
streamlit run app.py
```

The app opens automatically at `http://localhost:8501`.

---

## Project structure

```
rnaseq_app/
├── app.py                 # Main Streamlit application
├── analysis_utils.py      # Filtering, normalisation, DE functions
├── requirements.txt
├── README.md
└── data/
    ├── counts.csv
    └── metadata.csv
```

---

## What the app does

| Sidebar control | What changes |
|---|---|
| Min count threshold | Fewer/more low-expression genes removed |
| Min samples threshold | Stringency of filtering |
| Normalisation method | Distribution plots (DE always uses log-CPM) |
| FDR cutoff | Significance threshold for DE calls |
| Min \|log₂FC\| | Effect-size filter for DE calls |

### Tabs
- **Distributions** — histograms and violin plots before/after filtering
- **Volcano Plot** — −log₁₀(FDR) vs log₂FC with interactive hover
- **MA Plot** — fold change vs mean expression level
- **Top DE Genes** — sortable table + CSV download
- **Explanations** — conceptual guide for each analysis step

---

## Analysis method

- **Normalisation:** log₂(CPM + 1) for DE testing; CPM, raw, and size-factor available for display
- **DE test:** Paired t-test per gene (treatment vs control across donors)
- **Multiple testing:** Benjamini-Hochberg FDR correction
- **Significance:** FDR < threshold AND |log₂FC| ≥ threshold

> ⚠️ This is a simplified educational demo. A production analysis should use
> DESeq2 or limma-voom with proper dispersion modelling and batch correction.

---

## Deploying to Streamlit Cloud (share with supervisor)

1. Push this folder to a GitHub repository.
2. Go to [share.streamlit.io](https://share.streamlit.io) and connect the repo.
3. Set **Main file path** to `app.py`.
4. Deploy — you get a public URL in ~2 minutes.
