# OncoGene Explorer

An interactive Streamlit research dashboard for cancer transcriptomics. It supports differential-expression exploration, independent cohort comparison, exploratory cross-cohort evidence aggregation, GMT pathway enrichment, immune-signature screening, and source-traceable target-evidence matching.

## Quick start

```bash
python -m venv .venv
# macOS/Linux
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```

## Repository structure

- `app.py` — Streamlit dashboard
- `data/analysis_provenance.txt` — current data-bundling status and provenance notes
- `scripts/run_gse42568.R` — GEOquery + limma workflow for GSE42568 processed microarray data
- `templates/deg_upload_template.csv` — header-only input template
- `.github/workflows/ci.yml` — syntax and required-file smoke checks

## Add cohort data

The hosted app intentionally starts without cohort results because this repository does not yet bundle verified cohort DEG tables. Use **Study setup → Add cohort DEG CSV files** to upload one CSV per cohort.

Required fields (aliases may be supported by the app):

| Column | Description |
| --- | --- |
| `gene_id` | Gene symbol or other consistent gene identifier |
| `logFC` | Log fold change for the specified contrast |
| `P.Value` | Raw p-value |
| `adj.P.Val` | Multiple-testing-adjusted p-value |

One row should represent one gene. Confirm contrast direction, gene identifier conventions, duplicate mappings, and source provenance before uploading. The template contains only column headers; it does not contain example or fabricated gene results.

## GSE42568 workflow

The R script downloads the processed series matrix and GPL570 annotation from NCBI GEO, checks the expected sample-group counts (104 tumor, 17 normal), maps probes to gene symbols, aggregates multiple probes per gene by sample-wise median, and fits a limma tumor-versus-normal contrast with empirical Bayes moderation and Benjamini–Hochberg adjustment.

Install R and Bioconductor dependencies before running:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "Biobase"))
```

Then run from the repository root:

```bash
Rscript scripts/run_gse42568.R
```

Outputs are written to `data/GSE42568/`, including the all-genes DEG table, FDR < 0.05 table, sample metadata, and provenance notes. The script deliberately stops if sample labels or group counts do not match expectations; inspect its metadata-review output instead of overriding the check. This uses the GEO processed series matrix, not raw CEL reprocessing.

## Datasets and methods

- **GSE183947**: intended RNA-seq starting cohort; its verified DEG results are not bundled in this repository yet.
- **GSE42568**: Affymetrix GPL570 microarray cohort; analyze separately and harmonize identifiers before cross-cohort comparison.
- **Potential future cohort**: TCGA-BRCA with a documented tumor-versus-normal contrast and data version.

The app displays user-supplied precomputed results. It does not download GEO/TCGA data, process FASTQ files, or fit per-cohort differential-expression models at runtime. The optional pathway module expects a user-provided GMT gene-set file. The target module only matches user-supplied evidence and does not invent target/drug associations.

## Interpretation and limitations

Cohorts are analyzed separately. Missing gene mappings remain missing rather than being treated as zero effects or p-values of 1. Fisher's method combines raw p-values with BH correction across genes; it is not a formal effect-size meta-analysis. Direction consistency, candidate-prioritization scores, and immune-signature summaries are exploratory and do not validate biomarkers, establish causality, demonstrate pathway activation, or estimate cell abundance. Confirm sample metadata, pairing, contrast direction, normalization, gene identifiers, and source provenance before interpretation.

## Deploy

### Streamlit Community Cloud

Configure the app with repository `Aabida1/OncoGene-Explorer`, branch `main`, and entrypoint `app.py`. After deployment, upload verified cohort DEG CSVs in the sidebar; no cohort data are currently bundled.

### Render

Build command: `pip install -r requirements.txt`

Start command: `streamlit run app.py --server.port $PORT --server.address 0.0.0.0`

## Checks

GitHub Actions checks Python syntax and required project files. A passing CI run does not replace interactive testing or scientific validation of uploaded data.

## Disclaimer

Research and educational use only. Not for diagnosis, prognosis, or treatment decisions.
