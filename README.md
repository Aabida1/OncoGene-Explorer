# OncoGene Explorer

An interactive, Streamlit-based cancer transcriptomics research platform for differential gene-expression exploration, independent cohort comparison, exploratory consensus analysis, pathway enrichment, immune-signature exploration, and gene-level evidence reporting.

## Current focus
Breast-cancer transcriptomics, starting with the supplied GSE183947 precomputed differential-expression results and a separate GSE42568 Affymetrix workflow. Each cohort must be analyzed using a platform-appropriate method; raw expression matrices from different platforms must not be concatenated as if directly comparable.

## Features
- Interactive cohort-specific volcano plot and DEG tables with filters and CSV exports.
- Upload additional cohort-level DEG result tables and select the active cohort.
- Compare gene-level log fold changes across independent cohorts.
- Exploratory Fisher combined-p-value and Benjamini–Hochberg correction summaries, with direction consistency and coverage reported separately.
- GMT pathway over-representation analysis.
- Immune-signature summaries and gene-level evidence bridge using user-supplied evidence.
- Analysis provenance and limitations documentation.

## Quick start
Python 3.10+ is recommended; Python 3.11 is tested in CI.

```bash
python -m venv .venv
# macOS/Linux
source .venv/bin/activate
# Windows PowerShell: .venv\\Scripts\\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```
Then open the local URL printed by Streamlit, usually http://localhost:8501.

## Deploy
### Streamlit Community Cloud
1. Select **New app** and connect this repository.
2. Choose branch `main` and entrypoint `app.py`.
3. Deploy; dependencies are in `requirements.txt`.

### Render
Build command: `pip install -r requirements.txt`

Start command: `streamlit run app.py --server.port $PORT --server.address 0.0.0.0`

## Data and analysis notes
- `data/limma_tumor_vs_normal_all_genes.csv` contains supplied GSE183947 precomputed results; it is not raw read data.
- `data/sample_metadata.csv` contains supplied sample annotations; confirm sample mapping and the paired design before making biological claims.
- `scripts/run_gse42568.R` runs a separate limma workflow for the GSE42568 Affymetrix cohort and expects its source inputs to be downloaded and placed as described in the script.
- `templates/deg_upload_template.csv` shows the minimum upload schema. Additional cohort result CSVs should contain one row per gene and columns equivalent to `gene_id`, `logFC`, `P.Value`, and `adj.P.Val`.

## Cross-cohort interpretation
Cohorts are analyzed separately. Identifier mapping, contrast direction, gene duplication and study design should be checked before comparison. A gene missing from a cohort is not assigned zero expression or a p-value of 1. Fisher's method combines p-values; it is not a pooled effect-size estimate or a formal random-effects meta-analysis. Scores and immune signatures are exploratory and do not establish biomarker validity, pathway activation, or cell abundance.

## Checks
GitHub Actions runs Python syntax compilation and checks that the key app/data files are present on pushes and pull requests. This CI smoke check does not replace running the app interactively against real user-uploaded data or testing a deployment.

## Disclaimer
For research and educational use only. Not intended for diagnosis, treatment selection, or other clinical decisions.
