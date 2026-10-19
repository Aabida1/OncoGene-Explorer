# OncoGene Explorer

An interactive Streamlit research dashboard for cancer transcriptomics. It supports differential-expression exploration, independent cohort comparison, exploratory cross-cohort evidence aggregation, GMT pathway enrichment, immune-signature screening, and source-traceable target-evidence matching.

## Quick start
```bash
python -m venv .venv
# macOS/Linux: source .venv/bin/activate
# Windows: .venv\\Scripts\\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```

## Repository structure
- `app.py` — Streamlit dashboard
- `data/limma_tumor_vs_normal_all_genes.csv` — precomputed GSE183947 DEG result table
- `data/sample_metadata.csv` — sample annotations for the bundled starting cohort
- `data/analysis_provenance.txt` — processing notes
- `scripts/run_gse42568.R` — GSE42568 workflow starter for platform-appropriate limma analysis
- `templates/deg_upload_template.csv` — example input format
- `.github/workflows/ci.yml` — syntax and required-file smoke checks

## Data format for added cohorts
Upload one CSV per cohort, one row per gene, with columns equivalent to `gene_id`, `logFC`, `P.Value`, and `adj.P.Val`. Common aliases are supported. Use a consistent gene identifier system and confirm contrast coding before comparing results.

## Datasets and methods
- **GSE183947**: supplied processed-expression reanalysis, limma on `log2(FPKM + 1)`. FPKM values are not raw counts; DESeq2 was not used for this dataset.
- **GSE42568**: Affymetrix microarray cohort; analyze separately with platform-appropriate methods, and harmonize identifiers before cross-cohort comparison.
- **Potential future cohort**: TCGA-BRCA with a documented tumor-versus-normal contrast and data version.

The app displays precomputed results and does not download GEO/TCGA data, process FASTQ files, or fit per-cohort differential-expression models at runtime. The optional pathway module expects a user-provided GMT gene-set file. The target module only matches user-supplied evidence and does not invent target/drug associations.

## Interpretation
Cohorts are analyzed separately. Missing gene mappings remain missing rather than being treated as zero effects or p-values of 1. Fisher's method combines raw p-values, with BH correction across genes; it is not a formal effect-size meta-analysis. Direction consistency, candidate-prioritization scores, and immune-signature summaries are exploratory and do not validate biomarkers, establish causality, demonstrate pathway activation, or estimate cell abundance. Confirm sample metadata, pairing, contrast direction, normalization, gene identifiers, and source provenance before interpreting results.

## Deploy
### Streamlit Community Cloud
Create a new app from this repository, branch `main`, entrypoint `app.py`.

### Render
Build: `pip install -r requirements.txt`

Start: `streamlit run app.py --server.port $PORT --server.address 0.0.0.0`

## Checks
GitHub Actions checks Python syntax and required project files. This does not replace interactive testing or a successful hosted deployment.

## Disclaimer
Research and educational use only. Not for diagnosis, prognosis, or treatment decisions.
