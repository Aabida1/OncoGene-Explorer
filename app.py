from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

ROOT = Path(__file__).parent
DATA_FILE = ROOT / "data" / "limma_tumor_vs_normal_all_genes.csv"
META_FILE = ROOT / "data" / "sample_metadata.csv"
PROVENANCE_FILE = ROOT / "data" / "analysis_provenance.txt"

st.set_page_config(page_title="OncoGene Explorer", page_icon="🧬", layout="wide")
st.title("OncoGene Explorer")
st.caption("Cancer transcriptomics · exploratory research only")

ALIASES = {
    "gene_id": ["gene_id", "gene", "gene_symbol", "symbol", "genes", "gene_name"],
    "logFC": ["logfc", "log2foldchange", "log2fc", "log_fold_change", "logfoldchange"],
    "P.Value": ["p.value", "pvalue", "p_value", "pval"],
    "adj.P.Val": ["adj.p.val", "padj", "fdr", "qvalue", "adjusted_p_value", "adj_p_value"],
}

def normalize(df):
    df = df.copy()
    lookup = {str(c).strip().lower(): c for c in df.columns}
    rename = {}
    for target, aliases in ALIASES.items():
        for alias in aliases:
            if alias in lookup:
                rename[lookup[alias]] = target
                break
    df = df.rename(columns=rename)
    missing = [x for x in ALIASES if x not in df.columns]
    if missing:
        raise ValueError("Required columns missing: " + ", ".join(missing))
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    for c in ["logFC", "P.Value", "adj.P.Val"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=list(ALIASES))
    df = df[(df.gene_id != "") & df["P.Value"].between(0, 1) & df["adj.P.Val"].between(0, 1)]
    return df.drop_duplicates("gene_id", keep="first")

@st.cache_data(show_spinner=False)
def load_default(path):
    return normalize(pd.read_csv(path))

def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out

def compare_cohorts(cohorts):
    if len(cohorts) < 2:
        return pd.DataFrame()
    merged = None
    pcols, fccols = [], []
    for i, (_, df) in enumerate(cohorts.items()):
        part = df[["gene_id", "logFC", "P.Value"]].rename(columns={"logFC": f"logFC_{i}", "P.Value": f"p_{i}"})
        pcols.append(f"p_{i}")
        fccols.append(f"logFC_{i}")
        merged = part if merged is None else merged.merge(part, on="gene_id", how="outer")
    pmat = merged[pcols]
    k = pmat.notna().sum(axis=1)
    stat = -2 * np.log(pmat.clip(lower=1e-300, upper=1).fillna(1)).sum(axis=1)
    from scipy.stats import chi2
    p = chi2.sf(stat, 2 * k)
    p[k < 2] = np.nan
    valid = np.isfinite(p)
    fdr = np.full(len(merged), np.nan)
    fdr[valid] = bh_adjust(p[valid])
    signs = np.sign(merged[fccols].to_numpy(dtype=float))
    consensus = []
    for row in signs:
        v = row[np.isfinite(row)]
        consensus.append(bool(len(v) >= 2 and np.all(v == v[0])))
    merged["cohorts_present"] = k
    merged["cohorts_available"] = len(cohorts)
    merged["fisher_p"] = p
    merged["fisher_fdr"] = fdr
    merged["direction_consistent"] = consensus
    merged["median_logFC_descriptive"] = merged[fccols].median(axis=1, skipna=True)
    return merged.sort_values("fisher_fdr", na_position="last")

with st.sidebar:
    st.header("Analysis controls")
    uploads = st.file_uploader("Upload cohort DEG CSV files", type="csv", accept_multiple_files=True, help="Required fields: gene identifier, log fold change, raw p-value, adjusted p-value.")
    fdr_cutoff = st.slider("Adjusted p-value cutoff", 0.001, 0.20, 0.05, 0.001, format="%.3f")
    logfc_cutoff = st.slider("Minimum absolute log2 fold change", 0.0, 3.0, 1.0, 0.1)

cohorts = {}
if DATA_FILE.exists():
    try:
        cohorts["GSE183947"] = load_default(str(DATA_FILE))
    except Exception as exc:
        st.error(f"Could not load bundled GSE183947 data: {exc}")
else:
    st.error("Missing bundled result file: data/limma_tumor_vs_normal_all_genes.csv")
    st.stop()

for upload in uploads or []:
    label = Path(upload.name).stem
    if label.endswith("_all_genes"):
        label = label[:-10]
    if label in cohorts:
        label += " (uploaded)"
    try:
        cohorts[label] = normalize(pd.read_csv(upload))
    except Exception as exc:
        st.sidebar.error(f"{upload.name}: {exc}")

with st.sidebar:
    active_name = st.selectbox("Active cohort (controls volcano and DEG table)", list(cohorts), key="active_cohort")
    compare_names = st.multiselect("Cohorts for comparison", list(cohorts), default=list(cohorts))
active = cohorts[active_name].copy()
active["minus_log10_p"] = -np.log10(active["P.Value"].clip(lower=1e-300))
active["status"] = np.select(
    [(active["adj.P.Val"] < fdr_cutoff) & (active.logFC >= logfc_cutoff),
     (active["adj.P.Val"] < fdr_cutoff) & (active.logFC <= -logfc_cutoff)],
    ["Higher in contrast numerator", "Lower in contrast numerator"], default="Not passing both cutoffs")

st.markdown("Explore precomputed differential-expression results, compare independent cohorts, and investigate pathways using your own gene-set libraries.")
st.warning("Research use only. Confirm sample mapping, contrast direction, normalization and identifier mapping before biological interpretation.")
m = st.columns(4)
m[0].metric("Active cohort", active_name)
m[1].metric("Genes available", f"{len(active):,}")
m[2].metric("FDR-significant", f"{int((active['adj.P.Val'] < fdr_cutoff).sum()):,}")
m[3].metric("Pass FDR + logFC", f"{int(((active['adj.P.Val'] < fdr_cutoff) & (active.logFC.abs() >= logfc_cutoff)).sum()):,}")

tabs = st.tabs(["Overview", "DEG Results", "Cohort Comparison", "Pathway Biology", "Immune Signatures", "Target Evidence", "Samples & Methods"])
with tabs[0]:
    st.subheader(f"Volcano plot — {active_name}")
    fig = px.scatter(active, x="logFC", y="minus_log10_p", color="status", hover_name="gene_id", hover_data={"P.Value": ":.3g", "adj.P.Val": ":.3g", "logFC": ":.3f"}, opacity=0.75)
    fig.add_hline(y=-np.log10(fdr_cutoff), line_dash="dash", annotation_text="Raw-p guide only; use adjusted p-values for significance")
    fig.add_vline(x=logfc_cutoff, line_dash="dot")
    fig.add_vline(x=-logfc_cutoff, line_dash="dot")
    fig.update_layout(height=560, xaxis_title="log2 fold change", yaxis_title="−log10(raw p-value)")
    st.plotly_chart(fig, use_container_width=True, key=f"volcano_{active_name}_{fdr_cutoff}_{logfc_cutoff}")
    st.caption("The vertical guides show the selected fold-change threshold. The horizontal guide is only a visual reference because the y-axis uses raw p-values while filtering uses adjusted p-values.")
with tabs[1]:
    view = st.radio("Rows", ["Pass both cutoffs", "FDR-significant", "All genes"], horizontal=True)
    table = active.copy()
    if view == "Pass both cutoffs":
        table = table[(table["adj.P.Val"] < fdr_cutoff) & (table.logFC.abs() >= logfc_cutoff)]
    elif view == "FDR-significant":
        table = table[table["adj.P.Val"] < fdr_cutoff]
    query = st.text_input("Search gene identifier")
    if query:
        table = table[table.gene_id.str.contains(query, case=False, regex=False, na=False)]
    sort_by = st.selectbox("Sort by", ["adj.P.Val", "P.Value", "logFC", "gene_id"])
    table = table.sort_values(sort_by, ascending=sort_by in ["adj.P.Val", "P.Value", "gene_id"])
    st.write(f"{len(table):,} rows")
    st.dataframe(table, use_container_width=True, hide_index=True, height=500)
    st.download_button("Download displayed results", table.to_csv(index=False).encode(), "oncogene_DEG_results.csv", "text/csv")
    st.download_button("Download all genes in active cohort", cohorts[active_name].to_csv(index=False).encode(), "oncogene_all_genes.csv", "text/csv")
with tabs[2]:
    selected = {n: cohorts[n] for n in compare_names}
    if len(selected) < 2:
        st.info("Upload or select at least two cohorts to compare them.")
    else:
        result = compare_cohorts(selected)
        st.caption("Cohorts are matched by exact gene identifier. Missing values remain missing. Fisher's method combines raw p-values; it is not an effect-size meta-analysis.")
        st.dataframe(result[result.cohorts_present >= 2], use_container_width=True, hide_index=True, height=480)
        st.download_button("Download comparison", result.to_csv(index=False).encode(), "oncogene_cohort_comparison.csv", "text/csv")
with tabs[3]:
    st.subheader("Pathway over-representation analysis")
    gmt = st.file_uploader("Upload GMT gene-set file", type="gmt")
    gene_pick = st.radio("Input gene list", ["Pass both cutoffs", "FDR-significant"], horizontal=True)
    minimum = st.number_input("Minimum set size after intersection", min_value=3, max_value=500, value=5)
    if gmt:
        picked = active[active["adj.P.Val"] < fdr_cutoff]
        if gene_pick == "Pass both cutoffs":
            picked = picked[picked.logFC.abs() >= logfc_cutoff]
        selected_genes, universe = set(picked.gene_id), set(active.gene_id)
        records = []
        from scipy.stats import hypergeom
        for line in gmt.getvalue().decode("utf-8").splitlines():
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            gs = set(parts[2:]) & universe
            overlap = selected_genes & gs
            if len(gs) >= int(minimum) and overlap:
                records.append({"gene_set": parts[0], "description": parts[1], "set_size_in_universe": len(gs), "overlap_n": len(overlap), "overlap_genes": ", ".join(sorted(overlap)), "p_value": float(hypergeom.sf(len(overlap)-1, len(universe), len(gs), len(selected_genes)))})
        enrich = pd.DataFrame(records)
        if not enrich.empty:
            enrich["FDR"] = bh_adjust(enrich.p_value.to_numpy())
            enrich = enrich.sort_values("FDR")
            st.dataframe(enrich, use_container_width=True, hide_index=True)
            st.download_button("Download enrichment results", enrich.to_csv(index=False).encode(), "oncogene_pathway_enrichment.csv", "text/csv")
        else:
            st.info("No sets overlapped the selected genes with the chosen size cutoff.")
    else:
        st.info("Upload a GMT file from a resource whose licensing terms permit your intended use. No gene-set library is bundled.")
with tabs[4]:
    st.warning("Gene-list overlap is not sample-level immune-cell deconvolution or a measure of cell abundance.")
    signatures = {"Cytotoxic lymphocyte-associated": "CD8A;CD8B;GZMB;PRF1;NKG7;GNLY", "Interferon-response-associated": "STAT1;ISG15;IFIT1;IFIT3;MX1;CXCL10", "Antigen-presentation-associated": "B2M;HLA-A;HLA-B;HLA-C;TAP1;TAP2", "T-cell-associated": "CD3D;CD3E;TRAC;CD2;IL7R", "Macrophage-associated": "CD68;CSF1R;C1QA;C1QB;TYROBP", "Checkpoint-associated": "CD274;PDCD1;CTLA4;LAG3;TIGIT;HAVCR2"}
    sig_upload = st.file_uploader("Optional custom signature CSV (signature, genes)", type="csv")
    if sig_upload:
        try:
            sdf = pd.read_csv(sig_upload)
            if not {"signature", "genes"}.issubset(sdf.columns):
                st.error("Expected columns: signature, genes. Separate genes with semicolons or commas.")
            else:
                for _, row in sdf.iterrows():
                    signatures[str(row.signature)] = str(row.genes).replace(",", ";")
        except Exception as exc:
            st.error(f"Could not read custom signatures: {exc}")
    ids = set(active.gene_id)
    hits = set(active[(active["adj.P.Val"] < fdr_cutoff) & (active.logFC.abs() >= logfc_cutoff)].gene_id)
    summaries = []
    for name, text in signatures.items():
        genes = {g.strip() for g in text.split(";") if g.strip()} & ids
        overlap = genes & hits
        summaries.append({"signature": name, "genes_mapped": len(genes), "genes_passing_cutoffs": len(overlap), "overlap_genes": ", ".join(sorted(overlap))})
    st.dataframe(pd.DataFrame(summaries), use_container_width=True, hide_index=True)
    st.caption("Exact identifier text matching is used; harmonize gene symbols and Ensembl IDs before interpretation.")
with tabs[5]:
    st.write("Upload sourced target evidence to intersect it with candidate genes. This module does not fetch or invent target/drug evidence.")
    st.caption("Required columns: gene_id, evidence_type, source_name. Suggested: drug_or_compound, cancer_context, evidence_summary, source_url, evidence_level.")
    evidence_file = st.file_uploader("Evidence CSV", type="csv")
    if evidence_file:
        try:
            evidence = pd.read_csv(evidence_file)
            required = {"gene_id", "evidence_type", "source_name"}
            if not required.issubset(evidence.columns):
                st.error("Missing columns: " + ", ".join(sorted(required - set(evidence.columns))))
            else:
                candidate = set(active[(active["adj.P.Val"] < fdr_cutoff) & (active.logFC.abs() >= logfc_cutoff)].gene_id)
                evidence["candidate_overlap"] = evidence.gene_id.astype(str).str.strip().isin(candidate)
                st.dataframe(evidence, use_container_width=True, hide_index=True)
                st.download_button("Download evidence overlap", evidence.to_csv(index=False).encode(), "oncogene_target_evidence.csv", "text/csv")
        except Exception as exc:
            st.error(f"Could not read evidence: {exc}")
with tabs[6]:
    st.write(f"Active cohort: **{active_name}** · {len(active):,} genes")
    st.write("The dashboard displays supplied/precomputed differential-expression results; it does not fit a model or process FASTQ data at runtime.")
    if META_FILE.exists():
        try:
            meta = pd.read_csv(META_FILE)
            st.write(f"Bundled metadata rows: {len(meta):,}")
            st.dataframe(meta, use_container_width=True, hide_index=True)
        except Exception as exc:
            st.warning(f"Could not read metadata: {exc}")
    if PROVENANCE_FILE.exists():
        st.subheader("Provenance note")
        st.code(PROVENANCE_FILE.read_text(errors="replace"))
    st.markdown("**Important limits**\n- Confirm contrast direction and any sample pairing from source metadata.\n- Cross-cohort matching uses exact identifiers and needs harmonization.\n- The consensus output is exploratory p-value aggregation, not formal meta-analysis.\n- Immune outputs are gene-list overlap, not cell abundance estimates.\n- Not intended for clinical decision-making.")