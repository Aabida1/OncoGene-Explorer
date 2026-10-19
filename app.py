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
st.markdown("""
<style>
:root { --og-ink: #17263c; --og-muted: #64748b; --og-teal: #0f766e; --og-line: #dce6ed; }
.block-container {padding-top: 1.25rem; padding-bottom: 3rem; max-width: 1600px;}
[data-testid="stAppViewContainer"] {background: linear-gradient(180deg, rgba(13,148,136,.035) 0%, transparent 280px);}
[data-testid="stSidebar"] {border-right: 1px solid var(--og-line);}
[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] h2 {letter-spacing: -.02em;}
[data-testid="stMetric"] {background: var(--secondary-background-color); border: 1px solid var(--og-line); padding: 16px 18px; border-radius: 15px; box-shadow: 0 2px 8px rgba(15,23,42,.025);}
[data-testid="stMetricLabel"] {font-size: .83rem; color: var(--og-muted);}
[data-testid="stMetricValue"] {font-weight: 750; letter-spacing: -.035em;}
div[data-testid="stTabs"] {margin-top: 1.1rem;}
div[data-testid="stTabs"] button {font-weight: 650; padding: .7rem .8rem; border-radius: 9px 9px 0 0;}
[data-testid="stVerticalBlock"] > [data-testid="stHorizontalBlock"] {gap: 1rem;}
.og-hero {padding: 1.5rem 1.65rem; border-radius: 20px; border: 1px solid rgba(15,118,110,.22); background: linear-gradient(115deg, rgba(15,118,110,.12), rgba(59,130,246,.06) 65%, rgba(255,255,255,.2)); margin: .3rem 0 1.1rem 0;}
.og-eyebrow {text-transform: uppercase; letter-spacing: .13em; font-size: .72rem; font-weight: 800; color: var(--og-teal); margin-bottom: .55rem;}
.og-title {font-size: clamp(1.8rem, 3vw, 2.65rem); line-height: 1.1; letter-spacing: -.045em; font-weight: 800; color: var(--og-ink); margin: 0 0 .65rem 0;}
.og-subtitle {font-size: 1rem; line-height: 1.65; max-width: 900px; color: var(--og-muted); margin: 0;}
.og-chip {display:inline-block; border:1px solid var(--og-line); border-radius:999px; padding:.3rem .65rem; font-size:.76rem; font-weight:650; margin:.8rem .35rem 0 0; background:rgba(255,255,255,.55); color:var(--og-ink);}
[data-testid="stDataFrame"] {border: 1px solid var(--og-line); border-radius: 12px; overflow: hidden;}
[data-testid="stAlert"] {border-radius: 12px;}
hr {border-color: var(--og-line);}
@media (max-width: 800px) { .og-hero {padding: 1.15rem; border-radius: 15px;} }
</style>
""", unsafe_allow_html=True)

@st.cache_data
def load_results(path):
    df = pd.read_csv(path)
    required = {"gene_id", "logFC", "P.Value", "adj.P.Val"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Required columns missing: {', '.join(sorted(missing))}")
    for col in ["logFC", "P.Value", "adj.P.Val"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df = df.dropna(subset=["gene_id", "logFC", "P.Value", "adj.P.Val"])
    df = df[(df["P.Value"] >= 0) & (df["adj.P.Val"] >= 0)]
    df["neglog10_fdr"] = -np.log10(df["adj.P.Val"].clip(lower=1e-300))
    return df

def normalize_result_columns(df, cohort_name):
    """Accept common DEG result headers and standardize them for cross-cohort comparison."""
    df = df.copy()
    aliases = {
        "gene_id": ["gene_id", "gene", "gene_symbol", "symbol", "Gene", "Gene.symbol", "genes"],
        "logFC": ["logFC", "log2FoldChange", "log2FC", "log2_fold_change", "log2.fold.change"],
        "P.Value": ["P.Value", "pvalue", "p_value", "pval", "PValue", "p.value"],
        "adj.P.Val": ["adj.P.Val", "padj", "FDR", "fdr", "adjusted_p_value", "adj_p_value", "qvalue", "q_value"]
    }
    rename = {}
    for target, options in aliases.items():
        if target not in df.columns:
            found = next((c for c in options if c in df.columns), None)
            if found:
                rename[found] = target
    df = df.rename(columns=rename)
    required = {"gene_id", "logFC", "P.Value", "adj.P.Val"}
    if not required.issubset(df.columns):
        missing = ", ".join(sorted(required - set(df.columns)))
        raise ValueError(f"{cohort_name}: missing required columns ({missing}). See the upload template in the README.")
    for c in ["logFC", "P.Value", "adj.P.Val"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["gene_id"] = df["gene_id"].astype(str).str.strip()
    df = df.dropna(subset=["gene_id", "logFC", "P.Value", "adj.P.Val"])
    df = df[(df["gene_id"] != "") & (df["P.Value"].between(0, 1)) & (df["adj.P.Val"].between(0, 1))]
    # One row per gene per cohort; keep the row with the smallest adjusted p-value.
    df = df.sort_values("adj.P.Val").drop_duplicates("gene_id")
    df["cohort"] = cohort_name
    df["neglog10_fdr"] = -np.log10(df["adj.P.Val"].clip(lower=1e-300))
    return df


def combine_cohorts(cohorts):
    """Combine p-values with Fisher's method and report direction consistency, not a pooled effect."""
    from scipy.stats import chi2
    rows = []
    names = list(cohorts)
    all_genes = sorted(set().union(*(set(df["gene_id"]) for df in cohorts.values())))
    for gene in all_genes:
        present = []
        for name, df in cohorts.items():
            hit = df[df["gene_id"] == gene]
            if not hit.empty:
                r = hit.iloc[0]
                present.append((name, float(r["logFC"]), max(float(r["P.Value"]), 1e-300), float(r["adj.P.Val"])))
        if not present:
            continue
        stat = -2 * sum(np.log(x[2]) for x in present)
        fisher_p = float(chi2.sf(stat, 2 * len(present)))
        signs = [np.sign(x[1]) for x in present if x[1] != 0]
        same_direction = bool(signs) and (all(v > 0 for v in signs) or all(v < 0 for v in signs))
        rows.append({
            "gene_id": gene, "cohorts_present": len(present),
            "cohorts_tested": len(names), "coverage": len(present) / len(names),
            "direction_consistent": same_direction,
            "mean_logFC_descriptive": float(np.mean([x[1] for x in present])),
            "median_logFC_descriptive": float(np.median([x[1] for x in present])),
            "Fisher_combined_P": fisher_p,
            "cohorts": "; ".join(x[0] for x in present),
            "per_cohort_logFC": "; ".join(f"{x[0]}: {x[1]:.3f}" for x in present),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        # Benjamini-Hochberg correction across the combined gene-level Fisher p-values.
        pvals = out["Fisher_combined_P"].to_numpy()
        order = np.argsort(pvals)
        ranked = pvals[order] * len(pvals) / (np.arange(len(pvals)) + 1)
        adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
        fdr = np.empty_like(adjusted)
        fdr[order] = np.minimum(adjusted, 1.0)
        out["Fisher_combined_FDR"] = fdr
        out = out.sort_values(["Fisher_combined_FDR", "direction_consistent"], ascending=[True, False])
    return out


@st.cache_data
def load_metadata(path):
    df = pd.read_csv(path)
    expected = {"geo_accession", "title", "condition", "pair"}
    if expected.issubset(df.columns):
        return df
    return pd.DataFrame()

st.markdown("""
<div class="og-hero">
  <div class="og-eyebrow">Cancer bioinformatics · Research workspace</div>
  <div class="og-title">OncoGene Explorer</div>
  <p class="og-subtitle">From differential expression to reproducible cancer biology. Compare independent cohorts, map pathway-level signals, and assemble transparent evidence for candidate biomarkers.</p>
  <span class="og-chip">Breast cancer · Initial focus</span><span class="og-chip">Cross-cohort analysis</span><span class="og-chip">Evidence-led interpretation</span>
</div>
""", unsafe_allow_html=True)
with st.expander("About this workspace · data and interpretation", expanded=False):
    st.markdown("""
    - **Bundled starting cohort:** GEO GSE183947.
    - **Cross-cohort inputs:** upload precomputed gene-level differential-expression tables in the sidebar.
    - **Interpretation:** results are exploratory research evidence, not clinical diagnoses or validated treatment recommendations.
    - **Important:** cohorts are analyzed separately; uploaded gene identifiers and contrast directions must be harmonized before comparing studies.
    """)

if not DATA_FILE.exists():
    st.error("The bundled results CSV is missing. Ensure data/limma_tumor_vs_normal_all_genes.csv exists.")
    st.stop()
try:
    primary = normalize_result_columns(pd.read_csv(DATA_FILE), "GSE183947")
except Exception as exc:
    st.error(f"Could not load the bundled analysis file: {exc}")
    st.stop()

with st.sidebar:
    st.markdown("### Study setup")
    st.caption("Configure the evidence included in this session.")
    st.caption("Start with GSE183947. Add independent cohorts by uploading their precomputed DEG tables.")
    uploaded_files = st.file_uploader("Add cohort DEG CSV files", type=["csv"], accept_multiple_files=True,
        help="Each CSV must contain gene_id/gene symbol, logFC, P.Value and adj.P.Val (or common equivalents). Upload one file per cohort.")
    cohorts = {"GSE183947": primary}
    for i, file in enumerate(uploaded_files or []):
        label = Path(file.name).stem.strip() or f"Uploaded cohort {i+1}"
        if label in cohorts:
            label = f"{label}_{i+1}"
        try:
            cohorts[label] = normalize_result_columns(pd.read_csv(file), label)
        except Exception as exc:
            st.error(str(exc))
    st.divider()
    selected_names = st.multiselect(
        "Cohorts included in cross-study comparison",
        options=list(cohorts),
        default=list(cohorts),
        key="included_cohorts",
        help="Controls which cohorts contribute to cross-cohort comparison and consensus. This does not by itself choose the volcano plot cohort."
    )
    if not selected_names:
        st.warning("Select at least one cohort.")
        st.stop()

    # Keep the active-cohort widget synchronized when the included-cohort list changes.
    # Without this guard, Streamlit can retain a previous widget value that is no longer
    # present in the new options list, which makes the plot selection confusing/stale.
    default_active = "GSE183947" if "GSE183947" in selected_names else selected_names[0]
    if st.session_state.get("active_cohort_selector") not in selected_names:
        st.session_state["active_cohort_selector"] = default_active
    active_cohort_name = st.selectbox(
        "Dataset shown in volcano plot and DEG table",
        options=selected_names,
        index=selected_names.index(default_active),
        key="active_cohort_selector",
        help="Choose the single cohort shown in Overview and Results. Changing the comparison cohort list above does not automatically change this selection."
    )
    fdr = st.selectbox("Adjusted p-value (FDR) cutoff", [0.01, 0.05, 0.10], index=1)
    fc = st.select_slider("Minimum absolute log₂ fold change", options=[0.0, 0.5, 1.0, 1.5, 2.0], value=1.0)
    direction = st.selectbox("Direction", ["Both directions", "Upregulated", "Downregulated"])
    st.divider()
    st.caption("Filters are applied to precomputed results. The app does not refit the individual cohort models.")

cohorts = {name: cohorts[name] for name in selected_names}
combined = combine_cohorts(cohorts)
active = cohorts[active_cohort_name]
sig = active["adj.P.Val"] < fdr
filtered = active[sig & (active["logFC"].abs() >= fc)].copy()
filtered["regulation"] = np.where(filtered["logFC"] > 0, "Upregulated", "Downregulated")
if direction == "Upregulated":
    filtered = filtered[filtered["logFC"] > 0]
elif direction == "Downregulated":
    filtered = filtered[filtered["logFC"] < 0]
up = int(((active["adj.P.Val"] < fdr) & (active["logFC"] >= fc)).sum())
down = int(((active["adj.P.Val"] < fdr) & (active["logFC"] <= -fc)).sum())
col1, col2, col3, col4 = st.columns(4)
col1.metric("Selected cohorts", f"{len(cohorts)}")
col2.metric("Genes available", f"{len(active):,}")
col3.metric("FDR-significant genes", f"{int(sig.sum()):,}")
col4.metric("DEGs at current cutoffs", f"{int((sig & (active['logFC'].abs() >= fc)).sum()):,}")
if len(cohorts) > 1:
    consensus_n = int((combined["direction_consistent"] & (combined["Fisher_combined_FDR"] < fdr) & (combined["coverage"] == 1.0)).sum())
    st.caption(f"Cross-cohort candidates with consistent direction, full coverage and combined FDR < {fdr}: {consensus_n:,}. This is an exploratory evidence count, not a validated biomarker count.")
st.caption(f"Upregulated*: {up:,} · Downregulated*: {down:,} in {active['cohort'].iloc[0]} at the selected FDR and fold-change cutoffs. A positive logFC depends on the source study's contrast direction.")

overview, cohort_tab, consensus_tab, pathway_tab, priority_tab, target_tab, immune_tab, report_tab, gene_tab, results_tab, samples_tab, methods_tab = st.tabs(["Overview", "Cohorts", "Consensus", "Pathways", "Prioritization", "Targets", "Immune context", "Evidence report", "Gene explorer", "Results", "Samples", "Methods"])
with overview:
    st.subheader("Expression landscape")
    st.caption("Explore effect size against statistical evidence. Hover over points to inspect gene-level results.")
    plot = active.copy()
    plot["category"] = np.select(
        [(plot["adj.P.Val"] < fdr) & (plot["logFC"] >= fc),
         (plot["adj.P.Val"] < fdr) & (plot["logFC"] <= -fc)],
        ["Upregulated", "Downregulated"], default="Not significant")
    fig = px.scatter(plot, x="logFC", y="neglog10_fdr", color="category", hover_name="gene_id",
        hover_data={"P.Value": ":.3e", "adj.P.Val": ":.3e", "logFC": ":.3f", "neglog10_fdr": False},
        color_discrete_map={"Upregulated": "#D65A55", "Downregulated": "#4386C6", "Not significant": "#B8C0CC"},
        opacity=0.75, render_mode="webgl",
        labels={"logFC": "log₂ fold change (tumor vs normal, source contrast)", "neglog10_fdr": "−log₁₀(adjusted p-value)", "category": "Category"})
    fig.add_hline(y=-np.log10(fdr), line_dash="dash", line_color="#6B7280", annotation_text=f"FDR = {fdr}")
    fig.add_vline(x=fc, line_dash="dot", line_color="#6B7280")
    if fc > 0: fig.add_vline(x=-fc, line_dash="dot", line_color="#6B7280")
    fig.update_layout(height=540, margin=dict(l=15, r=15, t=35, b=15), legend_title_text="Gene category")
    st.plotly_chart(fig, use_container_width=True, key=f"volcano_plot_{active_cohort_name}_{fdr}_{fc}_{direction}")
    st.subheader("Loaded cohorts")
    cohort_summary = pd.DataFrame([{
        "Cohort": n,
        "Genes available": len(d),
        f"FDR < {fdr}": int((d["adj.P.Val"] < fdr).sum()),
        f"DEGs: FDR < {fdr} and |logFC| ≥ {fc:g}": int(((d["adj.P.Val"] < fdr) & (d["logFC"].abs() >= fc)).sum()),
        "Unique gene IDs": int(d["gene_id"].nunique()),
        "Median logFC (descriptive)": round(float(d["logFC"].median()), 3)
    } for n, d in cohorts.items()])
    st.dataframe(cohort_summary, use_container_width=True, hide_index=True)
    st.caption("Gene counts refer to rows remaining after validation and duplicate gene IDs are collapsed. Cross-cohort overlap can be smaller than either cohort's gene count.")
    st.info(f"Volcano plot dataset: **{active_cohort_name}**. Use the sidebar control **Dataset shown in volcano plot and DEG table** to switch datasets. The comparison cohort selector above it controls consensus inputs separately.")
    if len(cohorts) == 1:
        st.info("Only GSE183947 is bundled. Upload a separate cohort's gene-level DEG CSV in the sidebar to enable cross-cohort comparison. GSE42568 is a microarray cohort, so use platform-aware analysis and gene-ID harmonization before comparing it with RNA-seq cohorts.")
    else:
        st.success("Cross-cohort summary is active. Candidate ranking combines p-values and reports direction consistency; it does not pool fold changes into a formal random-effects meta-analysis.")

with cohort_tab:
    st.subheader("Per-cohort effect-size comparison")
    st.write("Genes are matched by the exact `gene_id` text. Standardize identifiers (preferably HGNC symbols or Ensembl IDs mapped consistently) before uploading.")
    gene_pool = sorted(set().union(*(set(d["gene_id"]) for d in cohorts.values())))
    chosen_genes = st.multiselect("Genes to compare", gene_pool, default=[] , max_selections=30, placeholder="Choose genes, e.g. BRCA1 or MKI67")
    if chosen_genes:
        rows=[]
        for gene in chosen_genes:
            for name, d in cohorts.items():
                h=d[d.gene_id==gene]
                if not h.empty:
                    r=h.iloc[0]
                    rows.append({"Gene":gene,"Cohort":name,"Status":"Present","logFC":r.logFC,"P.Value":r["P.Value"],"FDR":r["adj.P.Val"]})
                else:
                    rows.append({"Gene":gene,"Cohort":name,"Status":"Not found / not mapped","logFC":np.nan,"P.Value":np.nan,"FDR":np.nan})
        comp=pd.DataFrame(rows)
        if not comp.empty:
            fig2=px.bar(comp,x="Gene",y="logFC",color="Cohort",barmode="group",hover_data=["P.Value","FDR"],title="Effect direction and magnitude by cohort")
            fig2.add_hline(y=0,line_dash="dash",line_color="gray")
            st.plotly_chart(fig2,use_container_width=True)
            st.dataframe(comp,use_container_width=True,hide_index=True)
    else:
        st.info("Select genes above to compare their fold changes across cohorts.")

with consensus_tab:
    st.subheader("Consensus engine: cross-cohort candidate table")
    if len(cohorts) < 2:
        st.info("Upload at least one independent cohort DEG CSV to calculate cross-cohort evidence.")
    else:
        min_coverage = st.slider("Minimum fraction of selected cohorts containing the gene", 0.5, 1.0, 1.0, 0.1)
        consistent_only = st.checkbox("Require consistent fold-change direction", value=True)
        candidates=combined[combined.coverage >= min_coverage].copy()
        if consistent_only: candidates=candidates[candidates.direction_consistent]
        candidates=candidates[candidates.Fisher_combined_FDR < fdr]
        st.caption("Fisher combined p-values test aggregate evidence; BH FDR is applied across genes. Direction consistency is reported separately. These are exploratory cross-study candidates, not validated biomarkers. P-values should be independent or appropriately interpreted; this is not a formal effect-size meta-analysis.")
        st.download_button("Download consensus evidence report (CSV)", candidates.to_csv(index=False).encode("utf-8"), "oncogene_consensus_evidence.csv", "text/csv")
        st.write(f"Candidates passing the current settings: {len(candidates):,}")
        st.dataframe(candidates,use_container_width=True,hide_index=True)
        st.download_button("Download cross-cohort candidates (CSV)",candidates.to_csv(index=False).encode("utf-8"),"oncogene_cross_cohort_candidates.csv","text/csv")
        st.download_button("Download all combined evidence (CSV)",combined.to_csv(index=False).encode("utf-8"),"oncogene_cross_cohort_all_genes.csv","text/csv")

with pathway_tab:
    st.subheader("Pathway mapper: functional enrichment")
    st.write("Test whether the selected gene list is over-represented in biological gene sets such as GO Biological Process, KEGG, Reactome, or MSigDB Hallmark. Supply a GMT file to keep the analysis reproducible and avoid silently relying on a changing online database.")
    gmt_file = st.file_uploader("Upload a gene-set library (.gmt)", type=["gmt"], key="pathway_gmt")
    pathway_direction = st.radio("Gene set to analyze", ["Upregulated", "Downregulated", "All significant genes"], horizontal=True, key="pathway_direction")
    min_set_size = st.number_input("Minimum genes per pathway", min_value=3, max_value=100, value=5, step=1)
    max_set_size = st.number_input("Maximum genes per pathway", min_value=10, max_value=5000, value=500, step=10)
    if gmt_file is not None:
        gene_sets = {}
        for line in gmt_file.getvalue().decode("utf-8", errors="replace").splitlines():
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                gene_sets[parts[0]] = set(g.strip() for g in parts[2:] if g.strip())
        sig_df = active[active["adj.P.Val"] < fdr].copy()
        if pathway_direction == "Upregulated": sig_df = sig_df[sig_df["logFC"] >= fc]
        elif pathway_direction == "Downregulated": sig_df = sig_df[sig_df["logFC"] <= -fc]
        else: sig_df = sig_df[sig_df["logFC"].abs() >= fc]
        universe = set(active["gene_id"].astype(str))
        query_genes = set(sig_df["gene_id"].astype(str))
        from scipy.stats import hypergeom
        enrich_rows = []
        for pathway, members in gene_sets.items():
            members = members & universe
            if len(members) < int(min_set_size) or len(members) > int(max_set_size): continue
            overlap = members & query_genes
            if not overlap: continue
            pval = float(hypergeom.sf(len(overlap)-1, len(universe), len(members), len(query_genes)))
            enrich_rows.append({"pathway": pathway, "overlap_count": len(overlap), "gene_set_size_in_universe": len(members), "query_gene_count": len(query_genes), "overlap_genes": "; ".join(sorted(overlap)), "P.Value": pval, "gene_ratio": len(overlap)/max(len(query_genes),1), "background_ratio": len(members)/max(len(universe),1)})
        enr = pd.DataFrame(enrich_rows)
        if not enr.empty:
            pvals = enr["P.Value"].to_numpy(); order = np.argsort(pvals); ranked = pvals[order] * len(pvals) / (np.arange(len(pvals))+1); adjusted = np.minimum.accumulate(ranked[::-1])[::-1]; fdrs = np.empty_like(adjusted); fdrs[order] = np.minimum(adjusted,1.0); enr["FDR"] = fdrs
            enr = enr.sort_values("FDR")
            st.caption(f"Tested {len(enr):,} pathways with nonzero overlap; background universe: {len(universe):,} genes. Interpretation depends on consistent gene identifiers and a suitable background.")
            top_enr = enr.head(25).sort_values("FDR", ascending=True)
            figp = px.bar(top_enr, x="-log10(FDR)" if False else "FDR", y="pathway", orientation="h", hover_data=["overlap_count", "overlap_genes"], title="Top enriched gene sets")
            figp.update_xaxes(autorange="reversed", title="FDR (lower is stronger evidence)")
            figp.update_layout(height=max(420, 25*len(top_enr)), yaxis={"categoryorder":"total ascending"})
            st.plotly_chart(figp, use_container_width=True)
            st.dataframe(enr, use_container_width=True, hide_index=True)
            st.download_button("Download pathway enrichment (CSV)", enr.to_csv(index=False).encode("utf-8"), "oncogene_pathway_enrichment.csv", "text/csv")
        else:
            st.info("No pathways overlapped the selected gene list under the current size limits. Check identifier matching, thresholds, and the GMT file.")
    else:
        st.info("Upload a GMT gene-set file to run local over-representation analysis. Suitable libraries include GO, Reactome, KEGG, and MSigDB collections, subject to their terms of use.")

with priority_tab:
    st.subheader("Transparent biomarker-candidate prioritization")
    st.write("This is an exploratory ranking aid—not a clinical biomarker validation model. It rewards replicated direction, cohort coverage, combined statistical evidence, and effect magnitude while exposing the component metrics.")
    if len(cohorts) < 2:
        st.info("Upload at least one independent cohort to enable cross-cohort prioritization.")
    else:
        min_coverage_priority = st.slider("Minimum cohort coverage", 0.5, 1.0, 1.0, 0.1, key="priority_coverage")
        require_direction = st.checkbox("Require consistent direction across cohorts", value=True, key="priority_direction")
        rank = combined[combined["coverage"] >= min_coverage_priority].copy()
        if require_direction: rank = rank[rank["direction_consistent"]]
        if not rank.empty:
            # Transparent 0-100 heuristic: 40% statistical evidence, 30% cohort coverage, 20% direction consistency, 10% effect magnitude.
            evidence = (-np.log10(rank["Fisher_combined_FDR"].clip(lower=1e-300))).clip(0, 20) / 20
            magnitude = (rank["median_logFC_descriptive"].abs().clip(0, 2) / 2)
            rank["priority_score_0_100"] = (40*evidence + 30*rank["coverage"].clip(0,1) + 20*rank["direction_consistent"].astype(float) + 10*magnitude).round(2)
            rank["interpretation"] = np.select([rank["direction_consistent"] & (rank["coverage"] >= 0.99), rank["direction_consistent"], ~rank["direction_consistent"]], ["Replicated direction across all selected cohorts", "Consistent direction; incomplete cohort coverage", "Discordant direction; review cohort context"], default="Review evidence")
            rank = rank.sort_values(["priority_score_0_100", "Fisher_combined_FDR"], ascending=[False, True])
            st.warning("The score is a user-facing heuristic, not a validated predictive model. A high score does not establish diagnostic, prognostic, or therapeutic utility.")
            st.dataframe(rank[["gene_id", "priority_score_0_100", "interpretation", "cohorts_present", "cohorts_tested", "coverage", "direction_consistent", "median_logFC_descriptive", "Fisher_combined_P", "Fisher_combined_FDR", "per_cohort_logFC"]], use_container_width=True, hide_index=True)
            st.download_button("Download prioritized candidates (CSV)", rank.to_csv(index=False).encode("utf-8"), "oncogene_prioritized_candidates.csv", "text/csv")
        else:
            st.info("No candidates match these settings. Lower the coverage threshold or review the uploaded cohort tables.")


with target_tab:
    st.subheader("Therapeutic Target Evidence Bridge")
    st.write("Connect prioritized genes to external target and functional-dependency evidence. Upload a curated, source-traceable CSV; this module does not infer that a gene is druggable or that a drug will benefit patients.")
    template = pd.DataFrame(columns=["gene_id", "target_name", "drug_or_compound", "evidence_type", "cancer_context", "evidence_summary", "source_name", "source_url", "evidence_level"])
    st.download_button("Download target-evidence CSV template", template.to_csv(index=False).encode("utf-8"), "target_evidence_template.csv", "text/csv")
    st.markdown("**Useful external resources**")
    st.markdown("- [Open Targets Platform](https://platform.opentargets.org/) — target–disease evidence\n- [DepMap](https://depmap.org/portal/) — cancer-cell-line dependency evidence\n- [STRING](https://string-db.org/) — protein functional associations\n- [DGIdb](https://www.dgidb.org/) — gene–drug interaction evidence")
    target_file = st.file_uploader("Upload target evidence CSV", type=["csv"], key="target_evidence_csv")
    if target_file is not None:
        try:
            targets = pd.read_csv(target_file)
            required = {"gene_id", "evidence_type", "source_name"}
            missing = required - set(targets.columns)
            if missing:
                st.error("Missing required columns: " + ", ".join(sorted(missing)))
            else:
                targets["gene_id"] = targets["gene_id"].astype(str).str.strip()
                candidate_genes = set(combined.loc[combined.get("direction_consistent", pd.Series(dtype=bool)).astype(bool), "gene_id"].astype(str)) if not combined.empty else set()
                targets["matches_current_candidates"] = targets["gene_id"].isin(candidate_genes)
                if "source_url" in targets.columns:
                    targets["source_url"] = targets["source_url"].fillna("")
                st.metric("Evidence records", f"{len(targets):,}")
                st.metric("Records matching consistent candidate genes", f"{int(targets['matches_current_candidates'].sum()):,}")
                if "evidence_level" in targets.columns:
                    st.caption("Evidence levels are displayed as supplied; ensure they use a documented, consistent scheme.")
                st.dataframe(targets, use_container_width=True, hide_index=True)
                st.download_button("Download matched target evidence", targets.to_csv(index=False).encode("utf-8"), "oncogene_target_evidence_matched.csv", "text/csv")
                st.session_state["target_evidence_report"] = targets
        except Exception as exc:
            st.error(f"Could not read target evidence file: {exc}")
    else:
        st.info("No target evidence is bundled by default. Upload a sourced CSV to avoid presenting unsourced gene–drug claims as facts.")

with immune_tab:
    st.subheader("Tumor Microenvironment & Immune Signature Analysis")
    st.write("Score gene signatures against the uploaded differential-expression results using a transparent signed-rank approach. This is a gene-set-level proxy based on DE statistics, not a cell-infiltration estimate; sample-level expression data are needed for sample-wise scoring.")
    starter_signatures = {
        "T_cell_associated": ["CD3D", "CD3E", "CD3G", "TRAC", "CD247", "LCK", "TRBC1", "TRBC2"],
        "Cytotoxic_lymphocyte": ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "GZMH", "GNLY", "CTSW"],
        "Interferon_response": ["STAT1", "IRF1", "ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "CXCL10"],
        "Macrophage_associated": ["C1QA", "C1QB", "C1QC", "CSF1R", "CD68", "TYROBP", "LST1", "AIF1"],
        "Antigen_presentation": ["B2M", "TAP1", "TAP2", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1"],
        "Checkpoint_associated": ["CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1"]
    }
    sig_df = pd.DataFrame([{"signature": name, "genes": ";".join(genes), "source_note": "Starter gene set for exploratory screening; verify against a published signature before biological claims"} for name, genes in starter_signatures.items()])
    st.download_button("Download editable starter signatures CSV", sig_df.to_csv(index=False).encode("utf-8"), "oncogene_immune_signature_starters.csv", "text/csv")
    immune_file = st.file_uploader("Optional custom immune signature CSV (signature, genes)", type=["csv"], key="immune_signatures_csv")
    signatures = dict(starter_signatures)
    if immune_file is not None:
        try:
            custom = pd.read_csv(immune_file)
            if not {"signature", "genes"}.issubset(custom.columns):
                st.error("Custom signature CSV must have 'signature' and 'genes' columns; separate genes with semicolons or commas.")
            else:
                for _, row in custom.iterrows():
                    raw = str(row["genes"])
                    signatures[str(row["signature"])] = [g.strip() for g in raw.replace(",", ";").split(";") if g.strip()]
        except Exception as exc:
            st.error(f"Could not read signature file: {exc}")
    if st.button("Score immune signatures", key="score_immune_signatures"):
        expr = active.copy()
        expr["gene_key"] = expr["gene_id"].astype(str).str.upper().str.strip()
        expr = expr.drop_duplicates("gene_key").set_index("gene_key")
        out = []
        from scipy.stats import mannwhitneyu
        for name, genes in signatures.items():
            mapped = [g.upper() for g in genes if g.upper() in expr.index]
            if not mapped:
                out.append({"signature": name, "genes_in_signature": len(genes), "genes_matched": 0, "mean_logFC": np.nan, "median_logFC": np.nan, "upregulated_fraction": np.nan, "enrichment_p_value": np.nan, "matched_genes": ""})
                continue
            values = expr.loc[mapped, "logFC"].astype(float)
            # Directional enrichment: compare signed logFC of signature genes against all genes; exploratory only.
            try:
                pval = float(mannwhitneyu(values.to_numpy(), expr["logFC"].astype(float).to_numpy(), alternative="two-sided").pvalue)
            except Exception:
                pval = np.nan
            out.append({"signature": name, "genes_in_signature": len(genes), "genes_matched": len(mapped), "mean_logFC": float(values.mean()), "median_logFC": float(values.median()), "upregulated_fraction": float((values > 0).mean()), "enrichment_p_value": pval, "matched_genes": ";".join(mapped)})
        immune_results = pd.DataFrame(out)
        valid = immune_results["enrichment_p_value"].notna()
        immune_results["FDR"] = np.nan
        if valid.any():
            pv = immune_results.loc[valid, "enrichment_p_value"].to_numpy()
            order = np.argsort(pv); ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1); adj = np.minimum.accumulate(ranked[::-1])[::-1]
            corrected = np.empty_like(adj); corrected[order] = np.minimum(adj, 1.0)
            immune_results.loc[valid, "FDR"] = corrected
        st.session_state["immune_signature_results"] = immune_results
    if "immune_signature_results" in st.session_state:
        ir = st.session_state["immune_signature_results"]
        st.dataframe(ir, use_container_width=True, hide_index=True)
        st.download_button("Download immune signature results", ir.to_csv(index=False).encode("utf-8"), "oncogene_immune_signature_results.csv", "text/csv")
        st.warning("These starter signatures are illustrative and should be replaced or verified against published gene sets. Signature scoring from a DEG table is exploratory and does not estimate immune-cell abundance or establish immune activation.")
    else:
        st.info("Choose a cohort in the sidebar and click 'Score immune signatures'. Gene-symbol matching is required for these starter sets.")

with report_tab:
    st.subheader("Reproducible Evidence Report")
    st.write("Generate a downloadable Markdown report of the current cohort selection, thresholds, consensus evidence, and available pathway/immune/target results. Reports document analysis context; they do not constitute clinical validation.")
    report_lines = [
        "# OncoGene Explorer — Evidence Report", "", 
        "## Analysis configuration", f"- Cohorts selected: {', '.join(cohorts.keys())}",
        f"- Active cohort for single-cohort analyses: {active['cohort'].iloc[0]}", f"- FDR threshold: {fdr}", f"- Absolute log2 fold-change threshold: {fc}",
        f"- Genes in active cohort: {len(active)}", f"- FDR-significant genes in active cohort: {int((active['adj.P.Val'] < fdr).sum())}",
        "", "## Consensus analysis", "The application combines gene-level p-values with Fisher's method and applies Benjamini–Hochberg correction across genes. Direction consistency is evaluated separately. This is exploratory evidence aggregation, not a formal effect-size meta-analysis."
    ]
    if len(cohorts) > 1 and not combined.empty:
        top = combined.head(25)
        report_lines += ["", "### Top 25 combined-evidence candidates", "| Gene | Cohorts present | Direction consistent | Median logFC (descriptive) | Fisher FDR |", "|---|---:|---|---:|---:|"]
        for _, r in top.iterrows():
            report_lines.append(f"| {r['gene_id']} | {int(r['cohorts_present'])}/{int(r['cohorts_tested'])} | {bool(r['direction_consistent'])} | {r['median_logFC_descriptive']:.4g} | {r['Fisher_combined_FDR']:.4g} |")
    else:
        report_lines += ["", "Only one cohort is currently selected; cross-cohort consensus is not available."]
    report_lines += ["", "## Pathway and immune/target modules", "- Pathway enrichment requires a user-supplied GMT library; gene identifiers and universe affect results.", "- Immune signature scores are exploratory gene-set summaries from a DEG table, not cell-abundance estimates.", "- Therapeutic target evidence is only included when a user-supplied, source-traceable target-evidence CSV is uploaded.", "", "## Data provenance and limitations", "- The bundled GSE183947 results are processed-expression reanalysis results; consult the bundled provenance file and original GEO record.", "- Uploaded cohorts are assumed to have been analyzed with appropriate study/platform-aware methods.", "- Exact gene identifiers are matched; harmonize symbols/Ensembl IDs before combining cohorts.", "- RNA-seq and microarray effect estimates may differ by platform and design.", "- Results are research-only and are not diagnostic, prognostic, or treatment recommendations."]
    if PROVENANCE_FILE.exists():
        report_lines += ["", "### Bundled provenance", "```", PROVENANCE_FILE.read_text(encoding="utf-8").strip(), "```"]
    report_text = "\n".join(report_lines)
    st.download_button("Download reproducible evidence report (.md)", report_text.encode("utf-8"), "oncogene_evidence_report.md", "text/markdown")
    st.text_area("Report preview", report_text, height=360)


with gene_tab:
    st.subheader("Look up a gene across cohorts")
    query = st.text_input("Gene symbol / ID", placeholder="e.g., BRCA1, TP53, MKI67", key="gene_search").strip()
    if query:
        found=False
        for name, d in cohorts.items():
            matches=d[d["gene_id"].str.casefold() == query.casefold()]
            if not matches.empty:
                found=True
                for _, row in matches.iterrows():
                    with st.expander(name, expanded=True):
                        a,b,c=st.columns(3); a.metric("Gene",row["gene_id"]); b.metric("log₂ fold change",f"{row['logFC']:.3f}"); c.metric("Adjusted p-value",f"{row['adj.P.Val']:.3e}")
                        st.write("Direction in source contrast:","positive logFC" if row["logFC"]>0 else "negative logFC" if row["logFC"]<0 else "estimated change of zero")
        if len(cohorts)>1:
            match=combined[combined.gene_id.str.casefold()==query.casefold()]
            if not match.empty:
                st.subheader("Combined evidence (exploratory)"); st.dataframe(match,use_container_width=True,hide_index=True)
        if not found: st.warning("Gene/ID not found in the loaded result tables. Check identifier type and spelling.")
    st.caption("Differential expression alone does not validate a clinical biomarker.")

with results_tab:
    st.subheader(f"DEG results: {active_cohort_name}")
    st.caption("This table is based on precomputed differential-expression statistics. Changing filters does not refit limma/DESeq2.")
    view_mode = st.radio("Rows to display", ["Significant DEGs at current cutoffs", "All genes in active cohort"], horizontal=True, key="results_view_mode")
    result_base = filtered.copy() if view_mode == "Significant DEGs at current cutoffs" else active.copy()
    if view_mode == "All genes in active cohort":
        result_base["regulation"] = np.select(
            [result_base["logFC"] > 0, result_base["logFC"] < 0],
            ["Higher in tumor (if contrast is tumor vs normal)", "Lower in tumor (if contrast is tumor vs normal)"],
            default="No estimated change"
        )
    search_gene = st.text_input("Search gene symbol / ID", placeholder="e.g., BRCA1, TP53, ENSG...", key="results_gene_filter").strip()
    if search_gene:
        result_base = result_base[result_base["gene_id"].astype(str).str.contains(search_gene, case=False, regex=False)]
    sort_options = [c for c in ["adj.P.Val", "P.Value", "logFC", "gene_id"] if c in result_base.columns]
    sort_col = st.selectbox("Sort table by", sort_options, index=0, key="results_sort_col")
    ascending = st.checkbox("Ascending order", value=(sort_col in ["adj.P.Val", "P.Value", "gene_id"]), key="results_sort_ascending")
    columns = [c for c in ["gene_id","logFC","AveExpr","P.Value","adj.P.Val","t","B","regulation"] if c in result_base.columns]
    display = result_base[columns].sort_values(sort_col, ascending=ascending, na_position="last")
    st.write(f"Showing **{len(display):,}** rows. Download the full displayed table or all genes independently of the search filter.")
    st.dataframe(display, use_container_width=True, hide_index=True, height=520)
    st.download_button("Download displayed results (CSV)", display.to_csv(index=False).encode("utf-8"), "oncogene_displayed_results.csv", "text/csv")
    st.download_button("Download all active-cohort genes (CSV)", active.to_csv(index=False).encode("utf-8"), f"oncogene_{active_cohort_name}_all_genes.csv", "text/csv")
    st.download_button("Download current significant DEGs (CSV)", filtered.to_csv(index=False).encode("utf-8"), f"oncogene_{active_cohort_name}_filtered_degs.csv", "text/csv")

with samples_tab:
    st.subheader("Sample annotations")
    if META_FILE.exists():
        meta=load_metadata(META_FILE)
        if not meta.empty:
            st.write(f"Bundled GSE183947 sample annotation rows: {len(meta)}")
            st.dataframe(meta,use_container_width=True,hide_index=True)
            st.download_button("Download sample metadata (CSV)",meta.to_csv(index=False).encode("utf-8"),"oncogene_sample_metadata.csv","text/csv")
            st.warning("Confirm condition labels, sample counts, pairing, and sample-to-expression mapping against the original GEO records.")
        else: st.warning("Metadata file does not contain the expected columns.")
    st.info("Uploaded DEG tables do not include sample-level metadata. This app compares their precomputed gene-level statistics only.")

with methods_tab:
    st.subheader("Datasets and analysis")
    st.markdown("""
    ### Planned cohort strategy
    - **GSE183947** — bundled starting cohort; supplied limma results from log₂(FPKM + 1).
    - **GSE42568** — candidate independent validation cohort; this is an Affymetrix microarray study, not RNA-seq. Analyze with a platform-appropriate method and harmonize gene identifiers before cross-study comparison.
    - **TCGA-BRCA** — candidate large external cohort. Prefer a clearly defined tumor-versus-normal contrast and a reproducible data source/version.

    ### Biological interpretation modules
    - Therapeutic Target Evidence Bridge: source-traceable user-uploaded target evidence, matched to current candidate genes.
    - Tumor microenvironment signature screening: exploratory signed fold-change summaries; not sample-level deconvolution.
    - Downloadable Markdown evidence report with cohort configuration, consensus candidates, provenance, and limitations.
    - Optional local over-representation analysis from a user-supplied GMT gene-set library, with hypergeometric testing and Benjamini–Hochberg FDR.
    - Transparent candidate-prioritization heuristic exposing coverage, direction consistency, combined evidence, and descriptive effect magnitude.

    ### What this app computes
    - Optional matching of user-supplied therapeutic-target evidence to current consensus candidates.
    - Exploratory immune-signature summaries from a DEG table; these are not cell-infiltration estimates.
    - A downloadable report of current settings and available evidence.
    - Per-cohort volcano plots and filtered result tables from precomputed differential-expression statistics.
    - Gene-level comparison of log fold-change estimates across uploaded cohorts.
    - Exploratory cross-cohort evidence using Fisher's combined p-value, Benjamini–Hochberg FDR, and a separate direction-consistency flag.

    ### What this app does not do
    - It does **not** download GEO/TCGA data, process FASTQ files, fit limma/DESeq2 models, or perform a formal random-effects meta-analysis.
    - Fisher's method combines p-values, not effect sizes, and can be sensitive to cohort dependence and differing study designs. It should not be presented as a substitute for a properly specified meta-analysis.
    - Pathway enrichment depends on the uploaded GMT library, matching identifiers, gene universe, and statistical assumptions; enrichment is association, not proof of mechanism.
    - The prioritization score is a transparent heuristic and has not been trained or validated against clinical outcomes.
    - The app matches genes by exact identifier text. Map gene symbols/Ensembl IDs consistently and resolve duplicate mappings before uploading.

    ### Scientific limitations
    - Confirm original sample annotations, pairing, contrast direction, normalization, platform, and study methods before interpretation.
    - RNA-seq and microarray fold changes may not be directly comparable without harmonized analysis and careful interpretation.
    - Differential expression is not proof of causality, diagnostic performance, prognosis, or treatment response. Research/education only; not for clinical decision-making.
    """)
    if PROVENANCE_FILE.exists():
        with st.expander("Supplied GSE183947 analysis provenance"):
            st.code(PROVENANCE_FILE.read_text(encoding="utf-8"))

st.divider()
st.caption("OncoGene Explorer · Open-source research portfolio project · Verify source data before drawing conclusions")