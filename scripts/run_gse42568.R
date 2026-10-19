#!/usr/bin/env Rscript
# Reproducible GSE42568 processed-microarray differential-expression workflow.
# Output is a standardized gene-level DEG table for OncoGene Explorer.
# Does NOT combine expression matrices across studies.

options(stringsAsFactors = FALSE, timeout = 600)

required <- c("GEOquery", "limma", "Biobase")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install missing Bioconductor packages first: ", paste(missing, collapse = ", "),
       "\nRun: if (!requireNamespace('BiocManager')) install.packages('BiocManager'); BiocManager::install(c('GEOquery','limma','Biobase'))")
}

out_dir <- file.path("data", "GSE42568")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Downloading GSE42568 processed series matrix from NCBI GEO...")
gse_list <- GEOquery::getGEO("GSE42568", GSEMatrix = TRUE, getGPL = FALSE, destdir = out_dir)
if (length(gse_list) < 1) stop("No GEO ExpressionSet was returned for GSE42568.")
eset <- gse_list[[1]]
expr_probe <- Biobase::exprs(eset)
pheno <- Biobase::pData(eset)

if (ncol(expr_probe) != nrow(pheno)) stop("Expression/metadata sample dimensions do not match.")
if (is.null(colnames(expr_probe)) || !identical(colnames(expr_probe), rownames(pheno))) {
  common <- intersect(colnames(expr_probe), rownames(pheno))
  if (length(common) != ncol(expr_probe)) stop("Could not align expression columns with phenotype rows exactly.")
  expr_probe <- expr_probe[, common, drop = FALSE]
  pheno <- pheno[common, , drop = FALSE]
}

# GEO's official titles identify samples as 'Breast cancer' or 'Normal breast'.
text_fields <- intersect(c("title", "source_name_ch1", grep("characteristics_ch1", names(pheno), value = TRUE)), names(pheno))
if (!length(text_fields)) stop("No usable sample title/characteristics fields found in GEO metadata.")
sample_text <- apply(pheno[, text_fields, drop = FALSE], 1, function(z) paste(z, collapse = " | "))
condition <- ifelse(grepl("normal breast|normal tissue", sample_text, ignore.case = TRUE), "normal",
             ifelse(grepl("breast cancer|cancerous breast|tumou?r", sample_text, ignore.case = TRUE), "tumor", NA_character_))

meta <- data.frame(
  sample_id = rownames(pheno),
  title = if ("title" %in% names(pheno)) as.character(pheno$title) else sample_text,
  condition = condition,
  sample_text = sample_text,
  stringsAsFactors = FALSE
)
if (anyNA(meta$condition)) {
  write.csv(meta, file.path(out_dir, "sample_metadata_review.csv"), row.names = FALSE)
  stop("Could not assign condition for every sample. Review data/GSE42568/sample_metadata_review.csv; no DEG output was created.")
}
counts <- table(factor(meta$condition, levels = c("normal", "tumor")))
message("Sample groups: ", paste(names(counts), counts, sep = "=", collapse = ", "))
if (counts[["normal"]] != 17 || counts[["tumor"]] != 104) {
  write.csv(meta, file.path(out_dir, "sample_metadata_review.csv"), row.names = FALSE)
  stop("Unexpected sample counts for this GEO series (expected 17 normal and 104 tumor). Inspect metadata; no DEG output was created.")
}

# Use the GPL570 platform table for probe-to-symbol mapping.
message("Retrieving GPL570 probe annotation...")
gpl <- GEOquery::getGEO("GPL570", destdir = out_dir)
ann <- GEOquery::Table(gpl)
if (!"ID" %in% names(ann)) stop("GPL570 annotation lacks an ID column.")
symbol_col <- intersect(c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "gene_symbol", "Symbol"), names(ann))
if (!length(symbol_col)) {
  candidate <- grep("gene.?symbol|symbol", names(ann), ignore.case = TRUE, value = TRUE)
  if (!length(candidate)) stop("Could not identify gene-symbol field in GPL570 annotation. Inspect GPL570 table.")
  symbol_col <- candidate[[1]]
} else symbol_col <- symbol_col[[1]]
ann_small <- data.frame(probe_id = as.character(ann$ID), gene_id = as.character(ann[[symbol_col]]), stringsAsFactors = FALSE)
# Retain one unambiguous first symbol for annotation fields listing multiple symbols.
ann_small$gene_id <- trimws(sub("\\s*(///|//|;|\\|).*$", "", ann_small$gene_id))
ann_small$gene_id[is.na(ann_small$gene_id) | ann_small$gene_id %in% c("", "---", "NA")] <- NA_character_
ann_small <- ann_small[!is.na(ann_small$gene_id) & !duplicated(ann_small$probe_id), ]

probe_ids <- rownames(expr_probe)
map_idx <- match(probe_ids, ann_small$probe_id)
gene_ids <- ann_small$gene_id[map_idx]
keep <- !is.na(gene_ids) & nzchar(gene_ids) & apply(expr_probe, 1, function(x) all(is.finite(x)))
expr_probe <- expr_probe[keep, , drop = FALSE]
gene_ids <- gene_ids[keep]
if (nrow(expr_probe) < 1000) stop("Too few annotated probes remained; inspect platform mapping and input matrix.")

# Confirm the processed matrix is plausibly log2-scale. GEO series matrix values are typically GC-RMA/log2.
q <- as.numeric(stats::quantile(expr_probe, c(0.01, 0.5, 0.99), na.rm = TRUE))
if (q[3] > 100 || q[1] < -100) stop("Expression values are outside expected processed microarray range; review transformation before analysis.")

# Collapse multiple probes per gene using the median per sample, then fit limma at gene level.
# This avoids treating multiple probes for one gene as independent gene-level tests.
message("Collapsing probes to gene symbols and fitting limma model...")
gene_levels <- unique(gene_ids)
expr_gene <- vapply(gene_levels, function(g) apply(expr_probe[gene_ids == g, , drop = FALSE], 2, median, na.rm = TRUE), numeric(ncol(expr_probe)))
expr_gene <- t(expr_gene)
rownames(expr_gene) <- gene_levels
colnames(expr_gene) <- colnames(expr_probe)

condition <- factor(meta$condition, levels = c("normal", "tumor"))
design <- stats::model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
fit <- limma::lmFit(expr_gene, design)
contrast <- limma::makeContrasts(tumor - normal, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast), trend = TRUE, robust = TRUE)
res <- limma::topTable(fit2, number = Inf, sort.by = "none", adjust.method = "BH")
res$gene_id <- rownames(res)
res <- res[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
res <- res[order(res$adj.P.Val, res$P.Value), ]

# Metadata and provenance for auditability.
meta$condition <- as.character(meta$condition)
write.csv(meta, file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)
write.csv(res, file.path(out_dir, "GSE42568_limma_tumor_vs_normal_all_genes.csv"), row.names = FALSE)
write.csv(subset(res, adj.P.Val < 0.05), file.path(out_dir, "GSE42568_limma_tumor_vs_normal_FDR05.csv"), row.names = FALSE)
provenance <- c(
  "Dataset: GSE42568 (NCBI GEO)",
  "Platform: GPL570 Affymetrix Human Genome U133 Plus 2.0 Array",
  "Source: GEO processed series matrix; this is not a raw-CEL reprocessing workflow.",
  "Groups verified against GEO sample titles/characteristics; expected 104 breast cancer and 17 normal breast samples.",
  "Annotation: GPL570 platform table; ambiguous multi-symbol fields reduced to first listed symbol; unmapped probes removed.",
  "Preprocessing: gene-level median aggregation across probes per sample; input processed expression values retained on their supplied scale.",
  "Model: limma tumor versus normal, design ~ 0 + condition; empirical Bayes with trend=TRUE and robust=TRUE.",
  "Multiple testing: Benjamini-Hochberg FDR.",
  "Interpretation: independent unpaired cohort; not paired to GSE183947. Do not merge expression matrices across cohorts.",
  paste0("Samples: ", nrow(meta), "; normal=", sum(meta$condition == "normal"), "; tumor=", sum(meta$condition == "tumor")),
  paste0("Annotated genes tested: ", nrow(res)),
  paste0("Genes with FDR < 0.05: ", sum(res$adj.P.Val < 0.05, na.rm = TRUE)),
  paste0("Expression quantiles (1%, 50%, 99%): ", paste(signif(q, 4), collapse = ", ")),
  "Generated by scripts/run_gse42568.R; record R and Bioconductor package versions for publication-grade reproduction."
)
writeLines(provenance, file.path(out_dir, "analysis_provenance.txt"))
message("Finished. Standardized DEG file: ", file.path(out_dir, "GSE42568_limma_tumor_vs_normal_all_genes.csv"))
message("Upload that CSV to OncoGene Explorer as the GSE42568 cohort.")