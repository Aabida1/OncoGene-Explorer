# GSE42568 differential-expression workflow (Affymetrix GPL570)
# This script expects the processed series matrix and GPL annotation downloaded from GEO.
# Use limma on appropriately normalized microarray expression values, not DESeq2 on these values.

suppressPackageStartupMessages({library(GEOquery); library(limma)})

message('GSE42568 workflow starter. Review GEO sample characteristics and GPL570 annotation before setting group labels.')
message('Download GSE42568 series matrix and GPL570 platform annotation from NCBI GEO.')
message('This repository version deliberately does not infer tumour/normal labels automatically.')
stop('Before analysis, edit this script to point to verified local GEO files and confirm the sample group mapping in the metadata.')
