#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qpdf)
})

# ---------------------- Argument Parsing ----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("
Usage:
  Rscript step6_combine_all_plots.R <dir>

Example:
  Rscript step6_combine_all_plots.R /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots/temp_output/MEL_IFNG_eQTL_gene_fdr05 \\
", "\n")
  quit(save = "no", status = 1)
}

pdf_dir         <- args[[1]]
basename <- gsub(".*/","",pdf_dir)
output_dir <- dirname(dirname(pdf_dir))
output_file <- paste0(output_dir, "/combined_", basename, ".pdf")
pdf_files <- list.files(pdf_dir, pattern = "\\.pdf$", full.names = TRUE)

# Parse gene and snp from filename: plot_<gene>_<snp>.pdf
base <- basename(pdf_files)
m <- regexec("^plot_([^_]+)_([^_]+)\\.pdf$", base)
parts <- regmatches(base, m)

# keep only files that match the pattern
ok <- lengths(parts) == 3
pdf_files <- pdf_files[ok]
parts <- parts[ok]

# sort pdf files by gene name then by snp
gene <- vapply(parts, `[`, character(1), 2)
snp  <- vapply(parts, `[`, character(1), 3)
ord <- order(gene, snp)
pdf_sorted <- pdf_files[ord]

# generate output
qpdf::pdf_combine(input = pdf_sorted, output = output_file)
message(paste0("generated combined output ",output_file))
