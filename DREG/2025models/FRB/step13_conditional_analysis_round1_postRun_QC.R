# QTL modeling QC filtering for conditional_analysis_round1 output
Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(future.apply)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript step13_conditional_analysis_round1_postRun_QC.R <CELLTYPE> <condition> <QTLtype> <chunk_id>")
}

# take in arguments
ct <- args[1]              # MEL, KRT, FRB
condition_arg <- args[2]   # PBS, IFNG, IFNB, TNF
QTLtype_arg <- args[3]     # eQTL, reQTL
chunk_id <- args[4]        # "000", "001", etc

# define variables
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)

# NEW input path
model_stats_file <- paste0(
  dir,
  "/conditional_analysis_round1/results/",
  condition_arg, "/",
  QTLtype_arg, "/result_",
  chunk_id, ".tsv"
)

# # toy example for debugging
# ct <- "FRB"
# condition_arg <- "IFNB"
# QTLtype_arg <- "reQTL"
# chunk_id <- "219"
# dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)
# model_stats_file <- paste0(
#   dir,
#   "/conditional_analysis_round1/results/",
#   condition_arg, "/",
#   QTLtype_arg, "/result_",
#   chunk_id, ".tsv"
# )

message("====== QTL statistics QC start ======")
message("Model Statistics Chunk: ", model_stats_file)

# -----------------------------
# read new-format result file
# -----------------------------
if (!file.exists(model_stats_file)) {
  warning("Input file does not exist: ", model_stats_file)
  stats_all <- data.table()
} else {
  stats_all <- tryCatch(
    fread(model_stats_file, header = TRUE),
    error = function(e) {
      warning("Failed to read file: ", model_stats_file)
      data.table()
    }
  )
}

# handle empty input
if (nrow(stats_all) == 0) {
  message("Input file is empty. Writing empty QC output.")
  
  outdir <- paste0(
    dir,
    "/conditional_analysis_round1/results_QC/",
    condition_arg, "/", QTLtype_arg
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  outfile <- paste0(
    outdir,
    "/modeling_stats_postQC_",
    ct, "_", condition_arg, "_", QTLtype_arg, "_", chunk_id, ".txt"
  )
  
  # define expected output columns
  empty_dt <- data.table(
    snp = character(),
    gene = character(),
    dist = numeric(),
    beta = numeric(),
    p = numeric(),
    se = numeric(),
    n_geno0 = integer(),
    n_geno1 = integer(),
    n_geno2 = integer(),
    qc_genotypes = integer(),
    z = numeric()
  )
  
  fwrite(empty_dt,
         file = outfile,
         quote = FALSE, sep = "\t",
         row.names = FALSE, col.names = TRUE)
  
  message("Done (empty output): ", outfile)
  quit(save = "no")
}

required_cols <- c(
  "snp", "gene", "celltype", "condition", "QTLtype",
  "dist_to_TSS", "beta", "p", "se"
)
missing_cols <- setdiff(required_cols, colnames(stats_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
}

# sanity check: file content should match args
if (any(stats_all$celltype != ct)) {
  stop("Input file celltype does not match argument ct")
}
if (any(stats_all$condition != condition_arg)) {
  stop("Input file condition does not match argument condition_arg")
}
if (any(stats_all$QTLtype != QTLtype_arg)) {
  stop("Input file QTLtype does not match argument QTLtype_arg")
}

# =========================================================
# eQTL block: use qc string already present in new result file
# =========================================================
if (QTLtype_arg == "eQTL") {
  
  needed_qc_cols <- c("qc_genotypes", "qc_n_per_condition_by_genotype")
  missing_qc_cols <- setdiff(needed_qc_cols, colnames(stats_all))
  if (length(missing_qc_cols) > 0) {
    stop("Missing QC columns for eQTL input file: ", paste(missing_qc_cols, collapse = ", "))
  }
  
  # split qc string into 12 integer columns
  # order:
  # g0_PBS g0_IFNG g0_IFNB g0_TNF
  # g1_PBS g1_IFNG g1_IFNB g1_TNF
  # g2_PBS g2_IFNG g2_IFNB g2_TNF
  qc_mat <- tstrsplit(stats_all$qc_n_per_condition_by_genotype, "_", fixed = TRUE)
  qc_mat <- lapply(qc_mat, as.integer)
  
  conds <- c("PBS", "IFNG", "IFNB", "TNF")
  qc_names <- c(
    paste0("n_g0_", conds),
    paste0("n_g1_", conds),
    paste0("n_g2_", conds)
  )
  
  if (length(qc_mat) != 12) {
    stop("qc_n_per_condition_by_genotype does not split into 12 fields as expected.")
  }
  
  stats_all[, (qc_names) := qc_mat]
  
  stats_relevant <- copy(stats_all)
  stats_relevant[, dist := dist_to_TSS]
  
  col_g0 <- paste0("n_g0_", condition_arg)
  col_g1 <- paste0("n_g1_", condition_arg)
  col_g2 <- paste0("n_g2_", condition_arg)
  
  stats_relevant[, n_geno0 := .SD[[1]], .SDcols = col_g0]
  stats_relevant[, n_geno1 := .SD[[1]], .SDcols = col_g1]
  stats_relevant[, n_geno2 := .SD[[1]], .SDcols = col_g2]
  
  CPM_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
  CPM <- data.table::fread(CPM_FILE) %>%
    dplyr::select(-name, -gene) %>%
    dplyr::filter(final_gene != "ARMCX5-GPRASP2") %>%
    tibble::column_to_rownames("final_gene")
  
  pattern <- paste0(ct, "_", condition_arg, "_")
  if (condition_arg == "IFNG") {
    pattern <- paste0(ct, "_(", "IFN|IFNG", ")_")
  }
  selected_columns <- grep(pattern, names(CPM), value = TRUE)
  CPM <- CPM[, selected_columns, drop = FALSE]
  
  cpm_thresh <- 1
  min_frac <- 0.25
  n_samples <- ncol(CPM)
  min_n <- ceiling(n_samples * min_frac)
  
  expressed_genes <- rownames(CPM)[rowSums(CPM >= cpm_thresh, na.rm = TRUE) >= min_n]
}

# =========================================================
# reQTL block: preserve genotype-table join and re-count logic
# =========================================================
if (QTLtype_arg == "reQTL") {
  
  chunk_file <- paste0(dir, "/chunks/pairs_chunk_", chunk_id, ".tsv")
  GENO_FILE  <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
  META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
  CPM_FILE   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
  conds <- c("PBS","IFNG","IFNB","TNF")
  
  # load pairs from chunk, mainly to keep old workflow structure/context
  pairs <- fread(chunk_file, header = TRUE)
  colnames(pairs)[4] <- "gene"
  colnames(pairs)[11] <- "snp"
  
  # metadata: only PBS + selected stimulation
  selected_cond <- unique(c("PBS", condition_arg))
  meta <- fread(META_FILE) %>%
    dplyr::filter(celltype == ct) %>%
    dplyr::filter(condition %in% selected_cond)
  
  meta$condition <- factor(as.character(meta$condition), levels = conds)
  
  donor_cond_counts <- as.data.frame(table(meta$donor))
  paired_donors <- donor_cond_counts %>%
    dplyr::filter(Freq == 2) %>%
    pull(Var1) %>%
    as.character()
  
  # load genotype table
  geno_dt <- fread(GENO_FILE)
  fixed_cols <- c("CHROM","POS","ID","REF","ALT")
  
  # identify genotype sample columns
  old_sample_cols <- setdiff(colnames(geno_dt), fixed_cols)
  
  # clean sample names to donor IDs
  new_sample_cols <- old_sample_cols
  new_sample_cols <- gsub("skineQTL-", "", new_sample_cols)
  new_sample_cols <- gsub("^F0", "F", new_sample_cols)
  setnames(geno_dt, old_sample_cols, new_sample_cols)
  
  geno_sample_cols <- new_sample_cols
  
  # keep only donors present in metadata
  geno_keep_samples <- intersect(geno_sample_cols, unique(meta$donor))
  geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples), with = FALSE]
  
  # for reQTL: keep only paired donors
  geno_keep_samples_paired <- intersect(paired_donors, geno_keep_samples)
  geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples_paired), with = FALSE]
  
  # use new-format result file directly
  stats_relevant <- copy(stats_all)
  stats_relevant[, dist := dist_to_TSS]
  
  # join genotype table by SNP ID
  stats_relevant <- stats_relevant %>%
    left_join(geno_dt, by = c("snp" = "ID"))
  
  donor_cols <- intersect(geno_keep_samples_paired, colnames(stats_relevant))
  
  if (length(donor_cols) == 0) {
    stop("No paired donor genotype columns found after joining genotype table.")
  }
  
  # recompute genotype counts from paired genotype table
  stats_relevant <- stats_relevant %>%
    mutate(
      n_geno0 = rowSums(across(all_of(donor_cols), ~ .x == 0), na.rm = TRUE),
      n_geno1 = rowSums(across(all_of(donor_cols), ~ .x == 1), na.rm = TRUE),
      n_geno2 = rowSums(across(all_of(donor_cols), ~ .x == 2), na.rm = TRUE),
      qc_genotypes = (n_geno0 > 0) + (n_geno1 > 0) + (n_geno2 > 0)
    )
  
  # get expressed genes using paired donors only
  kept.samples <- meta %>%
    dplyr::filter(donor %in% donor_cols) %>%
    pull(sample)
  
  CPM <- data.table::fread(CPM_FILE) %>%
    dplyr::select(-name, -gene) %>%
    dplyr::filter(final_gene != "ARMCX5-GPRASP2") %>%
    tibble::column_to_rownames("final_gene") %>%
    dplyr::select(all_of(kept.samples))
  
  CPM.pbs <- CPM[, grep(paste0(ct, "_PBS_"), names(CPM), value = TRUE), drop = FALSE]
  
  pattern <- paste0(ct, "_", condition_arg, "_")
  if (condition_arg == "IFNG") {
    pattern <- paste0(ct, "_(", "IFN|IFNG", ")_")
  }
  selected_columns <- grep(pattern, names(CPM), value = TRUE)
  CPM.stim <- CPM[, selected_columns, drop = FALSE]
  
  cpm_thresh <- 1
  min_frac <- 0.25
  n_samples <- ncol(CPM.pbs)
  min_n <- ceiling(n_samples * min_frac)
  
  expressed_genes.pbs <- rownames(CPM.pbs)[rowSums(CPM.pbs >= cpm_thresh, na.rm = TRUE) >= min_n]
  expressed_genes.stim <- rownames(CPM.stim)[rowSums(CPM.stim >= cpm_thresh, na.rm = TRUE) >= min_n]
  expressed_genes <- union(expressed_genes.pbs, expressed_genes.stim)
}

# =========================================================
# apply QC filters
# =========================================================

# Rule A: require all three genotypes globally
stats_relevant.QCa <- stats_relevant %>% dplyr::filter(qc_genotypes == 3)

# Rule B: minimum donors per genotype is 3
# and total donors/samples >= 15
stats_relevant.QCab <- stats_relevant.QCa %>%
  dplyr::filter(pmin(n_geno0, n_geno1, n_geno2, na.rm = TRUE) >= 3) %>%
  dplyr::filter((n_geno0 + n_geno1 + n_geno2) >= 15)

# Rule C: gene is expressed
stats_relevant.QCabc <- stats_relevant.QCab %>%
  dplyr::filter(gene %in% expressed_genes)

# Rule D: z-score triage
selected_genes <- stats_relevant.QCabc %>%
  dplyr::mutate(z = beta / se) %>%
  dplyr::filter(abs(z) >= 2.5) %>%
  dplyr::pull(gene) %>%
  unique()

stats_relevant.QCabcd <- stats_relevant.QCabc %>%
  dplyr::filter(gene %in% selected_genes) %>%
  dplyr::mutate(z = beta / se)

# =========================================================
# write output
# =========================================================
outdir <- paste0(
  dir,
  "/conditional_analysis_round1/results_QC/",
  condition_arg, "/", QTLtype_arg
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- paste0(
  outdir,
  "/modeling_stats_postQC_",
  ct, "_", condition_arg, "_", QTLtype_arg, "_", chunk_id, ".txt"
)

data.table::fwrite(
  stats_relevant.QCabcd,
  file = outfile,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

message("Done QC filtering for: ", outfile)
