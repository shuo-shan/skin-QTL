#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

#### Helper functions ####
# ------------------------------#
# Helper: safe, fast scaling
# ------------------------------#
center_scale_vec <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sx <- stats::sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(NA_real_, length(x)))
  (x - mu) / sx
}

# Compute vector of correlations between r (length D) and each column of G (DxM)
# Returns numeric length M
cor_vec_mat <- function(r, G) {
  r_sc <- center_scale_vec(r)
  if (all(is.na(r_sc))) stop("Residual vector has zero variance or all NA.")
  # scale each SNP column (genotype) to mean0 sd1
  G_sc <- apply(G, 2, center_scale_vec)
  # apply() returns matrix with same dims only if M>1
  if (is.null(dim(G_sc))) G_sc <- matrix(G_sc, ncol = 1)
  
  # Correlation = mean(r_sc * g_sc) over donors with complete data
  # We'll do column-wise using matrix multiplication, but must handle NAs.
  # Simplest robust approach: restrict to donors complete for r and all SNPs is too strict.
  # Instead: handle per-SNP complete cases (still fast enough for ~1-5k SNPs).
  out <- rep(NA_real_, ncol(G_sc))
  for (j in seq_len(ncol(G_sc))) {
    ok <- is.finite(r_sc) & is.finite(G_sc[, j])
    if (sum(ok) < 4) { out[j] <- NA_real_; next }
    out[j] <- mean(r_sc[ok] * G_sc[ok, j])
  }
  out
}

# Adaptive permutation engine for one gene:
# returns list with p_hat, T_obs, perms_run, exceedances
adaptive_perm_gene <- function(r, G,
                               B0 = 200, B1 = 1000, B2 = 10000,
                               p_stop = 0.10, min_exceed_stop = 20,
                               seed = 1) {
  set.seed(seed)
  
  # Observed
  z_obs <- cor_vec_mat(r, G)
  T_obs <- max(abs(z_obs), na.rm = TRUE)
  
  # Stage runner
  run_stage <- function(B, b_prev = 0, T_prev = numeric(0)) {
    b <- b_prev
    T_perm <- T_prev
    D <- length(r)
    for (i in seq_len(B - length(T_perm))) {
      rp <- sample(r, size = D, replace = FALSE)
      z_p <- cor_vec_mat(rp, G)
      t_p <- max(abs(z_p), na.rm = TRUE)
      T_perm <- c(T_perm, t_p)
      if (is.finite(t_p) && t_p >= T_obs) b <- b + 1
    }
    list(b = b, T_perm = T_perm)
  }
  
  # Stage 0
  st0 <- run_stage(B0)
  b0 <- st0$b
  p0 <- (b0 + 1) / (B0 + 1)
  
  # Early stop if clearly null-ish
  # The "min_exceed_stop" rule prevents over-stopping due to randomness at small B
  if (b0 >= min_exceed_stop || p0 >= p_stop) {
    return(list(
      T_obs = T_obs, perms = B0, b = b0, p_hat = p0,
      stage = "stopped@B0", T_perm = st0$T_perm
    ))
  }
  
  # Stage 1
  st1 <- run_stage(B1, b_prev = b0, T_prev = st0$T_perm)
  b1 <- st1$b
  p1 <- (b1 + 1) / (B1 + 1)
  
  # If still not that small, stop
  if (p1 >= 0.01) {
    return(list(
      T_obs = T_obs, perms = B1, b = b1, p_hat = p1,
      stage = "stopped@B1", T_perm = st1$T_perm
    ))
  }
  
  # Stage 2 (optional for tiny p-values)
  st2 <- run_stage(B2, b_prev = b1, T_prev = st1$T_perm)
  b2 <- st2$b
  p2 <- (b2 + 1) / (B2 + 1)
  
  list(
    T_obs = T_obs, perms = B2, b = b2, p_hat = p2,
    stage = "stopped@B2", T_perm = st2$T_perm
  )
}

#### Main ####
# ------------------------------#
# Main: parse CLI args ####
# ------------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  cat("
Usage:
  Rscript adaptive_perm_gene_eQTL.R <ct> <condition> <QTLtype> <gene> <qc_gene_file> <geno_file> [out_prefix]

Example:
  Rscript adaptive_perm_gene_eQTL.R MEL TNF eQTL ERAP2 \\
    /.../results_QC/gene/modeling_stats_postQC_MEL_TNF_eQTL_ERAP2.txt \\
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
condition  <- args[[2]]
QTLtype    <- args[[3]]
gene       <- args[[4]]
qc_gene_file <- args[[5]]
geno_file  <- args[[6]]
vst_file   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
peer_file  <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_",ct,"_PBS-IFNG-IFNB-TNF.tsv")
pc_file    <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"
out_prefix <- if (length(args) >= 11) args[[11]] else paste(ct, condition, QTLtype, gene, sep = "_")

# toy example for debugging
ct         <- "MEL"
condition  <- "TNF"
QTLtype    <- "eQTL"
gene       <- "ERAP2"
qc_gene_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/results_QC/gene/modeling_stats_postQC_MEL_TNF_eQTL_ERAP2.txt"
geno_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/chunks/genotype_pairs_chunk_101.tsv"
vst_file   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
peer_file  <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_",ct,"_PBS-IFNG-IFNB-TNF.tsv")
pc_file    <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"
out_prefix <- if (length(args) >= 11) args[[11]] else paste(ct, condition, QTLtype, gene, sep = "_")

message("=== Adaptive permutation for ", ct, " ", condition, " ", QTLtype, " ", gene, " ===")

# ------------------------------#
# Step 1: SNP list for this gene ####
# ------------------------------#
qc <- read_tsv(qc_gene_file, col_names = FALSE, show_col_types = FALSE)
colnames(qc) <- c("snp","gene","dist","beta","p","se","n_geno0","n_geno1","n_geno2","qc_genotypes","z")
snps <- unique(qc$snp)
message("SNPs in QC file: ", length(snps))

# ------------------------------#
# Step 2: build donor-level expression for this gene in this condition ####
# ------------------------------#
# Read VST wide, keep gene and samples for this celltype+condition
vst <- read_tsv(vst_file, show_col_types = FALSE)

# Metadata
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(celltype == !!ct, condition == !!condition)

# Covariates
nPEER = 10
nGPC = 2

# PEERs (sample-level)
peer <- fread(peer_file) %>% dplyr::select(-c(PEER1, PEER2, PEER3, PEER4))
colnames(peer)[1] <- "sample"
peer <- peer[ , 1:(nPEER+1)]

# genotype PCs (donor-level)
pcs <- fread(pc_file) %>% dplyr::select(-c(FID, donor_num))
pcs <- pcs[ , 1:(nGPC+1)]

# Subset VST to this gene only
vst_gene <- vst %>% filter(gene == !!gene)

# Keep only columns that exist in meta$sample
sample_cols <- intersect(setdiff(names(vst_gene), "gene"), meta$sample)
if (length(sample_cols) == 0) stop("No overlapping TNF samples between VST and metadata for this ct/condition.")
vst_gene <- vst_gene %>% select(gene, all_of(sample_cols))

# Make long
expr_donor <- vst_gene %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expr") %>%
  inner_join(meta %>% select(sample, donor), by = "sample") %>%
  inner_join(peer, by = "sample") %>%
  left_join(pcs, by = "donor")

message("Donors with ", condition, " expression: ", nrow(expr_donor))
colnames(expr_donor)
# ------------------------------#
# Step 3: residualize expression w.r.t covariates (NO genotype) ####
# ------------------------------#
# Build formula: expr ~ PEERs + PCs
peer_use <- colnames(peer[ , -1])
pc_use <- colnames(pcs[ , -1])
covars <- c(peer_use, pc_use)

form <- as.formula(paste("expr ~", paste(covars, collapse = " + ")))
fit0 <- lm(form, data = expr_donor)
expr_donor$resid <- resid(fit0)

# Residual vector r aligned by donor
r <- expr_donor$resid
names(r) <- expr_donor$donor

# ------------------------------#
# Step 4: genotype matrix for those SNPs and donors ####
# ------------------------------#
geno <- read_tsv(geno_file, show_col_types = FALSE)

fixed_cols <- c("CHROM","POS","ID","REF","ALT")
# Identify genotype sample columns
old_sample_cols <- setdiff(colnames(geno), fixed_cols)

# Create NEW cleaned column names (donor IDs)
new_sample_cols <- old_sample_cols
new_sample_cols <- gsub("skineQTL-", "", new_sample_cols)
new_sample_cols <- gsub("^F0", "F", new_sample_cols)
setnames(geno, old_sample_cols, new_sample_cols)

# Now work with the renamed columns from here forward:
geno_sample_cols <- new_sample_cols
if (length(geno_sample_cols) == 0) stop("No sample GT columns found in GENOTYPE.")

# Keep only genotype donors that appear in expression data
donor_cols_present <- intersect(geno_sample_cols, unique(expr_donor$donor))
geno <- geno[, c(fixed_cols, donor_cols_present), with = FALSE]
if (length(donor_cols_present) < 5) stop("Too few donor genotype columns overlap with expression donors.")

# select SNPs in the gene's cis-window
geno_sub <- geno %>%
  filter(ID %in% snps) %>%
  select(ID, all_of(donor_cols_present))

message("SNPs found in genotype file: ", nrow(geno_sub))

# Build G: donors x SNPs
# Ensure donor order matches r
G <- geno_sub %>%
  arrange(match(ID, snps)) %>%
  select(-ID) %>%
  as.matrix()

# geno_sub columns are donors; ensure same order as donors in r
# (we already used donor_cols_present derived from expr donors; now align exactly)
G <- geno_sub %>%
  arrange(match(ID, snps)) %>%
  select(all_of(donors)) %>%
  as.matrix()

# transpose to donors x SNPs
G <- t(G)
storage.mode(G) <- "numeric"

# ------------------------------#
# Step 5: adaptive permutation ####
# ------------------------------#
res <- adaptive_perm_gene(
  r = r, G = G,
  B0 = 200, B1 = 1000, B2 = 10000,
  p_stop = 0.10, min_exceed_stop = 20,
  seed = 1
)

# ------------------------------#
# Step 6: write outputs
# ------------------------------#
summary_out <- tibble(
  celltype = ct,
  condition = condition,
  QTLtype = QTLtype,
  gene = gene,
  n_snps = ncol(G),
  n_donors = length(r),
  T_obs = res$T_obs,
  perms = res$perms,
  exceed = res$b,
  p_emp = res$p_hat,
  stage = res$stage
)

write_tsv(summary_out, paste0(out_prefix, ".perm_summary.tsv"))
message("Wrote: ", paste0(out_prefix, ".perm_summary.tsv"))

#write_tsv(tibble(T_perm = res$T_perm), paste0(out_prefix, ".Tperm.tsv"))
#message("Wrote: ", paste0(out_prefix, ".Tperm.tsv"))
#print(summary_out)
