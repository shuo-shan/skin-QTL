#!/usr/bin/env Rscript
# =============================================================================
# step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R
#
# PURPOSE:
#   For a given gene, test whether the PBS eQTL and a cytokine-stimulated eQTL
#   are driven by the same causal variant (H4) or distinct causal variants (H3).
#
#   This is the core analysis for my hypothesis that cytokine stimulation engages
#   different regulatory architecture than baseline expression. We use coloc.susie
#   because it handles multiple independent signals per locus (unlike standard
#   coloc.abf which assumes one signal per locus). coloc.abf is also run in
#   parallel as a sanity check and fallback.
#
# DESIGN NOTE — OVERLAPPING SAMPLES:
#   PBS and cytokine eQTLs are from the same donors. Standard coloc assumes
#   independent datasets, so sample overlap biases toward H4 (shared signal).
#   This means my H3 findings are CONSERVATIVE — the true distinct-signal rate
#   is likely even higher. I note this as a limitation in my methods.
#
# USAGE:
#   Rscript step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R <ct> <gene> <cytokine>
#
# EXAMPLES:
#   Rscript step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R FRB ERAP2 IFNG
#   Rscript step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R FRB IRF3 IFNG
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
  library(coloc)    # coloc 5.2.3 — provides coloc.susie and coloc.abf
  library(susieR)   # susie_rss — fine-mapping with summary statistics
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =============================================================================
# 0. CLI ARGUMENTS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("
Usage:
  Rscript step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R <ct> <gene> <cytokine>

Examples:
  Rscript step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R FRB ERAP2 IFNG
  Rscript step12_compare_PBSeQTL_cytokineeQTL_susie.coloc_perGene.R MEL ERAP2 IFNB
\n")
  quit(save = "no", status = 1)
}

ct       <- args[[1]]   # cell type, e.g. FRB (fibroblast) or MEL (melanocyte)
this_gene <- args[[2]]  # gene name, e.g. ERAP2
cytokine  <- args[[3]]  # stimulation condition, e.g. IFNG, IFNB, TNF

# # Uncomment for interactive testing:
# ct        <- "FRB"
# this_gene <- "IRF3"
# cytokine  <- "IFNG"

# =============================================================================
# 1. GLOBAL PATHS
# =============================================================================
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)

# Each gene belongs to a chunk — this maps gene -> chunk ID for fast file lookup
chunk_id_lookup <- data.table::fread(paste0(dir, "/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[chunk_id_lookup$gene == this_gene, ]$chunk)
chunk_id <- sprintf("%03d", chunk_id)
rm(chunk_id_lookup)

pair_file       <- paste0(dir, "/chunks/pairs_chunk_", chunk_id, ".tsv")
geno_file       <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
vst_file        <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
modelstats_file <- paste0(dir, "/results/result_", chunk_id, ".tsv")
meta_file       <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"

# Covariate files — needed to residualize genotypes before computing LD
# PEER factors are run jointly across all conditions (PBS+IFNG+IFNB+TNF),
# so this file contains one row per sample across all conditions.
# PEER1-4 are condition indicator dummies (not real latent factors) — skip them.
# Real PEER factors start at PEER5. I use 10 PEERs (PEER5-PEER14) + 2 genotype PCs.
PEER_FILE     <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_", ct, "_PBS-IFNG-IFNB-TNF.tsv")
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"
nPEER <- 10   # number of PEER factors used in eQTL modeling (PEER5-PEER14)
nGPC  <- 2    # number of genotype PCs used in eQTL modeling

# Output directory — one subfolder per cytokine for easy organization
out_dir    <- paste0(dir, "/coloc_susie/", cytokine)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_prefix <- paste("coloc_susie", this_gene, sep = "_")

# =============================================================================
# 2. UTILITY FUNCTIONS
# =============================================================================

# fix_p_zeros: replace exact zero p-values with half the minimum nonzero value.
# Zero p-values arise from numerical underflow in the linear model (p too small
# to represent as float). Coloc.abf uses p-values internally and cannot handle
# zeros — this prevents log(0) crashes without distorting the signal.
fix_p_zeros <- function(p) {
  p <- as.numeric(p)
  if (!any(is.finite(p) & p > 0, na.rm = TRUE)) return(p)
  min_nonzero <- min(p[p > 0], na.rm = TRUE)
  if (!is.finite(min_nonzero)) return(p)
  p[p == 0] <- min_nonzero / 2
  p
}

# =============================================================================
# 3. LOAD SHARED DATA
# =============================================================================
# Load pairs, genotypes, and metadata once — shared across PBS and cytokine

pairs <- fread(pair_file, header = TRUE) %>%
  dplyr::filter(gene_name == this_gene) %>%
  dplyr::mutate(
    SNP = stringr::str_to_lower(SNP_ID),
    key = paste0(gene_name, "_", SNP_ID)
  )

# Genotype matrix — SNPs x donors, dosage coded 0/1/2
# Rows = SNPs, columns = donor IDs. Used for LD computation.
genotype_all <- fread(geno_file) %>%
  dplyr::mutate(ID = stringr::str_to_lower(ID)) %>%
  dplyr::filter(ID %in% pairs$SNP_ID)

meta_all <- readr::read_tsv(meta_file, show_col_types = FALSE)

# Genomic coordinates — used for Manhattan plot x-axis
chr            <- unique(pairs$gene_chr)
ciswindow_left  <- unique(pairs$gene_start) - 500000
ciswindow_right <- unique(pairs$gene_start) + 500000

# =============================================================================
# 4. PBS eQTL SUMMARY STATISTICS
# =============================================================================
kept_condition <- "PBS"

meta_PBS   <- meta_all %>%
  dplyr::filter(celltype == ct, condition == kept_condition) %>%
  dplyr::select(sample, donor, condition)
donors_PBS <- meta_PBS$donor

# sdY: standard deviation of expression for this gene in PBS samples.
# Required by coloc for type="quant" when using beta/varbeta input format.
# This is preferable to the pvalue+MAF+N format as it is more accurate.
y_PBS        <- data.table::fread(vst_file, showProgress = FALSE) %>%
  dplyr::filter(gene == this_gene) %>%
  dplyr::select(all_of(meta_PBS$sample)) %>%
  t() %>% as.data.frame() %>% dplyr::pull(V1)
sdY_eqtl_PBS <- stats::sd(y_PBS, na.rm = TRUE)

# Load modeling stats and extract PBS-specific columns.
# IMPORTANT: grep patterns use "^eQTL_" prefix to avoid matching "reQTL_*" columns
# that also exist in the result file (reQTL = response eQTL, different model).
modelstats_all <- data.table::fread(modelstats_file, showProgress = FALSE)

fixed_cols         <- c("snp", "gene", "dist")
stats_cols         <- setdiff(colnames(modelstats_all), fixed_cols)
selected_PBS_cols  <- grep(paste0("^eQTL.*", kept_condition, "$"), stats_cols, value = TRUE)

modelstats_PBS <- modelstats_all[, c(fixed_cols, selected_PBS_cols), with = FALSE] %>%
  dplyr::mutate(
    snp = stringr::str_to_lower(snp),
    key = paste0(gene, "_", snp)
  ) %>%
  dplyr::filter(key %in% pairs$key)

# Identify beta/se/p columns — strict "^eQTL_" prefix prevents reQTL contamination
beta_col_PBS <- grep("^eQTL_beta_", colnames(modelstats_PBS), value = TRUE)
se_col_PBS   <- grep("^eQTL_se_",   colnames(modelstats_PBS), value = TRUE)
p_col_PBS    <- grep("^eQTL_p_",    colnames(modelstats_PBS), value = TRUE)

df_eqtl_PBS <- modelstats_PBS %>%
  dplyr::transmute(
    rsid            = snp,
    gene            = gene,
    dist            = as.numeric(dist),
    beta_qtl_PBS    = as.numeric(.data[[beta_col_PBS]]),
    varbeta_qtl_PBS = as.numeric(.data[[se_col_PBS]])^2,
    p_qtl_PBS       = fix_p_zeros(as.numeric(.data[[p_col_PBS]]))
  )

# =============================================================================
# 5. CYTOKINE eQTL SUMMARY STATISTICS
# =============================================================================
kept_condition <- cytokine

meta_cyto   <- meta_all %>%
  dplyr::filter(celltype == ct, condition == kept_condition) %>%
  dplyr::select(sample, donor, condition)
donors_cyto <- meta_cyto$donor

# Note: cytokine sample sizes are smaller than PBS due to cell death upon
# stimulation. This reduces power for fine-mapping but does not invalidate
# the analysis — it just means credible sets will be wider.
y_cyto        <- data.table::fread(vst_file, showProgress = FALSE) %>%
  dplyr::filter(gene == this_gene) %>%
  dplyr::select(all_of(meta_cyto$sample)) %>%
  t() %>% as.data.frame() %>% dplyr::pull(V1)
sdY_eqtl_cyto <- stats::sd(y_cyto, na.rm = TRUE)

selected_cyto_cols <- grep(paste0("^eQTL.*", kept_condition, "$"), stats_cols, value = TRUE)

modelstats_cyto <- modelstats_all[, c(fixed_cols, selected_cyto_cols), with = FALSE] %>%
  dplyr::mutate(
    snp = stringr::str_to_lower(snp),
    key = paste0(gene, "_", snp)
  ) %>%
  dplyr::filter(key %in% pairs$key)

beta_col_cyto <- grep("^eQTL_beta_", colnames(modelstats_cyto), value = TRUE)
se_col_cyto   <- grep("^eQTL_se_",   colnames(modelstats_cyto), value = TRUE)
p_col_cyto    <- grep("^eQTL_p_",    colnames(modelstats_cyto), value = TRUE)

df_eqtl_cyto <- modelstats_cyto %>%
  dplyr::transmute(
    rsid             = snp,
    gene             = gene,
    dist             = as.numeric(dist),
    beta_qtl_cyto    = as.numeric(.data[[beta_col_cyto]]),
    varbeta_qtl_cyto = as.numeric(.data[[se_col_cyto]])^2,
    p_qtl_cyto       = fix_p_zeros(as.numeric(.data[[p_col_cyto]]))
  )

# =============================================================================
# 6. GENOTYPE MATRIX AND VARIANT TABLE
# =============================================================================
fixed_col <- c("CHROM", "POS", "ID", "REF", "ALT")
donor_col <- setdiff(colnames(genotype_all), fixed_col)

# geno_mat: SNPs x donors, all donors in the chunk (not yet condition-filtered)
# Condition-specific subsetting happens during residualization below.
geno_mat <- genotype_all %>%
  dplyr::select(all_of(donor_col)) %>%
  dplyr::mutate(dplyr::across(everything(), as.numeric)) %>%
  as.matrix()
rownames(geno_mat) <- stringr::str_to_lower(genotype_all$ID)

alt_af  <- rowMeans(geno_mat, na.rm = TRUE) / 2
maf_qtl <- pmin(alt_af, 1 - alt_af)
n_called <- rowSums(!is.na(geno_mat))

maf_tbl <- tibble::tibble(
  rsid       = names(alt_af) %>% stringr::str_to_lower(),
  AF_ALT_QTL = unname(alt_af),
  MAF_QTL    = unname(maf_qtl),
  N_CALLED   = unname(n_called)
)

variant_tbl <- genotype_all %>%
  dplyr::select(POS, ID, REF, ALT) %>%
  dplyr::mutate(
    rsid = stringr::str_to_lower(ID),
    REF  = stringr::str_to_upper(REF),
    ALT  = stringr::str_to_upper(ALT)
  ) %>%
  dplyr::left_join(maf_tbl, by = "rsid")

# Quick summary
cat(sprintf("[%s | %s | %s] Data loaded.\n", ct, this_gene, cytokine))
cat(sprintf("  PBS SNPs: %d | Cyto SNPs: %d\n", nrow(df_eqtl_PBS), nrow(df_eqtl_cyto)))

# =============================================================================
# 7. INTERSECT SNP SET
# =============================================================================
# Only run coloc on SNPs present in BOTH conditions.
# SNPs absent from cytokine due to MAF/allele count QC filtering (from smaller
# sample size) are excluded. This is a known limitation — PBS lead SNPs in rare
# LD blocks may be absent from the cytokine matrix, making low-r2 between lead
# SNPs ambiguous (could be distinct biology OR just SNP unavailability).
# Restricting to shared SNPs is the conservative and correct approach for coloc.

shared_rsids <- intersect(df_eqtl_PBS$rsid, df_eqtl_cyto$rsid)

cat(sprintf("  Shared SNPs: %d (PBS=%d, cyto=%d, dropped=%d)\n",
            length(shared_rsids),
            nrow(df_eqtl_PBS), nrow(df_eqtl_cyto),
            nrow(df_eqtl_PBS) + nrow(df_eqtl_cyto) - 2 * length(shared_rsids)))

df_PBS <- df_eqtl_PBS  %>% dplyr::filter(rsid %in% shared_rsids) %>% dplyr::arrange(rsid)
df_cyto <- df_eqtl_cyto %>% dplyr::filter(rsid %in% shared_rsids) %>% dplyr::arrange(rsid)
vt      <- variant_tbl  %>% dplyr::filter(rsid %in% shared_rsids) %>% dplyr::arrange(rsid)

# Enforce identical SNP ordering across all three data frames — required for
# correct LD matrix alignment with summary statistics
stopifnot(identical(df_PBS$rsid, df_cyto$rsid))
stopifnot(identical(df_PBS$rsid, vt$rsid))

# =============================================================================
# 8. RESIDUALIZE GENOTYPES AND COMPUTE LD MATRICES
# =============================================================================
# WHY RESIDUALIZE:
#   My eQTL model includes PEER factors and genotype PCs as covariates.
#   The LD matrix passed to susie_rss should reflect the effective correlation
#   structure AFTER covariate correction — i.e., computed from residualized
#   genotypes. Using raw genotype correlation causes z-score vs LD mismatch,
#   which makes SuSiE's prior variance estimator blow up.
#
# WHY CONDITION-SPECIFIC LD:
#   PBS and cytokine conditions have different donor subsets (cell death causes
#   dropout in stimulated conditions). The LD matrix must be computed from the
#   exact donors used in each condition's eQTL model. Using PBS donors to
#   compute IFNG LD was the root cause of earlier SuSiE failures.
#
# WHY PEER FACTORS ARE CONDITION-SPECIFIC:
#   PEER was run jointly on all conditions, but each sample gets its own PEER
#   values. PEER1-4 are condition indicator dummies (skip them). Real PEER
#   latent factors start at PEER5. I use PEER5-PEER14 (10 factors) + 2 geno PCs.

peers <- data.table::fread(PEER_FILE)
colnames(peers)[1] <- "sample"
all_peer_cols <- setdiff(colnames(peers), "sample")
peer_use      <- all_peer_cols[5:(5 + nPEER - 1)]  # PEER5 through PEER14

gpcs   <- data.table::fread(GENO_PCS_FILE)
pc_use <- setdiff(colnames(gpcs), c("FID", "donor", "donor_num"))[seq_len(nGPC)]

cat(sprintf("  Covariates: %s + %s\n",
            paste(peer_use, collapse = ","),
            paste(pc_use, collapse = ",")))

# build_covar_matrix: assemble donor x covariate matrix for one condition.
# Includes intercept column (required for QR projection).
# Only donors present in BOTH metadata AND PEER file are retained —
# this matches the exact set used in the eQTL model.
build_covar_matrix <- function(condition_name) {
  meta_cond <- meta_all %>%
    dplyr::filter(celltype == ct, condition == condition_name) %>%
    dplyr::select(sample, donor)
  
  peer_cond <- as.data.frame(peers) %>%
    dplyr::filter(sample %in% meta_cond$sample) %>%
    dplyr::select(sample, all_of(peer_use)) %>%
    dplyr::left_join(meta_cond, by = "sample")
  
  n_dropped <- nrow(meta_cond) - nrow(peer_cond)
  if (n_dropped > 0) {
    cat(sprintf("  [%s] %d sample(s) dropped (not in PEER file)\n",
                condition_name, n_dropped))
  }
  
  pc_cond <- as.data.frame(gpcs) %>%
    dplyr::select(donor, all_of(pc_use)) %>%
    dplyr::filter(donor %in% peer_cond$donor)
  
  covar <- peer_cond %>%
    dplyr::left_join(pc_cond, by = "donor") %>%
    dplyr::select(donor, all_of(peer_use), all_of(pc_use))
  
  # Warn if NAs snuck in (shouldn't happen if PEER and PC files are complete)
  n_na <- sum(is.na(covar))
  if (n_na > 0) {
    cat(sprintf("  WARNING: %d NAs in covariate matrix for %s — dropping rows\n",
                n_na, condition_name))
    covar <- covar[complete.cases(covar), ]
  }
  
  mat <- cbind(intercept = 1, as.matrix(covar[, -1]))
  rownames(mat) <- covar$donor
  mat
}

# residualize_geno_for_condition: project out covariates from genotype dosages.
# Uses QR decomposition for numerical stability (better than solve() at small N).
# Result is genotype residuals that match the effective LD in the eQTL model.
residualize_geno_for_condition <- function(condition_name, snp_ids) {
  covar_mat    <- build_covar_matrix(condition_name)
  donors_avail <- intersect(rownames(covar_mat), colnames(geno_mat))
  
  cat(sprintf("  [%s] Residualizing LD using %d donors\n",
              condition_name, length(donors_avail)))
  
  G <- t(geno_mat[snp_ids, donors_avail, drop = FALSE])  # donors x SNPs
  C <- covar_mat[donors_avail, , drop = FALSE]            # donors x covariates
  
  Q       <- qr.Q(qr(C))
  G_resid <- G - Q %*% (t(Q) %*% G)   # residuals after projecting out C
  
  t(G_resid)   # back to SNPs x donors
}

# clean_LD: compute correlation matrix from residualized genotypes, enforce
# symmetry, and replace NAs with 0 (can occur for near-monomorphic SNPs).
clean_LD <- function(geno_resid) {
  LD <- cor(t(geno_resid), use = "pairwise.complete.obs")
  LD[is.na(LD)] <- 0
  (LD + t(LD)) / 2   # enforce exact symmetry (floating point can break this)
}

geno_resid_PBS  <- residualize_geno_for_condition("PBS",    shared_rsids)
geno_resid_cyto <- residualize_geno_for_condition(cytokine, shared_rsids)

LD_PBS  <- clean_LD(geno_resid_PBS)
LD_cyto <- clean_LD(geno_resid_cyto)

# Eigenvalue check — negative eigenvalues indicate non-positive-definite LD.
# Values around -1e-13 are floating point noise and are fine.
# Values below -1e-5 would indicate a real problem.
eigen_min_PBS  <- min(eigen(LD_PBS,  symmetric = TRUE, only.values = TRUE)$values)
eigen_min_cyto <- min(eigen(LD_cyto, symmetric = TRUE, only.values = TRUE)$values)
cat(sprintf("  LD min eigenvalue — PBS: %.3e | %s: %.3e\n",
            eigen_min_PBS, cytokine, eigen_min_cyto))

# Lambda check — median(chi-sq) / 0.4549. Under the null this should be ~1.
# Lambda > 2 after residualization suggests residual covariate mismatch.
# Lambda being high due to a strong real signal (e.g. ERAP2) is expected and fine.
z_PBS  <- df_PBS$beta_qtl_PBS   / sqrt(df_PBS$varbeta_qtl_PBS)
z_cyto <- df_cyto$beta_qtl_cyto / sqrt(df_cyto$varbeta_qtl_cyto)
cat(sprintf("  Lambda — PBS: %.3f | %s: %.3f\n",
            median(z_PBS^2) / 0.4549, cytokine, median(z_cyto^2) / 0.4549))

# =============================================================================
# 9. ASSEMBLE COLOC DATASET OBJECTS
# =============================================================================
# Using beta + varbeta + sdY + N format (preferred over p-value + MAF + N).
# N = actual donors in the LD matrix, NOT metadata count. These differ because
# some donors are dropped when their PEER factors are missing. Using the wrong N
# directly causes SuSiE's prior variance estimator to blow up.

dataset_PBS <- list(
  beta     = df_PBS$beta_qtl_PBS,
  varbeta  = df_PBS$varbeta_qtl_PBS,
  N        = ncol(geno_resid_PBS),    # actual donors after PEER filtering
  sdY      = sdY_eqtl_PBS,
  snp      = df_PBS$rsid,
  position = vt$POS,
  type     = "quant",
  LD       = LD_PBS
)

dataset_cyto <- list(
  beta     = df_cyto$beta_qtl_cyto,
  varbeta  = df_cyto$varbeta_qtl_cyto,
  N        = ncol(geno_resid_cyto),   # actual donors after PEER filtering
  sdY      = sdY_eqtl_cyto,
  snp      = df_cyto$rsid,
  position = vt$POS,
  type     = "quant",
  LD       = LD_cyto
)

cat(sprintf("  N for SuSiE — PBS: %d | %s: %d\n",
            dataset_PBS$N, cytokine, dataset_cyto$N))

# =============================================================================
# 10. SuSiE + COLOC (with full error recovery for bsub job arrays)
# =============================================================================
# The entire SuSiE + coloc block is wrapped in tryCatch so that a single bad
# gene doesn't kill the whole job array. On any error, a failure record is
# written and the job exits cleanly (status 0 so bsub marks it as succeeded).

tryCatch({
  
  # ---------------------------------------------------------------------------
  # 10a. SuSiE fine-mapping
  # ---------------------------------------------------------------------------
  # SuSiE parameters:
  #   L = 10: maximum number of independent signals per locus. Most eQTL loci
  #           have 1-3, but 10 is a safe upper bound. Results are stable to this
  #           choice as long as L >= true number of signals.
  #   estimate_residual_variance = FALSE: required because my LD matrix is from
  #           residualized genotypes, not the exact linear algebra inverse of the
  #           model. Fixing sigma^2 = 1 is appropriate since my betas and SEs
  #           already come from a correctly specified model.
  #   scaled_prior_variance = 0.2: prior SD on effect sizes ~ 0.45 * sdY.
  #           Slightly more permissive than default (0.1) to handle strong eQTLs
  #           like ERAP2 where beta can be ~2 on VST scale.
  #   min_abs_corr = 0.5: purity threshold for credible sets. SNPs in a credible
  #           set must have min pairwise |r| >= 0.5. Impure sets are unreliable
  #           and are dropped by coloc.susie. Can lower to 0.3 if many genes
  #           fail to produce credible sets, at cost of reliability.
  #
  # TWO-STAGE APPROACH:
  #   Stage 1: estimate prior variance (preferred — adapts to signal strength)
  #   Stage 2: fix prior variance (fallback for very strong signals like ERAP2
  #            where estimation becomes numerically unstable at small N)
  
  susie_L            <- 10
  susie_min_abs_corr <- 0.5
  
  run_susie_safe <- function(dataset, label) {
    
    # Stage 1
    result <- tryCatch({
      susie_rss(
        bhat                       = dataset$beta,
        shat                       = sqrt(dataset$varbeta),
        R                          = dataset$LD,
        n                          = dataset$N,
        L                          = susie_L,
        estimate_residual_variance = FALSE,
        scaled_prior_variance      = 0.2,
        min_abs_corr               = susie_min_abs_corr,
        verbose                    = FALSE
      )
    }, error = function(e) NULL)
    
    # Stage 2 — only reached if Stage 1 failed
    if (is.null(result)) {
      cat(sprintf("  [%s] Stage 1 failed, trying fixed prior variance...\n", label))
      result <- tryCatch({
        susie_rss(
          bhat                       = dataset$beta,
          shat                       = sqrt(dataset$varbeta),
          R                          = dataset$LD,
          n                          = dataset$N,
          L                          = susie_L,
          estimate_residual_variance = FALSE,
          estimate_prior_variance    = FALSE,  # fix at scaled_prior_variance
          scaled_prior_variance      = 0.2,
          min_abs_corr               = susie_min_abs_corr,
          verbose                    = FALSE
        )
      }, error = function(e) {
        cat(sprintf("  SuSiE failed for %s: %s\n", label, conditionMessage(e)))
        NULL
      })
    }
    result
  }
  
  cat("  Running SuSiE on PBS...\n")
  susie_PBS  <- run_susie_safe(dataset_PBS,  "PBS")
  cat(sprintf("  Running SuSiE on %s...\n", cytokine))
  susie_cyto <- run_susie_safe(dataset_cyto, cytokine)
  
  # Report credible set counts and convergence
  report_susie <- function(s, label) {
    if (is.null(s)) {
      cat(sprintf("  %s: SuSiE failed — no credible sets\n", label))
      return(invisible(NULL))
    }
    n_cs <- if (is.null(s$sets$cs)) 0L else length(s$sets$cs)
    cat(sprintf("  %s: %d credible set(s) | converged: %s\n",
                label, n_cs, isTRUE(s$converged)))
  }
  report_susie(susie_PBS,  "PBS")
  report_susie(susie_cyto, cytokine)
  
  # ---------------------------------------------------------------------------
  # 10b. coloc.susie
  # ---------------------------------------------------------------------------
  # coloc.susie runs pairwise coloc between each PBS credible set and each
  # cytokine credible set. This is the key advantage over standard coloc.abf —
  # it handles multiple independent signals per locus correctly.
  #
  # Output: one row per credible set pair tested.
  # I report the pair with the highest PP.H4 as the "best" result, which is
  # conservative for H3 claims (I'm not cherry-picking the H3 pairs).
  
  result_summary  <- NULL
  result_cs_pairs <- NULL
  
  if (!is.null(susie_PBS) && !is.null(susie_cyto)) {
    dataset_PBS$susie.fit  <- susie_PBS
    dataset_cyto$susie.fit <- susie_cyto
    
    coloc_result <- tryCatch(
      coloc.susie(dataset_PBS, dataset_cyto),
      error = function(e) {
        cat(sprintf("  coloc.susie failed: %s\n", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(coloc_result)) {
      result_summary  <- coloc_result$summary
      result_cs_pairs <- coloc_result$results
      cat("  coloc.susie completed.\n")
      print(result_summary)
    }
  }
  
  # ---------------------------------------------------------------------------
  # 10c. coloc.abf (fallback + sanity check)
  # ---------------------------------------------------------------------------
  # coloc.abf does not use LD and assumes one causal variant per locus.
  # It is used as:
  #   (1) Fallback when SuSiE fails (weak signal or numerical issues)
  #   (2) Sanity check against coloc.susie — if they agree on H3 vs H4,
  #       the conclusion is robust. If they disagree sharply, flag the gene.
  #
  # coloc.abf does not need LD — remove it to avoid confusing the function.
  
  dataset_PBS_noLD          <- dataset_PBS
  dataset_PBS_noLD$LD       <- NULL
  dataset_PBS_noLD$susie.fit <- NULL
  dataset_cyto_noLD          <- dataset_cyto
  dataset_cyto_noLD$LD       <- NULL
  dataset_cyto_noLD$susie.fit <- NULL
  
  coloc_abf <- tryCatch(
    coloc.abf(dataset_PBS_noLD, dataset_cyto_noLD),
    error = function(e) {
      cat(sprintf("  coloc.abf failed: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(coloc_abf)) {
    cat("  coloc.abf summary:\n")
    print(coloc_abf$summary)
  }
  
  # ---------------------------------------------------------------------------
  # 10d. Compile output row
  # ---------------------------------------------------------------------------
  n_cs_PBS  <- if (!is.null(susie_PBS)  && !is.null(susie_PBS$sets$cs))
    length(susie_PBS$sets$cs)  else 0L
  n_cs_cyto <- if (!is.null(susie_cyto) && !is.null(susie_cyto$sets$cs))
    length(susie_cyto$sets$cs) else 0L
  
  # Best coloc.susie pair = row with highest PP.H4.
  # Reporting highest PP.H4 is conservative for H3 claims — I am not
  # selecting the pair that looks most like H3.
  if (!is.null(result_summary) && nrow(result_summary) > 0) {
    best_row    <- result_summary[which.max(result_summary$PP.H4.abf), ]
    PP_H0_susie <- best_row$PP.H0.abf
    PP_H1_susie <- best_row$PP.H1.abf
    PP_H2_susie <- best_row$PP.H2.abf
    PP_H3_susie <- best_row$PP.H3.abf
    PP_H4_susie <- best_row$PP.H4.abf
    hit1_susie  <- as.character(best_row$hit1)
    hit2_susie  <- as.character(best_row$hit2)
    n_cs_pairs  <- nrow(result_summary)
  } else {
    PP_H0_susie <- PP_H1_susie <- PP_H2_susie <- PP_H3_susie <- PP_H4_susie <- NA_real_
    hit1_susie  <- hit2_susie  <- NA_character_
    n_cs_pairs  <- 0L
  }
  
  if (!is.null(coloc_abf)) {
    s         <- coloc_abf$summary
    PP_H0_abf <- s["PP.H0.abf"]; PP_H1_abf <- s["PP.H1.abf"]
    PP_H2_abf <- s["PP.H2.abf"]; PP_H3_abf <- s["PP.H3.abf"]
    PP_H4_abf <- s["PP.H4.abf"]
  } else {
    PP_H0_abf <- PP_H1_abf <- PP_H2_abf <- PP_H3_abf <- PP_H4_abf <- NA_real_
  }
  
  susie_PBS_converged  <- isTRUE(susie_PBS$converged)
  susie_cyto_converged <- isTRUE(susie_cyto$converged)
  
  # Reliability flags — use these downstream to stratify results:
  #   flag_low_N: N < 60 in either condition — fine-mapping less reliable
  #   flag_low_snp_overlap: < 500 shared SNPs — possible ascertainment issue
  #   flag_susie_failed: SuSiE returned NULL — coloc.susie not run
  #   flag_no_credible_sets: SuSiE ran but found no CS — gene has weak signal
  #   flag_susie_not_converged: IBSS didn't converge — results less reliable
  output_row <- tibble::tibble(
    celltype                 = ct,
    gene                     = this_gene,
    cytokine                 = cytokine,
    N_PBS                    = ncol(geno_resid_PBS),
    N_cyto                   = ncol(geno_resid_cyto),
    n_snps_PBS               = nrow(df_eqtl_PBS),
    n_snps_cyto              = nrow(df_eqtl_cyto),
    n_snps_shared            = length(shared_rsids),
    n_cs_PBS                 = n_cs_PBS,
    n_cs_cyto                = n_cs_cyto,
    n_cs_pairs_tested        = n_cs_pairs,
    susie_PP_H0              = PP_H0_susie,
    susie_PP_H1              = PP_H1_susie,
    susie_PP_H2              = PP_H2_susie,
    susie_PP_H3              = PP_H3_susie,
    susie_PP_H4              = PP_H4_susie,
    susie_hit1               = hit1_susie,
    susie_hit2               = hit2_susie,
    abf_PP_H0                = PP_H0_abf,
    abf_PP_H1                = PP_H1_abf,
    abf_PP_H2                = PP_H2_abf,
    abf_PP_H3                = PP_H3_abf,
    abf_PP_H4                = PP_H4_abf,
    flag_low_N               = (ncol(geno_resid_PBS) < 60) | (ncol(geno_resid_cyto) < 60),
    flag_low_snp_overlap     = length(shared_rsids) < 500,
    flag_susie_failed        = is.null(susie_PBS) | is.null(susie_cyto),
    flag_no_credible_sets    = (n_cs_PBS == 0) | (n_cs_cyto == 0),
    flag_susie_not_converged = !susie_PBS_converged | !susie_cyto_converged,
    skip_reason              = NA_character_
  )
  
  readr::write_tsv(output_row,
                   file.path(out_dir, paste0(out_prefix, "_coloc_summary.tsv")))
  cat(sprintf("  Summary written to: %s\n",
              file.path(out_dir, paste0(out_prefix, "_coloc_summary.tsv"))))
  
  # Save full credible set pair table — useful for inspecting which CS pairs
  # are H3 vs H4 when a gene has multiple independent signals
  if (!is.null(result_summary) && nrow(result_summary) > 0) {
    result_summary %>%
      dplyr::mutate(celltype = ct, gene = this_gene, cytokine = cytokine,
                    .before = 1) %>%
      readr::write_tsv(file.path(out_dir, paste0(out_prefix, "_coloc_cs_pairs.tsv")))
  }
  
  # ---------------------------------------------------------------------------
  # 10e. Diagnostic plot (non-fatal — plotting failure does not abort the job)
  # ---------------------------------------------------------------------------
  tryCatch({
    
    # make_cs_membership: label each SNP with its credible set name or "none"
    # Used to color Manhattan plots by credible set membership
    make_cs_membership <- function(susie_fit, rsids) {
      cs_vec <- rep("none", length(rsids))
      if (is.null(susie_fit) || is.null(susie_fit$sets$cs)) return(cs_vec)
      for (cs_name in names(susie_fit$sets$cs)) {
        idx <- susie_fit$sets$cs[[cs_name]]
        cs_vec[idx] <- cs_name
      }
      cs_vec
    }
    
    # Extended color palette to cover up to L=10 credible sets
    cs_colors <- c(
      "none" = "grey80",
      "L1"   = "#E41A1C", "L2"  = "#377EB8", "L3"  = "#4DAF4A",
      "L4"   = "#FF7F00", "L5"  = "#984EA3", "L6"  = "#A65628",
      "L7"   = "#F781BF", "L8"  = "#999999", "L9"  = "#66C2A5",
      "L10"  = "#FC8D62"
    )
    
    plot_df <- tibble::tibble(
      rsid      = df_PBS$rsid,
      pos       = vt$POS,
      logp_PBS  = -log10(df_PBS$p_qtl_PBS),
      logp_cyto = -log10(df_cyto$p_qtl_cyto),
      beta_PBS  = df_PBS$beta_qtl_PBS,
      beta_cyto = df_cyto$beta_qtl_cyto,
      cs_PBS    = make_cs_membership(susie_PBS,  df_PBS$rsid),
      cs_cyto   = make_cs_membership(susie_cyto, df_cyto$rsid)
    )
    
    p_manhat_PBS <- ggplot(plot_df, aes(x = pos, y = logp_PBS, color = cs_PBS)) +
      geom_point(size = 0.8, alpha = 0.7) +
      scale_color_manual(values = cs_colors, name = "CS (PBS)") +
      labs(title = sprintf("%s | %s | PBS eQTL", ct, this_gene),
           x = NULL, y = expression(-log[10](p))) +
      theme_bw(base_size = 11) + theme(legend.position = "right")
    
    p_manhat_cyto <- ggplot(plot_df, aes(x = pos, y = logp_cyto, color = cs_cyto)) +
      geom_point(size = 0.8, alpha = 0.7) +
      scale_color_manual(values = cs_colors, name = sprintf("CS (%s)", cytokine)) +
      labs(title = sprintf("%s | %s | %s eQTL", ct, this_gene, cytokine),
           x = "Position", y = expression(-log[10](p))) +
      theme_bw(base_size = 11) + theme(legend.position = "right")
    
    p_scatter <- ggplot(plot_df, aes(x = beta_PBS, y = beta_cyto, size = logp_PBS)) +
      geom_point(alpha = 0.4, color = "steelblue") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = 0, color = "grey70") +
      geom_vline(xintercept = 0, color = "grey70") +
      scale_size_continuous(name = expression(-log[10](p[PBS])), range = c(0.3, 3)) +
      labs(title = "Effect size concordance",
           x = "Beta (PBS)", y = sprintf("Beta (%s)", cytokine)) +
      theme_bw(base_size = 11)
    
    # Posterior probability bar chart
    # Shows coloc.susie (best CS pair) and coloc.abf side by side
    if (!is.na(PP_H3_susie)) {
      pp_df <- tibble::tibble(
        hypothesis = factor(
          c("H0\n(no signal)", "H1\n(PBS only)", "H2\n(cyto only)",
            "H3\n(distinct)", "H4\n(shared)"),
          levels = c("H0\n(no signal)", "H1\n(PBS only)", "H2\n(cyto only)",
                     "H3\n(distinct)", "H4\n(shared)")
        ),
        PP_susie = c(PP_H0_susie, PP_H1_susie, PP_H2_susie, PP_H3_susie, PP_H4_susie),
        PP_abf   = c(PP_H0_abf,   PP_H1_abf,   PP_H2_abf,   PP_H3_abf,   PP_H4_abf)
      ) %>%
        tidyr::pivot_longer(cols = c(PP_susie, PP_abf),
                            names_to = "method", values_to = "PP") %>%
        dplyr::mutate(method = dplyr::recode(method,
                                             "PP_susie" = "coloc.susie (best CS pair)",
                                             "PP_abf"   = "coloc.abf (sanity check)"))
      
      p_posteriors <- ggplot(pp_df, aes(x = hypothesis, y = PP, fill = method)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = c("coloc.susie (best CS pair)" = "#2166AC",
                                     "coloc.abf (sanity check)"   = "#D1E5F0")) +
        scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
        labs(title = "Coloc posteriors", x = NULL,
             y = "Posterior probability", fill = NULL) +
        theme_bw(base_size = 11) + theme(legend.position = "bottom")
      
    } else {
      # SuSiE failed — show coloc.abf only
      pp_df <- tibble::tibble(
        hypothesis = factor(
          c("H0\n(no signal)", "H1\n(PBS only)", "H2\n(cyto only)",
            "H3\n(distinct)", "H4\n(shared)"),
          levels = c("H0\n(no signal)", "H1\n(PBS only)", "H2\n(cyto only)",
                     "H3\n(distinct)", "H4\n(shared)")
        ),
        PP = c(PP_H0_abf, PP_H1_abf, PP_H2_abf, PP_H3_abf, PP_H4_abf)
      )
      p_posteriors <- ggplot(pp_df, aes(x = hypothesis, y = PP)) +
        geom_col(fill = "#D1E5F0") +
        scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
        labs(title = "Coloc posteriors (coloc.abf only — SuSiE failed)",
             x = NULL, y = "Posterior probability") +
        theme_bw(base_size = 11)
    }
    
    combined_plot <- (p_manhat_PBS / p_manhat_cyto / (p_scatter | p_posteriors)) +
      patchwork::plot_annotation(
        title   = sprintf("%s | %s | PBS vs %s coloc", ct, this_gene, cytokine),
        caption = sprintf(
          "N_PBS=%d | N_%s=%d | shared SNPs=%d | CS_PBS=%d | CS_%s=%d | converged: PBS=%s %s=%s",
          ncol(geno_resid_PBS), cytokine, ncol(geno_resid_cyto),
          length(shared_rsids),
          n_cs_PBS, cytokine, n_cs_cyto,
          susie_PBS_converged, cytokine, susie_cyto_converged
        )
      )
    
    out_plot_file <- file.path(out_dir, paste0(out_prefix, "_coloc_diagnostic.pdf"))
    ggsave(out_plot_file, combined_plot, width = 10, height = 12)
    cat(sprintf("  Diagnostic plot written to: %s\n", out_plot_file))
    
  }, error = function(e) {
    # Plotting failure is non-fatal — TSV output already written above
    cat(sprintf("  WARNING: plotting failed (non-fatal): %s\n", conditionMessage(e)))
  })
  
  cat(sprintf("\n[DONE] %s | %s | PBS vs %s\n", ct, this_gene, cytokine))
  cat(sprintf(
    "  coloc.susie: H3=%.3f, H4=%.3f | coloc.abf: H3=%.3f, H4=%.3f\n",
    ifelse(is.na(PP_H3_susie), -1, PP_H3_susie),
    ifelse(is.na(PP_H4_susie), -1, PP_H4_susie),
    ifelse(is.na(PP_H3_abf),   -1, PP_H3_abf),
    ifelse(is.na(PP_H4_abf),   -1, PP_H4_abf)
  ))
  
}, error = function(e) {
  # =============================================================================
  # TOP-LEVEL ERROR HANDLER
  # =============================================================================
  # If anything above fails unexpectedly, write a failure record and exit cleanly.
  # exit status 0 means bsub marks the job as succeeded — this prevents a single
  # bad gene from failing the entire job array and blocking downstream aggregation.
  # Check skip_reason == "ERROR:..." in the aggregated results to find failures.
  cat(sprintf("\n[ERROR] %s | %s | %s: %s\n",
              ct, this_gene, cytokine, conditionMessage(e)))
  tibble::tibble(
    celltype    = ct,
    gene        = this_gene,
    cytokine    = cytokine,
    skip_reason = paste("ERROR:", conditionMessage(e))
  ) %>%
    readr::write_tsv(
      file.path(out_dir, paste0(out_prefix, "_coloc_summary.tsv"))
    )
  quit(save = "no", status = 0)
})
