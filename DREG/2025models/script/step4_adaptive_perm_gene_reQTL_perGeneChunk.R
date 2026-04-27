#!/usr/bin/env Rscript
# script for running adaptive permutation for reQTLs

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(future.apply)
  library(magrittr)
  library(purrr)
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

# Compute vector of correlations between r (length D, vector of residuals)
# and each column of G (genotype matrix, rows (D): donors, columns (M): SNPs)
# Returns numeric length M
cor_vec_mat <- function(r, G) {
  r_sc <- center_scale_vec(r)
  if (all(is.na(r_sc))) stop("Residual vector has zero variance or all NA.")
  # scale each SNP column (genotype) to mean0 sd1
  G_sc <- apply(G, 2, center_scale_vec)
  # apply() returns matrix with same dims only if M>1
  if (is.null(dim(G_sc))) G_sc <- matrix(G_sc, ncol = 1)
  
  # Pearson Correlation = mean(r_sc * g_sc) over donors with complete data, * means dot product
  # We'll do column-wise using matrix multiplication, but must handle NAs.
  # Simplest robust approach: restrict to donors complete for r and all SNPs is too strict.
  # Instead: handle per-SNP complete cases (still fast enough for ~1-5k SNPs). 
  # sum(ok) < 4 = not enough data to trust this SNP
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
# Adaptive permutation engine for one gene, now with optional B3.
# returns list with p_hat, T_obs, perms_run, exceedances, stage, (optional) T_perm
adaptive_perm_gene <- function(r, G,
                               B0 = 200, B1 = 1000, B2 = 10000,
                               B3 = NULL,                 # e.g. 250000 or 500000 to enable
                               seed = 1,
                               keep_T_perm = FALSE) {
  
  set.seed(seed)
  
  # Observed statistic
  z_obs <- cor_vec_mat(r, G)
  T_obs <- max(abs(z_obs), na.rm = TRUE)
  
  # Stage runner: adds permutations up to target B
  # - If keep_T_perm=FALSE: does NOT store permutation statistics (saves memory)
  # - If keep_T_perm=TRUE: stores T_perm so you can plot the null later
  run_stage <- function(B, b_prev = 0, T_prev = numeric(0), n_prev = 0L) {
    b <- b_prev
    D <- length(r)
    
    # how many perms already done?
    n_done <- if (keep_T_perm) length(T_prev) else n_prev
    n_to_add <- B - n_done
    if (n_to_add <= 0) {
      return(list(b = b, T_perm = if (keep_T_perm) T_prev else NULL, n_done = n_done))
    }
    
    if (keep_T_perm) {
      T_perm <- c(T_prev, rep(NA_real_, n_to_add))
      start_idx <- n_done + 1L
    } else {
      T_perm <- NULL
    }
    
    for (i in seq_len(n_to_add)) {
      rp <- sample(r, size = D, replace = FALSE)
      z_p <- cor_vec_mat(rp, G)
      t_p <- max(abs(z_p), na.rm = TRUE)
      
      if (keep_T_perm) {
        T_perm[start_idx + i - 1L] <- t_p
      }
      
      if (is.finite(t_p) && is.finite(T_obs) && t_p >= T_obs) {
        b <- b + 1L
      }
    }
    
    list(b = b, T_perm = T_perm, n_done = B)
  }
  
  # Helper to standardize return payload
  make_out <- function(perms, b, p_hat, stage, T_perm = NULL, next_B = NA_integer_, timing = NULL) {
    out <- list(
      T_obs  = T_obs,
      perms  = perms,
      b      = b,
      p_hat  = p_hat,
      stage  = stage,
      next_B = next_B,
      timing = timing
    )
    if (keep_T_perm) out$T_perm <- T_perm
    out
  }
  
  # ------------------ Stage 0 ------------------
  t0_start <- Sys.time()
  st0 <- run_stage(B0, b_prev = 0L, T_prev = numeric(0), n_prev = 0L)
  t0_sec <- as.numeric(difftime(Sys.time(), t0_start, units = "secs"))
  
  b0 <- st0$b
  p0 <- (b0 + 1) / (B0 + 1)
  timing0 <- list(stage = "B0", seconds = t0_sec, sec_per_perm = t0_sec / B0)
  
  # stop early if any exceedance
  if (b0 > 0) {
    return(make_out(B0, b0, p0, "stopped@B0", T_perm = st0$T_perm, timing = timing0))
  }
  
  # ------------------ Stage 1 ------------------
  t1_start <- Sys.time()
  st1 <- run_stage(B1,
                   b_prev = b0,
                   T_prev = if (keep_T_perm) st0$T_perm else numeric(0),
                   n_prev = st0$n_done)
  t1_sec <- as.numeric(difftime(Sys.time(), t1_start, units = "secs"))
  
  b1 <- st1$b
  p1 <- (b1 + 1) / (B1 + 1)
  timing1 <- list(stage = "B1", seconds = t1_sec, sec_per_perm = t1_sec / (B1 - B0))
  
  if (b1 > 0) {
    return(make_out(B1, b1, p1, "stopped@B1", T_perm = st1$T_perm,
                    timing = list(B0 = timing0, B1 = timing1)))
  }
  
  # ------------------ Stage 2 ------------------
  t2_start <- Sys.time()
  st2 <- run_stage(B2,
                   b_prev = b1,
                   T_prev = if (keep_T_perm) st1$T_perm else numeric(0),
                   n_prev = st1$n_done)
  t2_sec <- as.numeric(difftime(Sys.time(), t2_start, units = "secs"))
  
  b2 <- st2$b
  p2 <- (b2 + 1) / (B2 + 1)
  timing2 <- list(stage = "B2", seconds = t2_sec, sec_per_perm = t2_sec / (B2 - B1))
  
  # If B3 is not requested (or invalid), stop here but flag b2==0 genes
  if (is.null(B3) || is.na(B3) || B3 <= B2) {
    stage_label <- if (b2 == 0) "needToRunB3" else "stopped@B2"
    next_B <- if (b2 == 0 && !(is.null(B3) || is.na(B3))) B3 else NA_integer_
    
    return(make_out(B2, b2, p2, stage_label, T_perm = st2$T_perm, next_B = next_B,
                    timing = list(B0 = timing0, B1 = timing1, B2 = timing2)))
  }
  
  # If any exceedance at B2, stop
  if (b2 > 0) {
    return(make_out(B2, b2, p2, "stopped@B2", T_perm = st2$T_perm,
                    timing = list(B0 = timing0, B1 = timing1, B2 = timing2)))
  }
  
  # ------------------ Stage 3 ------------------
  t3_start <- Sys.time()
  st3 <- run_stage(B3,
                   b_prev = b2,
                   T_prev = if (keep_T_perm) st2$T_perm else numeric(0),
                   n_prev = st2$n_done)
  t3_sec <- as.numeric(difftime(Sys.time(), t3_start, units = "secs"))
  
  b3 <- st3$b
  p3 <- (b3 + 1) / (B3 + 1)
  timing3 <- list(stage = "B3", seconds = t3_sec, sec_per_perm = t3_sec / (B3 - B2))
  
  return(make_out(B3, b3, p3, "stopped@B3", T_perm = st3$T_perm,
                  timing = list(B0 = timing0, B1 = timing1, B2 = timing2, B3 = timing3)))
}




#### Main ####
# ------------------------------#
# Main: parse CLI args ####
# ------------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("
Usage:
  Rscript adaptive_perm_gene_eQTL.R <ct> <condition> <QTLtype> <chunk_id> [out_prefix]

Example:
  Rscript adaptive_perm_gene_eQTL.R MEL TNF eQTL 000 \\
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
condition  <- args[[2]]
QTLtype    <- args[[3]]
chunk_id       <- args[[4]] # 000 to 120
this_condition <- condition

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
pair_file <- paste0(dir,"/permutation/chunk/pair_chunk_",chunk_id,".txt")
geno_file  <- paste0(dir,"/permutation/chunk/genotype_pair_chunk_",chunk_id,".txt")
modelstats_file <- paste0(dir,"/permutation/model_stats/model_stats_",chunk_id,".txt")

vst_file   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
peer_file  <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_",ct,"_PBS-IFNG-IFNB-TNF.tsv")
pc_file    <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"
out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")

n_workers <- 2
base_seed <- 42

# toy example for debugging
ct         <- "MEL"
condition  <- "IFNB"
QTLtype    <- "reQTL"
chunk_id <- "034"
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
pair_file <- paste0(dir,"/permutation/chunk/pair_chunk_",chunk_id,".txt")
geno_file  <- paste0(dir,"/permutation/chunk/genotype_pair_chunk_",chunk_id,".txt")
modelstats_file <- paste0(dir,"/permutation/model_stats/model_stats_",chunk_id,".txt")
this_condition <- condition
vst_file   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
peer_file  <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_",ct,"_PBS-IFNG-IFNB-TNF.tsv")
pc_file    <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"
out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")

n_workers <- 2
base_seed <- 42
message("=== Adaptive permutation for ", ct, " ", condition, " ", QTLtype, ", chunk: ", chunk_id, " ===")

#### ----Load data sets -------- ####
#### --- Load chunk files --------- #
# load SNP:gene pairs
pairs <- fread(pair_file, header=F) %>% 
  magrittr::set_colnames(c("SNP","gene"))
pairs$key <- paste0(pairs$gene,"_",pairs$SNP)

# load genotype of all SNPs in chunk
genotype_all <- fread(geno_file)

# load modeling stats for the QTL of interest of chunk
modelstats_all <- fread(modelstats_file)
fixed_col <- c("snp","gene","dist")
stats_col <- setdiff(colnames(modelstats_all), fixed_col)
pattern <- paste0("^",QTLtype,".*",condition,"$")
selected_stats_col <- grep(pattern, stats_col, value=TRUE)
modelstats_subset <- modelstats_all[ , c(fixed_col, selected_stats_col), with = FALSE] %>%
  dplyr::mutate(key = paste0(gene,"_",snp)) %>%
  dplyr::filter(key %in% pairs$key)

#### --- Load general data files ------------- #
# Metadata (for reQTL, select both PBS and stimulant conditions)
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(celltype == ct, condition %in% c("PBS",!!condition))

# Covariates
nPEER = 10
nGPC = 2

# PEERs (sample-level)
peer <- fread(peer_file) %>% 
  dplyr::select(-c(PEER1, PEER2, PEER3, PEER4))
colnames(peer)[1] <- "sample"
peer <- peer[ , 1:(nPEER+1)]
peer <- left_join(peer, meta[,c("sample","donor","condition")]) %>%
  arrange(donor)
paired_donors <- table(peer[,c(donor,condition)]) %>% as.data.frame %>% dplyr::filter(Freq==2) %>% pull(Var1) %>% as.character()

peer_cols <- setdiff(colnames(peer), c("sample", "donor", "condition"))
peer.pbs <- peer %>%
  dplyr::filter(donor %in% paired_donors, condition == "PBS") %>%
  dplyr::select(donor, dplyr::all_of(peer_cols)) %>%
  dplyr::arrange(donor)

peer.stim <- peer %>%
  dplyr::filter(donor %in% paired_donors, condition == this_condition) %>%
  dplyr::select(donor, dplyr::all_of(peer_cols)) %>%
  dplyr::arrange(donor)
peer.pbs_tbl  <- dplyr::as_tibble(peer.pbs)
peer.stim_tbl <- dplyr::as_tibble(peer.stim)

# safety: confirm donors align 1:1
stopifnot(identical(peer.pbs$donor, peer.stim$donor))
peer.delta <- peer.pbs_tbl
peer.delta[peer_cols] <- peer.stim_tbl[peer_cols] - peer.pbs_tbl[peer_cols]


# genotype PCs (donor-level)
#pcs <- fread(pc_file) %>% dplyr::select(-c(FID, donor_num))
#pcs <- pcs[ , 1:(nGPC+1)]

# Read VST wide, keep gene and samples for this celltype+condition
vst <- read_tsv(vst_file, show_col_types = FALSE) %>% 
  dplyr::filter(gene %in% unique(pairs$gene)) %>%
  column_to_rownames("gene") %>%
  t() %>% as.data.frame %>%
  rownames_to_column("sample") %>%
  left_join(meta[,c("sample","donor","condition")]) %>%
  dplyr::filter(donor %in% paired_donors) %>%
  arrange(donor)

vst.pbs <- vst[which(vst$condition=="PBS"),] %>%
  dplyr::select(-c("sample","condition")) %>%
  tibble::remove_rownames() %>%
  dplyr::as_tibble()

vst.stim <- vst[which(vst$condition==this_condition),] %>%
  dplyr::select(-c("sample","condition")) %>%
  tibble::remove_rownames() %>%
  dplyr::as_tibble()

vst_cols <- setdiff(colnames(vst.pbs),"donor")
stopifnot(identical(vst.pbs$donor, vst.stim$donor))
vst.delta <- vst.pbs
vst.delta[vst_cols] <- vst.stim[vst_cols] - vst.pbs[vst_cols]
vst.delta <- vst.delta %>% dplyr::select(donor, dplyr::everything())

#### --- Process for every gene --------- ####
run_one_gene <- function(this_gene,
                         pairs, genotype_all, vst.delta, meta, peer.delta,
                         ct, condition, QTLtype,
                         B0=200, B1=1000, B2=10000, B3=NULL,
                         seed=1, keep_T_perm=FALSE) {
  
  message(this_gene)
  # ---- SNP list for this gene ####
  snps <- pairs %>% dplyr::filter(gene == this_gene) %>% pull(SNP) %>% unique()
  
  # ---- Step 2: build donor-level expression for this gene in this condition ####
  # Subset VST to this gene only
  vst_gene <- vst.delta %>% dplyr::select(all_of(c("donor",this_gene))) %>%
    set_colnames(c("donor","delta_expr"))
  
  # Make long
  expr_tbl <- vst_gene %>%
    inner_join(peer.delta, by="donor")
  donors <- expr_tbl$donor
  
  message("Paired Donors with ", condition, " expression: ", nrow(expr_tbl))
  
  # ---- Step 3: residualize expression w.r.t covariates (NO genotype) ####
  # The pragmatic approach for reQTL minP across SNPs is:
  # Build null formula: delta expr ~ delta PEERs
  # Test: interaction term via permuting donor genotypes within the paired design
  peer_use <- setdiff(colnames(peer.delta), "donor")
  covars <- peer_use
  form0 <- as.formula(paste0("delta_expr ~ ", paste(covars, collapse=" + ")))
  fit0 <- lm(form0, data = expr_tbl)
  r <- resid(fit0)
  names(r) <- expr_tbl$donor
  
  # safety checks
  stopifnot(length(r) == length(unique(donors)))
  stopifnot(!any(is.na(r)))
  # Guard: zero variance residuals
  if (!is.finite(stats::sd(r)) || stats::sd(r) == 0) {
    return(tibble(celltype=ct, condition=condition, QTLtype=QTLtype, gene=this_gene,
                  n_snps=length(snps), n_donors=length(donors),
                  T_obs=NA_real_, perms=NA_integer_, exceed=NA_integer_,
                  p_emp=NA_real_, stage="resid_zero_variance", next_B=NA_integer_))
  }
  
  if (length(donors) < 5) {
    return(tibble(
      celltype=ct, condition=condition, QTLtype=QTLtype, gene=this_gene,
      n_snps=length(snps), n_donors=length(donors),
      T_obs=NA_real_, perms=NA_integer_, exceed=NA_integer_,
      p_emp=NA_real_, stage="too_few_donors", next_B=NA_integer_
    ))
  }
  
  # ---- Step 4: genotype matrix G donors x SNPs ----####
  fixed_cols <- c("CHROM","POS","ID","REF","ALT")
  geno_sample_cols <- setdiff(colnames(genotype_all), fixed_cols)
  
  # Keep only genotype donors that appear in expression data
  donors <- names(r)
  donors_use <- intersect(donors, geno_sample_cols)
  
  if (length(donors_use) < 5) {
    return(tibble(
      celltype=ct, condition=condition, QTLtype=QTLtype, gene=this_gene,
      n_snps=length(snps), n_donors=length(donors),
      T_obs=NA_real_, perms=NA_integer_, exceed=NA_integer_,
      p_emp=NA_real_, stage="no_geno_overlap", next_B=NA_integer_
    ))
  }
  
  # Reorder donors_use to match donors (i.e., residual order)
  donors_use <- donors[donors %in% donors_use]
  
  # select SNPs in the gene's cis-window
  geno_sub <- genotype_all %>%
    dplyr::filter(ID %in% snps) %>%
    dplyr::select(ID, all_of(donors_use))
  
  # Ensure SNP order matches 'snps'
  geno_sub <- geno_sub %>%
    dplyr::mutate(.ord = match(ID, snps)) %>%
    dplyr::arrange(.ord) %>%
    dplyr::select(-.ord)
  
  # Build G: donors x SNPs
  # align SNP order + donor order
  # ensure donor column order matches donors
  # geno_sub columns after ID are donor_cols_present; reorder to donors
  G_snp_by_donor <- as.matrix(dplyr::select(geno_sub, -ID))
  storage.mode(G_snp_by_donor) <- "numeric"
  
  # Transpose to donors x SNPs
  G <- t(G_snp_by_donor)
  
  # Also align r to donors_use (critical)
  r_use <- r[donors_use]
  stopifnot(identical(rownames(G), donors_use) || TRUE)  # matrix typically has no rownames
  
  # ---- Step 5: adaptive permutation ----####
  # make seed gene-specific so workers are reproducible but not identical
  gene_seed <- seed + (sum(utf8ToInt(this_gene)) %% 1000000)
  
  res <- adaptive_perm_gene(
    r = r_use, G = G,
    B0 = B0, B1 = B1, B2 = B2,
    B3 = B3,
    seed = gene_seed,
    keep_T_perm = keep_T_perm
  )
  
  # ---- Step 6: summary row ----
  tibble(
    celltype   = ct,
    condition  = condition,
    QTLtype    = QTLtype,
    gene       = this_gene,
    n_snps     = ncol(G),
    n_donors   = length(r_use),
    T_obs      = res$T_obs,
    perms      = res$perms,
    exceed     = res$b,
    p_emp      = res$p_hat,
    stage      = res$stage,
    next_B     = res$next_B
  )
}


# ----------- MAIN ---------------------------- ####
# Main: Run adaptive permutation on all genes in chunk across 5 workers
# ------------------------------#
start_time <- Sys.time()
out_file <- file.path(paste0(dir, "/permutation/result/", condition,"/",QTLtype,"/",out_prefix, ".perm_summary.tsv"))

future::plan(future::multicore, workers = n_workers)

gene_lst <- unique(pairs$gene)

# Optional: keep a deterministic seed stream across workers
# (recommended when you care about reproducibility)
results_list <- future.apply::future_lapply(
  gene_lst,
  function(g) {
    run_one_gene(
      this_gene = g,
      pairs = pairs,
      genotype_all = genotype_all,
      vst.delta = vst.delta,
      meta = meta,
      peer.delta = peer.delta,
      ct = ct,
      condition = condition,
      QTLtype = QTLtype,
      B0 = 200, B1 = 1000, B2 = 10000, B3 = NULL,
      seed = base_seed,
      keep_T_perm = FALSE
    )
  },
  future.seed = TRUE
)

summary_out <- dplyr::bind_rows(results_list)

end_time <- Sys.time()

cat("Start:", format(start_time), "\n")
cat("End:  ", format(end_time), "\n")
cat("Elapsed:", as.numeric(difftime(end_time, start_time, units = "secs")), "seconds\n")


# Write ONCE (safe)
readr::write_tsv(summary_out, out_file)
message("Wrote: ", out_file)
