#!/usr/bin/env Rscript
# Mixed-effects eQTL / reQTL per SNP-gene pair (sample-level PEERs)
# Inputs (paths are CLI-args or defaults below):
#   1) chunk_file: TSV with columns: snp   gene
#   2) output_file: TSV to write results
# Environment files (edit paths or pass via args 3..6):
#   - VST.txt (wide): gene_id | <sample...>
#   - sample_metadata.txt: sample, donor, condition
#   - PEERs.txt: sample, PEER1..PEERk  (SAMPLE-LEVEL)
#   - genotype_pcs.txt: donor, PC1, PC2, ... (DONOR-LEVEL)
#   - GENOTYPE: TSV with columns: CHR, POS, START, END, ID, REF, ALT, <sample1>, <sample2>, ...
#               GTs in 0,1,2 format
#
# Model (one fit per pair):
#   VST ~ geno * condition + PEERs + (1 | donor)
#   eQTL per condition via emtrends; reQTL (cond vs PBS) via contrasts
#

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
  stop("Usage: Rscript run_QTL_chunk.R <CELLTYPE> <nPEER> <nGPC> <chunk_id>")
}

ct <- args[1] # MEL, KRT, FRB
nPEER <- as.integer(args[2]) # number of PEER factors to include
nGPC <- as.integer(args[3]) # number of genotype PCs to include
chunk_id  <- args[4]

dir="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL"
chunk_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
GENO_FILE <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv") # columns: CHR POS ID REF ALT sample1 sample2 ...
output_file <- paste0(dir,"/results/result_",chunk_id,".tsv")


# toy example for debugging
ct <- "MEL"
nPEER=8
nGPC=2
chunk_id="000"
dir="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL"
chunk_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
#GENO_FILE <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv") # columns: CHR POS ID REF ALT sample1 sample2 ...
GENO_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots/genotype_debug.tsv"
output_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/results/test_result_000.tsv"

## ------------------ Paths ------------------
# Expression + metadata
VST_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
PEER_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_MEL_PBS-IFNG-IFNB-TNF.tsv"             # columns: sample, PEER1, PEER2, ...
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"      # columns: donor, PC1, PC2, ...

message("====== QTL chunk start ======")
message("Chunk:  ", chunk_file)
message("Output: ", output_file)

## ---------- 0) Settings ----------
# Parallel workers (10 workers is pre-determined through benchmarking)
plan(multicore, workers = 10)
set.seed(1)

conds_all <- c("PBS","IFNG","IFNB","TNF")  # allowed condition labels

## ---------- 1) Load once ---------- ####
# Load pairs
message("[1/7] Loading pairs from chunk ...")
pairs <- fread(chunk_file, header = TRUE)
colnames(pairs)[4] <- "gene"
colnames(pairs)[11] <- "snp"
# true negative but identified as false positive in TNF reQTL (p=0, beta=-235)
pairs <- rbind(pairs, as.list(rep(NA, ncol(pairs))), fill = TRUE)
gene="CALD1"; snp="rs759854"
pairs[nrow(pairs),11] <- snp
pairs[nrow(pairs),4] <- gene
# true positive and also detected positive in IFNG reQTL 
pairs <- rbind(pairs, as.list(rep(NA, ncol(pairs))), fill = TRUE)
gene="GBP3"; snp="rs72726528"
pairs[nrow(pairs),11] <- snp
pairs[nrow(pairs),4] <- gene
npairs <- nrow(pairs)
message("Pairs in chunk: ", npairs)

message("[2/7] Loading VST ...")
expr_wide <- fread(VST_FILE) %>%
  dplyr::filter(gene %in% pairs$gene) %>%
  dplyr::select(c("gene", contains(ct)))

message("[3/7] Loading sample metadata ...")
meta <- fread(META_FILE) %>% dplyr::filter(celltype==ct)
stopifnot(all(c("sample","donor","condition") %in% names(meta)))

# Keep VST columns that have metadata (for this celltype)
expr_samples <- intersect(colnames(expr_wide)[-1], meta$sample)
expr_wide <- as.data.table(expr_wide)[, c("gene_id", expr_samples), with = FALSE]
meta <- meta[sample %in% expr_samples]
meta$condition <- factor(as.character(meta$condition), levels = conds_all)

# Melt VST to long (gene_id, sample, VST) and merge metadata
message("[4/7] Reshaping VST and merging metadata ...")
expr_long <- melt(as.data.table(expr_wide), id.vars = "gene_id",
                  variable.name = "sample", value.name = "VST")
expr_long <- merge(expr_long, meta, by = "sample", all.x = TRUE, all.y = FALSE)
rm(expr_wide); invisible(gc())

# PEER factors at sample level
message("[5/7] Loading PEER (sample-level) ...")
peers <- fread(PEER_FILE)
colnames(peers)[1] <- "sample"
stopifnot("sample" %in% names(peers))
peer_cols <- setdiff(colnames(peers), "sample")
stopifnot(length(peer_cols) > 0)

# Merge PEERs onto sample rows directly (sample-level PEER merge)
expr_long <- merge(expr_long, peers, by = "sample", all.x = TRUE)
# Drop rows without PEERs if any remain
expr_long <- expr_long[complete.cases(expr_long[, ..peer_cols]), ]

# Genotype PCs (per donor)
message("[6/7] Loading genotype PCs ...")
if (file.exists(GENO_PCS_FILE)) {
  gpcs <- fread(GENO_PCS_FILE)
  stopifnot("donor" %in% names(gpcs))
  pc_cols <- setdiff(colnames(gpcs), "donor")
  expr_long <- merge(expr_long, gpcs, by = "donor", all.x = TRUE)
} else {
  pc_cols <- character(0)
  warning("GENO_PCS_FILE not found; proceeding without genotype PCs.")
}

# Load genotype table
message("[7/7] Loading GENOTYPE table ...")
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")

# Identify genotype sample columns
old_sample_cols <- setdiff(colnames(geno_dt), fixed_cols)

# Create NEW cleaned column names (donor IDs)
new_sample_cols <- old_sample_cols
new_sample_cols <- gsub("skineQTL-", "", new_sample_cols)
new_sample_cols <- gsub("^F0", "F", new_sample_cols)
setnames(geno_dt, old_sample_cols, new_sample_cols)

# Now work with the renamed columns from here forward:
geno_sample_cols <- new_sample_cols
if (length(geno_sample_cols) == 0) stop("No sample GT columns found in GENOTYPE.")

# Keep only genotype donors that appear in expression data
geno_keep_samples <- intersect(geno_sample_cols, unique(expr_long$donor))
if (length(geno_keep_samples) == 0) stop("No overlapping donor IDs between GENOTYPE and expression samples.")
geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples), with = FALSE]

## ---------- 2) Prep fast genotype access ---------- ####
# Key by SNP ID for O(1) row lookup
data.table::setkey(geno_dt, ID)

# Keep vector of donor columns present in genotype table
geno_donors <- setdiff(colnames(geno_dt), fixed_cols)  # donors in genotype table

# Small helper to fetch a named donor->geno vector (numeric 0/1/2) for one SNP
fetch_geno_for_snp <- function(snp_id, donor_vec) {
  row <- geno_dt[list(snp_id)]
  if (nrow(row) == 0L) return(rep(NA_real_, length(donor_vec)))
  g <- as.numeric(row[, ..geno_donors])
  names(g) <- geno_donors
  g[donor_vec]
}

## ---------- 3) Build model formula ---------- ####
peer_cols <- setdiff(colnames(peers), c("sample","PEER1","PEER2","PEER3","PEER4"))
pc_cols <- setdiff(colnames(gpcs), c("FID","donor","donor_num"))
# subset to number of factors to include
peer_use <- peer_cols[seq_len(nPEER)]
pc_use <- pc_cols[seq_len(nGPC)]
cov_terms <- c(peer_use, pc_use)
cov_part <- if (length(cov_terms) > 0) paste("+", paste(cov_terms, collapse = " + ")) else ""
# final formula string (random intercept by donor)
fml_str <- paste0("VST ~ geno * condition ", cov_part, " + (1 | donor)")
message("Model formula: ", fml_str)

## ---------- 4) Per-pair fit function ---------- ####
conds_all <- c("PBS","IFNG","IFNB","TNF")
# true negative but identified as false positive in TNF reQTL (p=0, beta=-235)
gene="CALD1"; snp="rs759854"
# true positive and also detected positive in IFNG reQTL 
gene="GBP3"; snp="rs72726528"

# ---------- QC CHECKS ----------

# (A) genotype support check — returns TRUE if some genotype has < 3 donors
qc_genotype_support <- function(df) {
  geno_counts <- df %>%
    group_by(geno) %>%
    summarize(n_donors = n_distinct(donor), .groups = "drop")
  min(geno_counts$n_donors) < 3
}

# (B) leave-one-donor-out sensitivity — returns SD of delta_IFNG across donors
qc_lodo_sensitivity <- function(df, fml, orig_IFNG_delta) {
  donors <- unique(df$donor)
  vals <- c()
  
  for (d in donors) {
    df_lodo <- df[df$donor != d, ]
    fit_lodo <- tryCatch(
      lmer(fml, data=df_lodo, REML=FALSE,
           control=lmerControl(optimizer="bobyqa", calc.derivs=FALSE)),
      error=function(e) NULL
    )
    if (is.null(fit_lodo)) next
    
    tr_lodo <- tryCatch(emmeans::emtrends(fit_lodo, ~ condition, var="geno"),
                        error=function(e) NULL)
    if (is.null(tr_lodo)) next
    
    re_lodo <- tryCatch(emmeans::contrast(tr_lodo, "trt.vs.ctrl", ref="PBS"),
                        error=function(e) NULL)
    if (is.null(re_lodo)) next
    
    tdf <- as.data.frame(summary(re_lodo))
    tgt <- subset(tdf, contrast == "IFNG - PBS")
    if (nrow(tgt) == 1) vals <- c(vals, tgt$estimate)
  }
  
  if (length(vals) < 2) return(NA_real_)
  sd(vals, na.rm = TRUE)
}

# (C) donor influence check — return max Cook's distance + donor name
qc_leverage <- function(fit) {
  infl <- tryCatch(influence.ME::influence(fit, group="donor"),
                   error=function(e) NULL)
  if (is.null(infl)) return(list(max_cooks=NA_real_, high_leverage=NA, top_donor=NA))
  
  cd <- influence.ME::cooks.distance(infl, sort = FALSE)
  max_c <- max(cd, na.rm=TRUE)
  list(
    max_cooks = max_c,
    high_leverage = max_c > 1,
    top_donor = names(which.max(cd))
  )
}

fit_one <- function(snp, gene) {
  cols_needed <- c("sample","donor","condition","VST", peer_cols, pc_cols)
  df <- tryCatch(expr_long[gene_id == gene, ..cols_needed], error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  
  df$condition <- factor(as.character(df$condition))
  present_conds <- intersect(conds_all, unique(df$condition))
  df <- df[condition %in% present_conds]
  if (nrow(df) < 20) return(NULL)
  
  if ("PBS" %in% present_conds)
    df$condition <- stats::relevel(df$condition, ref="PBS")
  
  df$geno <- fetch_geno_for_snp(snp, df$donor)
  if (all(is.na(df$geno))) return(NULL)
  
  ug <- na.omit(unique(df$geno))
  if (length(ug) < 2) return(NULL)  # need ≥2 genotype groups
  
  mac <- min(table(df$geno[!is.na(df$geno)]))
  if (is.infinite(mac) || mac < 3) return(NULL)
  
  # ---------- fit model ----------
  fml <- as.formula(fml_str)
  fit <- tryCatch(
    lmer(fml, data=df, REML=FALSE,
         control=lmerControl(optimizer="bobyqa", calc.derivs=FALSE)),
    error=function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  # ---------- eQTL slopes ----------
  tr <- tryCatch(emmeans::emtrends(fit, ~ condition, var="geno"),
                 error=function(e) NULL)
  if (is.null(tr)) return(NULL)
  tr_summ <- as.data.frame(summary(tr))
  
  # ---------- reQTL deltas ----------
  re_tbl <- NULL
  if ("PBS" %in% tr_summ$condition && length(unique(tr_summ$condition)) > 1) {
    re_con <- tryCatch(emmeans::contrast(tr, "trt.vs.ctrl", ref="PBS"),
                       error=function(e) NULL)
    if (!is.null(re_con)) {
      re_tbl <- as.data.frame(summary(re_con))
      re_tbl$cond <- sub(" - PBS$", "", re_tbl$contrast)
    }
  }
  
  # ---------- output table ----------
  out <- data.frame(
    snp = snp, gene = gene,
    eQTL_beta_PBS = NA_real_, eQTL_p_PBS = NA_real_,
    eQTL_beta_IFNG = NA_real_, eQTL_p_IFNG = NA_real_,
    eQTL_beta_IFNB = NA_real_, eQTL_p_IFNB = NA_real_,
    eQTL_beta_TNF = NA_real_,  eQTL_p_TNF  = NA_real_,
    reQTL_dbeta_IFNG = NA_real_, reQTL_p_IFNG = NA_real_,
    reQTL_dbeta_IFNB = NA_real_, reQTL_p_IFNB = NA_real_,
    reQTL_dbeta_TNF  = NA_real_, reQTL_p_TNF  = NA_real_
  )
  
  for (cond in intersect(conds_all, tr_summ$condition)) {
    r <- tr_summ[tr_summ$condition == cond, ]
    out[[paste0("eQTL_beta_", cond)]] <- r$geno.trend
    out[[paste0("eQTL_p_", cond)]]    <- r$p.value
  }
  
  orig_delta <- if (!is.null(re_tbl))
    re_tbl$estimate[re_tbl$cond == "IFNG"] else NA_real_
  
  if (!is.null(re_tbl)) {
    for (cond in intersect(c("IFNG","IFNB","TNF"), re_tbl$cond)) {
      r <- re_tbl[re_tbl$cond == cond, ]
      out[[paste0("reQTL_dbeta_", cond)]] <- r$estimate
      out[[paste0("reQTL_p_", cond)]]     <- r$p.value
    }
  }
  
  # ---------- QC checks ----------
  out$low_genotype_support <- qc_genotype_support(df)
  
  out$lodo_sd <- qc_lodo_sensitivity(df, fml, orig_delta)
  out$lodo_unstable <- !is.na(out$lodo_sd) &&
    out$lodo_sd > abs(orig_delta) * 0.25
  
  lev <- qc_leverage(fit)
  out$max_cooks <- lev$max_cooks
  out$high_leverage <- lev$high_leverage
  out$top_leverage_donor <- lev$top_donor
  
  out
}

  # ---------- output table ----------
  out <- data.frame(
    snp = snp, gene = gene,
    eQTL_beta_PBS = NA_real_, eQTL_p_PBS = NA_real_,
    eQTL_beta_IFNG = NA_real_, eQTL_p_IFNG = NA_real_,
    eQTL_beta_IFNB = NA_real_, eQTL_p_IFNB = NA_real_,
    eQTL_beta_TNF = NA_real_,  eQTL_p_TNF  = NA_real_,
    reQTL_dbeta_IFNG = NA_real_, reQTL

# ## ---------- 5) Parallel over pairs: benchmark timing with different n pairs ####
# TEST_N <- 1000
# npairs <- TEST_N
# message("Fitting ", npairs, " pairs in parallel (", future::nbrOfWorkers(), " workers) ...")
# t_start <- Sys.time()
# res_list <- future.apply::future_lapply(
#   seq_len(npairs),
#   function(i) fit_one(pairs$snp[i], pairs$gene[i]),
#   future.seed = TRUE
# )
# t_end <- Sys.time()
# message("Time for fitting ", npairs, " pairs: ", round(difftime(t_end, t_start, units="secs"),2), " seconds")


## ---------- 5) Parallel over pairs ---------- ####
message("Fitting ", npairs, " pairs in parallel (", future::nbrOfWorkers(), " workers) ...")
t_start <- Sys.time()
res_list <- future.apply::future_lapply(
  seq_len(npairs),
  function(i) fit_one(pairs$snp[i], pairs$gene[i]),
  future.seed = TRUE
)

res <- data.table::rbindlist(res_list, fill = TRUE)
data.table::fwrite(res, output_file, sep = "\t")
message("Wrote ", nrow(res), " rows to ", output_file)
t_end <- Sys.time()
message("====== QTL chunk done ======")
message("Time for fitting ", npairs, " pairs: ", round(difftime(t_end, t_start, units="mins"),2), "mins")