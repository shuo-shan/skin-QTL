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
# true negative but identified as false positive in TNF reQTL (p=0, beta=-235)
#gene="CALD1"; snp="rs759854"
this_gene="LIPG"; this_snp="rs8096411"
# true positive and also detected positive in IFNG reQTL 
#gene="GBP3"; snp="rs72726528"

dir="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL"
setwd(dir)
chunk_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
#GENO_FILE <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv") # columns: CHR POS ID REF ALT sample1 sample2 ...
GENO_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots/genotype_debug.tsv"
output_file <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/debug/debug_",this_gene,"_",this_snp,".tsv")

## ------------------ Paths ------------------
# Expression + metadata
CPM_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
VST_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
PEER_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_MEL_PBS-IFNG-IFNB-TNF.tsv"             # columns: sample, PEER1, PEER2, ...
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"      # columns: donor, PC1, PC2, ...

message("====== QTL fitting start ======")
message("Gene: ",this_gene,"\n", "SNP: ",this_snp,"\n","Output: ", output_file)

## ---------- 0) Settings ----------
# Parallel workers (10 workers is pre-determined through benchmarking)
plan(multicore, workers = 10)
set.seed(1)
conds_all <- c("PBS","IFNG","IFNB","TNF")  # allowed condition labels

## ---------- 1) Load once ---------- ####
message("[2/7] Loading expression table ...")
expr_wide <- fread(VST_FILE) %>%
  dplyr::filter(gene==this_gene) %>%
  dplyr::select(c("gene", contains(ct)))

expr_wide <- fread(CPM_FILE) %>%
  dplyr::filter(gene==this_gene) %>%
  dplyr::select(c("gene", contains(ct)))


message("[3/7] Loading sample metadata ...")
meta <- fread(META_FILE) %>% dplyr::filter(celltype==ct)
# Keep RNA expression columns that have metadata (for this celltype)
expr_samples <- intersect(colnames(expr_wide)[-1], meta$sample)
expr_wide <- as.data.table(expr_wide)[, c("gene", expr_samples), with = FALSE]
meta <- meta[sample %in% expr_samples]
meta$condition <- factor(as.character(meta$condition), levels = conds_all)


# Melt expression table to long (gene_id, sample, expression) and merge metadata
message("[4/7] Reshaping expression table and merging metadata ...")
expr_long <- melt(as.data.table(expr_wide), id.vars = "gene",
                  variable.name = "sample", value.name = "expression")
expr_long <- merge(expr_long, meta, by = "sample", all.x = TRUE, all.y = FALSE)
rm(expr_wide); invisible(gc())


# PEER factors at sample level
message("[5/7] Loading PEER (sample-level) ...")
peers <- fread(PEER_FILE)
colnames(peers)[1] <- "sample"
peer_cols <- setdiff(colnames(peers), "sample")
# Merge PEERs onto sample rows directly (sample-level PEER merge)
expr_long <- merge(expr_long, peers, by = "sample", all.x = TRUE)
# Drop rows without PEERs if any remain
expr_long <- expr_long[complete.cases(expr_long[, ..peer_cols]), ]


# Genotype PCs (per donor)
message("[6/7] Loading genotype PCs ...")
gpcs <- fread(GENO_PCS_FILE)
pc_cols <- setdiff(colnames(gpcs), "donor")
expr_long <- merge(expr_long, gpcs, by = "donor", all.x = TRUE)

# Genotype table
message("[7/7] Loading genotype table ...")
# before this step, in cluster, folder/debug/, run get_genotype.sh
geno_dt <- fread(paste0(dir,"/debug/filtered_",this_snp,".txt")) %>%
  dplyr::select(-c(QUAL,FILTER,INFO,FORMAT))
colnames(geno_dt) <- sapply(colnames(geno_dt), function(x) strsplit(x, "_")[[1]][1])
colnames(geno_dt)[1] <- "CHROM"
geno_cols <- colnames(geno_dt)[6:ncol(geno_dt)]
geno_dt[, (geno_cols) := lapply(.SD, function(x) {
  # map GT strings to 0/1/2
  x <- gsub("0/0", "0", x)
  x <- gsub("0/1", "1", x)
  x <- gsub("1/0", "1", x)
  x <- gsub("1/1", "2", x)
  as.integer(x)
}), .SDcols = geno_cols]

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

# Keep only genotype donors that appear in expression data
geno_keep_samples <- intersect(geno_sample_cols, unique(expr_long$donor))
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
fml_str <- paste0("expression ~ geno * condition ", cov_part, " + (1 | donor)")
message("Model formula: ", fml_str)

## ---------- 4) Per-pair fit function ---------- ####
conds_all <- c("PBS","IFNG","IFNB","TNF")

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
  snp = this_snp
  gene = this_gene
  cols_needed <- c("sample","donor","condition","expression", peer_cols, pc_cols)
  df <- tryCatch(expr_long[gene == this_gene, ..cols_needed], error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  
  df$condition <- factor(as.character(df$condition))
  present_conds <- intersect(conds_all, unique(df$condition))
  df <- df[condition %in% present_conds]
  df$condition <- stats::relevel(df$condition, ref="PBS")
  df$geno <- fetch_geno_for_snp(snp, df$donor)
  
  ug <- na.omit(unique(df$geno))
  if (length(ug) < 2) return(NULL)  # need ≥2 genotype groups
  
  # require at least 5 donors in each genotype group in each condition
  tab_gc <- table(df$geno, df$condition)
  min_per_cond <- apply(tab_gc, 2, min)
  cutoff <- 5
  keep_conds <- names(min_per_cond)[min_per_cond >= cutoff]
  if (!"PBS" %in% keep_conds) {
    return(NULL)  # baseline too thin → skip this SNP–gene pair
  }
  if (length(keep_conds) < 2) {
    return(NULL)
  }
  df_sub <- df[df$condition %in% keep_conds, ]
  
  
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
  tr_summ <- as.data.frame(summary(tr, infer = c(TRUE, TRUE)))
  
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
    snp = this_snp, gene = this_gene,
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
)


#### Leave One Out test of Donor Stability ####
library(lme4)
library(emmeans)

lodo_reQTL <- function(df_sub, fml_str, test_cond = "IFNG") {
  
  # ---- helper: fit model + extract beta & p for PBS vs test_cond ----
  get_beta_p <- function(dat) {
    # need PBS & test_cond present
    if (!all(c("PBS", test_cond) %in% unique(dat$condition))) {
      return(c(beta = NA_real_, p = NA_real_))
    }
    
    fml <- as.formula(fml_str)
    fit <- tryCatch(
      lmer(fml, data = dat, REML = FALSE,
           control = lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(c(beta = NA_real_, p = NA_real_))
    
    tr <- tryCatch(
      emtrends(fit, ~ condition, var = "geno"),
      error = function(e) NULL
    )
    if (is.null(tr)) return(c(beta = NA_real_, p = NA_real_))
    
    re_con <- tryCatch(
      contrast(tr, "trt.vs.ctrl", ref = "PBS"),
      error = function(e) NULL
    )
    if (is.null(re_con)) return(c(beta = NA_real_, p = NA_real_))
    
    re_tbl <- as.data.frame(summary(re_con))
    
    # row corresponding to PBS vs test_cond
    row <- re_tbl[grepl(test_cond, re_tbl$contrast), , drop = FALSE]
    if (!nrow(row)) return(c(beta = NA_real_, p = NA_real_))
    
    c(beta = row$estimate, p = row$p.value)
  }
  
  donors <- sort(unique(df_sub$donor))
  
  # results table
  res <- data.frame(
    donor_left_out = c("NONE", donors),
    beta_IFNG = NA_real_,
    p_IFNG    = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # full data
  res[1, c("beta_IFNG", "p_IFNG")] <- get_beta_p(df_sub)
  
  # leave-one-donor-out
  for (i in seq_along(donors)) {
    d <- donors[i]
    dat_i <- df_sub[df_sub$donor != d, ]
    res[i + 1, c("beta_IFNG", "p_IFNG")] <- get_beta_p(dat_i)
  }
  
  res
}
lodo_res <- lodo_reQTL(df, fml_str, "IFNG")

#### plots ####
# ---------------------- Plot Function ----------------------
# ---------------------- Plot Helper (PBS vs condition panel) ----------------------
make_subpanel <- function(df, ct, stim, snp, gene){
  
  sub <- df %>% filter(celltype == ct, condition %in% c("PBS", stim))
  if(nrow(sub) == 0) return(NULL)
  
  # Count donors per genotype + condition
  counts <- sub %>%
    group_by(genotype, dosage, condition) %>%
    summarize(n = n_distinct(donor), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
    mutate(strip_label = paste0(genotype, "\nPBS=", PBS, " ", stim, "=", !!sym(stim)))
  
  # join counts into sub
  sub <- sub %>%
    left_join(counts, by = c("genotype", "dosage"))
  
  # ---- enforce facet order by numeric dosage ----
  sub <- sub %>%
    arrange(dosage) %>%   # ensures 0 → 1 → 2 before converting to factor
    mutate(
      genotype_label = factor(strip_label, levels = unique(strip_label))
    )
  
  # ensure PBS appears left on x-axis
  sub$condition <- factor(sub$condition, levels = c("PBS", stim))
  
  
  # ---- plotting ----
  st <- get_stats(snp, gene, ct, stim)
  subtitle <- paste0(ct, " — PBS vs ", stim, "  (", fmt_stat(st$beta, st$p), ")")
  
  ggplot(sub, aes(x = condition, y = CPM, group = donor)) +
    geom_line(alpha = 0.25, color = "grey70") +
    geom_point(aes(color = condition), size = 2.2, alpha = 0.9) +
    scale_color_manual(values = c(
      "PBS" = "gray35",
      "IFNG" = "#0072B2",
      "IFNB" = "#E69F00",
      "TNF" = "#D55E00"
    )) +
    facet_wrap(~ genotype_label, scales = "fixed", nrow = 1) +
    labs(title = subtitle, y = "CPM", x = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 3),
      strip.text = element_text(size = 8, lineheight = 0.8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.major.x = element_blank(),   # already removing vertical major grid
      panel.grid.major.y = element_blank(),   # ← removes horizontal major grid
      panel.grid.minor.y = element_blank()    # ← removes horizontal minor grid
    )
  
}


# ---------------------- Main Plot for a SNP : gene ----------------------
make_plot <- function(snp, gene){
  
  df <- expr_long %>% filter(gene_id == gene)
  geno <- get_genotype_df(snp)
  df <- df %>% inner_join(geno, by="donor")
  
  ct_order <- c("FRB","KRT","MEL")
  stim_order <- c("IFNG","IFNB","TNF")
  
  panels <- list()
  for(ct in ct_order){
    for(stim in stim_order){
      p <- make_subpanel(df, ct, stim, snp, gene)
      if(!is.null(p)) panels[[paste(ct, stim)]] <- p
    }
  }
  
  combined <- wrap_plots(panels, ncol=3)
  
  # Apply global title WITHOUT overriding sub-panel title sizes
  combined <- combined + plot_annotation(
    title = paste0(snp, " : ", gene)
  )
  
  combined & theme(
    plot.title = element_text(size=14, face="bold", hjust=0.5)  # affects ONLY main title now
  )
}

# ---------------------- Output ----------------------
# message("[6] Writing PDF: ", OUT_PDF)
# pdf(OUT_PDF, width=16, height=8)
# #for(i in seq_len(nrow(pairs))){
# for(i in seq_len(10)){
#   print(make_plot(pairs$snp[i], pairs$gene[i]))
# }
# dev.off()
# message("Done")