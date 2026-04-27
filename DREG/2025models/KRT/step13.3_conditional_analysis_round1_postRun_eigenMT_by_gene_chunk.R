#!/usr/bin/env Rscript
# script for controlling within-gene FWER of eQTL or reQTL modeling stats (pval) using eigenMT after conditioning on lead SNP

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
  library(corpcor)
})

#### Input ####
# ------------------------------#
# Main: parse CLI args ####
# ------------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("
Usage:
  Rscript step13.3_conditional_analysis_round1_postRun_eigenMT_by_gene_chunk.R <ct> <condition> <QTLtype> <chunk_id>

Example:
  Rscript step13.3_conditional_analysis_round1_postRun_eigenMT_by_gene_chunk.R KRT IFNG reQTL 076
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
condition  <- args[[2]]
QTLtype    <- args[[3]]
chunk_id       <- args[[4]] # 000 to 120
this_condition <- condition

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
pair_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
geno_file  <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv")
modelstats_file <- paste0(dir,"/conditional_analysis_round1/results_QC/",this_condition,"/",QTLtype,"/modeling_stats_postQC_",ct,"_",condition,"_",QTLtype,"_",chunk_id,".txt")
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")
n_workers <- 2
base_seed <- 42
# 
# # toy example for debugging
# ct         <- "FRB"
# condition  <- "IFNB"
# QTLtype    <- "eQTL"
# chunk_id <- "034"
# this_condition <- condition
# dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
# pair_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
# geno_file  <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv")
# modelstats_file <- paste0(dir,"/conditional_analysis_round1/results_QC/",this_condition,"/",QTLtype,"/modeling_stats_postQC_",ct,"_",condition,"_",QTLtype,"_",chunk_id,".txt")
# meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
# out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")
# n_workers <- 2
# base_seed <- 42

message("=== eigenMT for ", ct, " ", condition, " ", QTLtype, ", chunk: ", chunk_id, " ===")

#### ----Load data sets -------- ####
#### --- Load chunk files --------- ####
# load modeling stats for the QTL of interest of chunk
if (!file.exists(modelstats_file)) {
  warning("Input file does not exist: ", modelstats_file)
  modelstats_subset <- data.table()
} else {
  modelstats_subset <- tryCatch(
    fread(modelstats_file, header = TRUE),
    error = function(e) {
      warning("Failed to read file: ", modelstats_file)
      data.table()
    }
  )
}

# handle empty input
if (nrow(modelstats_subset) == 0) {
  message("Input file is missing/empty. Writing empty eigenMT output.")
  
  out_dir <- paste0(dir, "/conditional_analysis_round1/eigenMT/", condition, "/", QTLtype)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(out_prefix, ".eigenMT.tsv"))
  
  empty_dt <- data.table(
    celltype = character(),
    condition = character(),
    QTLtype = character(),
    gene = character(),
    n_snps_total = integer(),
    n_snps_used = integer(),
    n_donors_ld = integer(),
    pmin = numeric(),
    lead_snp = character(),
    Meff = numeric(),
    p_gene_eigenMT = numeric(),
    stage = character()
  )
  
  fwrite(
    empty_dt,
    file = out_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
  message("Wrote empty output: ", out_file)
  quit(save = "no")
}

modelstats_subset <- modelstats_subset %>%
  dplyr::mutate(key = paste0(gene, "_", snp))

table(modelstats_subset$gene)

# load SNP:gene pairs
pairs <- fread(pair_file, header=T) 
table(pairs$gene_name)

pairs <- pairs %>% 
  dplyr::mutate(key = paste0(gene_name,"_",SNP_ID)) %>%
  dplyr::filter(key %in% modelstats_subset$key)

# load genotype of all SNPs that passed QC in chunk
genotype_all <- fread(geno_file) %>%
  dplyr::filter(ID %in% unique(modelstats_subset$snp))

donors_with_genotype <- fread(geno_file) %>% dplyr::select(-c(CHROM, POS, ID, REF, ALT)) %>% colnames() %>% unique()

#### --- Load metadata to define LD donor set --- ####
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  dplyr::filter(celltype==ct) %>%
  dplyr::select(sample, donor, condition)
donors_ld = unique(meta$donor)

# ----------- eigenMT function for 1 gene ------ ####
# compute eigenMT stats for one gene
run_one_gene_eigenMT <- function(this_gene,
                                 pairs, genotype_all, modelstats_subset, p_col="p",
                                 donors_ld,
                                 ct, condition, QTLtype) {
  
  # SNP list for this gene from pairs
  snps <- pairs %>% dplyr::filter(gene_name == this_gene) %>% dplyr::pull(SNP_ID) %>% unique()
  
  # Pull SNP-level p-values for this gene
  ms_g <- modelstats_subset %>%
    dplyr::filter(gene == this_gene, snp %in% snps) %>%
    dplyr::select(snp, gene, dplyr::all_of(p_col))
  
  pvals <- suppressWarnings(as.numeric(ms_g[[p_col]]))
  names(pvals) <- ms_g$snp
  pvals <- pvals[is.finite(pvals) & !is.na(pvals)]
  
  # Build LD matrix from genotype_all for SNPs with p-values
  fixed_cols <- c("CHROM","POS","ID","REF","ALT")
  geno_sample_cols <- setdiff(colnames(genotype_all), fixed_cols)
  donors_use <- intersect(donors_ld, geno_sample_cols)
  
  snps_use <- intersect(names(pvals), genotype_all$ID)

  geno_sub <- genotype_all %>%
    dplyr::filter(ID %in% snps_use) %>%
    dplyr::select(ID, dplyr::all_of(donors_use))
  
  # enforce SNP order (match snps_use)
  geno_sub <- geno_sub[match(snps_use, geno_sub$ID), , drop = FALSE]
  
  G_snp_by_donor <- as.matrix(dplyr::select(geno_sub, -ID))
  storage.mode(G_snp_by_donor) <- "numeric"
  
  # drop SNPs with zero variance across donors (causes NA correlations)
  sds <- apply(G_snp_by_donor, 1, stats::sd, na.rm = TRUE)
  keep <- which(is.finite(sds) & sds > 0)
  
  if (length(keep) == 0) {
    return(tibble(
      celltype = ct, condition = condition, QTLtype = QTLtype, gene = this_gene,
      n_snps_total = length(snps), n_snps_used = 0L, n_donors_ld = length(donors_use),
      pmin = NA_real_, lead_snp = NA_character_, Meff = NA_real_, p_gene_eigenMT = NA_real_,
      stage = "no_var_snps"
    ))
  }
  
  # Final SNP list after genotype+variance filtering
  snps_final <- geno_sub$ID[keep]
  
  # Recompute pmin and lead_snp on FINAL SNP set
  pvals_final <- pvals[snps_final]
  pmin <- min(pvals_final, na.rm = TRUE)
  lead <- names(which.min(pvals_final))
  
  # Restrict genotype matrix to final SNP set
  G_snp_by_donor <- G_snp_by_donor[keep, , drop = FALSE]
  
  # If <2 SNPs remain, Meff ~ 1 and Sidak reduces to pmin
  if (nrow(G_snp_by_donor) < 2) {
    Meff <- 1
    p_gene <- 1 - (1 - pmin)^Meff
    return(tibble(
      celltype = ct, condition = condition, QTLtype = QTLtype, gene = this_gene,
      n_snps_total = length(snps), n_snps_used = nrow(G_snp_by_donor), n_donors_ld = length(donors_use),
      pmin = pmin, lead_snp = lead, Meff = Meff, p_gene_eigenMT = p_gene,
      stage = "lt2_var_snps"
    ))
  }
  
  # correlation across SNPs using donors as observations
  # # Option 1: Apply Ledoit-Wolf Shrinkage for regularization
  # this is horrible because it causes more problem than solving the problem
  # G_transposed <- t(G_snp_by_donor)
  # R_shrunk <- corpcor::cor.shrink(G_transposed, verbose = FALSE)
  # lambda_used <- attr(R_shrunk, "lambda")
  # ev <- eigen(R_shrunk, symmetric = TRUE, only.values = TRUE)$values
  # Meff <- sum(pmin(ev, 1))
  # p_gene <- 1 - (1 - pmin)^Meff
  # 
  # # Option 2: Apply tiny jitter to identity matrix for shrinkage to ensure matrix is definite
  # # this artificially inflates Meff, not good.
  # R <- stats::cor(t(G_snp_by_donor), use = "pairwise.complete.obs")
  # diag(R) <- diag(R) + 1e-6
  # ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  # ev <- pmax(ev, 0) 
  # Meff <- sum(pmin(ev, 1))
  # p_gene <- 1 - (1 - pmin)^Meff
  
  # Option 3: no regularization step
  # SNP x donor -> cor across SNP => cor(t(G))
  R <- suppressWarnings(stats::cor(t(G_snp_by_donor), use = "pairwise.complete.obs"))
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values # without regularization, this step results in negative eigen values.
  ev <- pmax(ev, 0) # turn negative eigenvalues to 0
  Meff <- sum(pmin(ev, 1))
  p_gene <- 1 - (1 - pmin)^Meff

  tibble(
    celltype = ct,
    condition = condition,
    QTLtype = QTLtype,
    gene = this_gene,
    n_snps_total = length(snps),
    n_snps_used = nrow(G_snp_by_donor),
    n_donors_ld = length(donors_use),
    pmin = pmin,
    lead_snp = lead,
    Meff = Meff,
    p_gene_eigenMT = p_gene,
    stage = "ok"
  )
}

# ----------- MAIN ---------------------------- ####
# Main: Run eigenMT on all genes in chunk across 2 workers
# ------------------------------#
#### --- Process for every gene (chunk) --------- ####
start_time <- Sys.time()

# output folder and filename
out_dir <- paste0(dir, "/conditional_analysis_round1/eigenMT/", condition,"/", QTLtype)
out_file <- file.path(out_dir, paste0(out_prefix, ".eigenMT.tsv"))

future::plan(future::multisession, workers = n_workers)

gene_lst <- unique(pairs$gene_name)

results_list <- future.apply::future_lapply(
  gene_lst,
  function(g) {
    run_one_gene_eigenMT(
      this_gene = g,
      pairs = pairs,
      genotype_all = genotype_all,
      modelstats_subset = modelstats_subset,
      p_col = "p",
      donors_ld = donors_ld,
      ct = ct,
      condition = condition,
      QTLtype = QTLtype
    )
  },
  future.seed = TRUE
)

summary_out <- dplyr::bind_rows(results_list)

end_time <- Sys.time()

cat("Start:", format(start_time), "\n")
cat("End:  ", format(end_time), "\n")
cat("Elapsed:", as.numeric(difftime(end_time, start_time, units = "secs")), "seconds\n")

readr::write_tsv(summary_out, out_file)
message("Wrote: ", out_file)
