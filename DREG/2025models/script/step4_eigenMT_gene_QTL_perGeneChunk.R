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
  Rscript step4_eigenMT_gene_QTL_perGeneChunk.R <ct> <condition> <QTLtype> <chunk_id>

Example:
  Rscript step4_eigenMT_gene_QTL_perGeneChunk.R MEL IFNG reQTL 076
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
condition  <- args[[2]]
QTLtype    <- args[[3]]
chunk_id       <- args[[4]] # 000 to 120
this_condition <- condition

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
pair_file <- paste0(dir,"/data/chunk/pair_chunk_",chunk_id,".txt")
geno_file  <- paste0(dir,"/data/chunk/genotype_pair_chunk_",chunk_id,".txt")
modelstats_file <- paste0(dir,"/permutation/model_stats/model_stats_",chunk_id,".txt")
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")
n_workers <- 2
base_seed <- 42

# # toy example for debugging
# ct         <- "MEL"
# condition  <- "IFNG"
# QTLtype    <- "reQTL"
# chunk_id <- "034"
# dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
# pair_file <- paste0(dir,"/data/chunk/pair_chunk_",chunk_id,".txt")
# geno_file  <- paste0(dir,"/data/chunk/genotype_pair_chunk_",chunk_id,".txt")
# modelstats_file <- paste0(dir,"/permutation/model_stats/model_stats_",chunk_id,".txt")
# meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
# out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")
# n_workers <- 2
# base_seed <- 42
# this_condition <- condition

message("=== eigenMT for ", ct, " ", condition, " ", QTLtype, ", chunk: ", chunk_id, " ===")

#### ----Load data sets -------- ####
#### --- Load chunk files --------- ####
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

#### --- Load metadata to define LD donor set --- ####
meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  dplyr::filter(celltype == ct, condition %in% unique(c("PBS", this_condition))) %>%
  dplyr::select(sample, donor, condition)

# Define donors for eQTL and paired donors for reQTL

if (QTLtype=="reQTL") {
  donor_counts <- meta %>%
    dplyr::count(donor, condition, name = "n") %>%
    dplyr::filter(n >= 1)
  
  paired_donors <- donor_counts %>%
    dplyr::filter(condition %in% c("PBS", this_condition)) %>%
    dplyr::group_by(donor) %>%
    dplyr::summarise(n_cond = dplyr::n_distinct(condition), .groups = "drop") %>%
    dplyr::filter(n_cond == 2) %>%
    dplyr::pull(donor) %>%
    as.character()
  
  if (length(paired_donors) < 5) {
    stop("Too few paired donors for LD: ", length(paired_donors))
  }
  donors_ld = paired_donors
} else {
  donors_ld = unique(meta$donor)
}

# ----------- eigenMT function for 1 gene ------ ####
# compute eigenMT stats for one gene
run_one_gene_eigenMT <- function(this_gene,
                                 pairs, genotype_all, modelstats_subset, p_col,
                                 donors_ld,
                                 ct, condition, QTLtype) {
  
  # SNP list for this gene from pairs
  snps <- pairs %>% dplyr::filter(gene == this_gene) %>% dplyr::pull(SNP) %>% unique()
  
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
  
  # keep donor order stable (not crucial for cor, but good hygiene)
  donors_use <- intersect(donors_ld, geno_sample_cols)
  donors_use <- donors_ld[donors_ld %in% donors_use]
  
  if (length(donors_use) < 5) {
    return(tibble(
      celltype = ct, condition = condition, QTLtype = QTLtype, gene = this_gene,
      n_snps_total = length(snps), n_snps_used = length(pvals), n_donors_ld = length(donors_use),
      pmin = pmin, lead_snp = lead, Meff = NA_real_, p_gene_eigenMT = NA_real_,
      stage = "too_few_ld_donors"
    ))
  }
  
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
out_dir <- paste0(dir, "/eigenMT/", condition,"/", QTLtype)
out_file <- file.path(out_dir, paste0(out_prefix, ".eigenMT.tsv"))

future::plan(future::multicore, workers = n_workers)

gene_lst <- unique(pairs$gene)

results_list <- future.apply::future_lapply(
  gene_lst,
  function(g) {
    run_one_gene_eigenMT(
      this_gene = g,
      pairs = pairs,
      genotype_all = genotype_all,
      modelstats_subset = modelstats_subset,
      p_col = paste0(QTLtype,"_p_",condition),
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