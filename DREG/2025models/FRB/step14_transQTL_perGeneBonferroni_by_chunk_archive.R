#!/usr/bin/env Rscript
# script for transQTL controlling within-gene FWER of eQTL or reQTL modeling stats (pval) using eigenMT

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
  Rscript step14_transQTL_perGeneBonferroni_by_chunk.R <ct> <condition> <QTLtype> <chunk_id>

Example:
  Rscript step14_transQTL_perGeneBonferroni_by_chunk.R KRT IFNG reQTL 076
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
condition  <- args[[2]]
QTLtype    <- args[[3]]
chunk_id       <- args[[4]] # 000 to 120
this_condition <- condition

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
pair_file <- paste0(dir,"/transQTL/eGene_QTL_pairs/",ct,"_",condition,"_",QTLtype,"_trans_pairs.txt")
modelstats_file <- paste0(dir,"/transQTL/results/",condition,"/",QTLtype,"/result_",chunk_id,".tsv")
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")
n_workers <- 2
base_seed <- 42

# # toy example for debugging
# ct         <- "FRB"
# condition  <- "IFNG"
# QTLtype    <- "eQTL"
# chunk_id <- "002"
# this_condition <- condition
# dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
# pair_file <- paste0(dir,"/transQTL/eGene_QTL_pairs/",ct,"_",condition,"_",QTLtype,"_trans_pairs.txt")
# modelstats_file <- paste0(dir,"/transQTL/results/",condition,"/",QTLtype,"/result_",chunk_id,".tsv")
# meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
# out_prefix <- paste(ct, condition, QTLtype, chunk_id, sep = "_")
# n_workers <- 2
# base_seed <- 42


message("=== within-gene bonferroni correction for ", ct, " ", condition, " ", QTLtype, ", chunk: ", chunk_id, " ===")

#### ----Load data sets -------- ####
#### --- Load chunk files --------- ####
# load modeling stats for the QTL of interest of chunk
modelstats_subset <- fread(modelstats_file) %>%
  dplyr::mutate(key = paste0(gene,"_",snp))
table(modelstats_subset$gene)

# load SNP:gene pairs
pairs <- fread(pair_file, header=T) 
#table(pairs$gene)

pairs <- pairs %>% 
  dplyr::mutate(key = paste0(gene,"_",snp)) %>%
  dplyr::filter(key %in% modelstats_subset$key)ls

# ----------- bonferroni function for 1 gene ------ ####
# compute bonferroni stats for one gene
run_one_gene_bonferroni <- function(this_gene,
                                 pairs, modelstats_subset,
                                 ct, condition, QTLtype) {
  
  # SNP list for this gene from pairs
  snps <- pairs %>% dplyr::filter(gene == this_gene) %>% dplyr::pull(snp) %>% unique()
  ntests <- length(snps)
  
  # Pull SNP-level p-values for this gene
  ms_g <- modelstats_subset %>%
    dplyr::filter(gene == this_gene, snp %in% snps) %>%
    dplyr::select(snp, gene, p, SNPtag,cis_trans_category) %>%
    dplyr::mutate(padj = pmin(p * ntests, 1))
  
  # lead SNP = SNP with smallest raw p
  lead_idx <- which.min(ms_g$p)

  tibble::tibble(
    celltype = ct,
    condition = condition,
    QTLtype = QTLtype,
    gene = this_gene,
    n_snps_total = ntests,
    lead_snp = ms_g$snp[lead_idx],
    pmin = ms_g$p[lead_idx],
    p_gene_bonferroni = ms_g$padj[lead_idx],
    lead_snp_tag = ms_g$SNPtag[lead_idx],
    lead_snp_cis_trans_category = ms_g$cis_trans_category[lead_idx]
  )
}

# ----------- MAIN ---------------------------- ####
# Main: Run eigenMT on all genes in chunk across 2 workers
# ------------------------------#
#### --- Process for every gene (chunk) --------- ####
start_time <- Sys.time()

# output folder and filename
out_dir <- paste0(dir, "/transQTL/bonferroni/", condition,"/", QTLtype)
out_file <- file.path(out_dir, paste0(out_prefix, ".bonferroni.tsv"))

future::plan(future::multisession, workers = n_workers)

gene_lst <- unique(pairs$gene)

results_list <- future.apply::future_lapply(
  gene_lst,
  function(g) {
    run_one_gene_bonferroni(
      this_gene = g,
      pairs = pairs,
      modelstats_subset = modelstats_subset,
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
