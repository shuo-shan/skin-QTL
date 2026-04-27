# QTL modeling QC filtering.
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
  stop("Usage: Rscript step3_QTL_result_QC.R <CELLTYPE> <condition> <QTLtype> <chunk_id>")
}

# take in arguments
ct <- args[1] # MEL, KRT, FRB
condition <- args[2] # PBS, IFNG, IFNB, TNF
QTLtype <- args[3] # eQTL, reQTL
chunk_id  <- args[4] # "000", "001", etc

# define variables
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
model_stats_file <- paste0(dir,"/results/result_",chunk_id,".tsv")

# # toy example for debugging
# ct <- "MEL"
# condition="IFNG"
# QTLtype="reQTL"
# chunk_id="006"
# dir=paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
# model_stats_file <- paste0(dir,"/results/result_",chunk_id,".tsv")

message("====== QTL statistics QC start ======")
message("Model Statistics Chunk: ", model_stats_file)

#### eQTL: load statistics table --> stats_relevant variable ####
if (QTLtype == "eQTL") {
  stats_all <- fread(model_stats_file, header = TRUE)
  # split qc string into 12 integer columns
  qc_mat <- tstrsplit(stats_all$qc_n_per_condition_by_genotype, "_", fixed = TRUE)
  qc_mat <- lapply(qc_mat, as.integer)
  # name them as n_g{0,1,2}_{COND}
  conds <- c("PBS","IFNG","IFNB","TNF")
  qc_names <- c(
    paste0("n_g0_", conds),
    paste0("n_g1_", conds),
    paste0("n_g2_", conds)
  )
  # add to table
  stats_all[, (qc_names) := qc_mat]
  
  # Select relevant columns from QTL statistics result table
  pattern <- paste0("^", QTLtype, "_.*_", condition, "$")
  qtl_cols <- grep(pattern, names(stats_all), value = TRUE)
  qtl_cols <- qtl_cols[!grepl("baseline", qtl_cols)]
  n_genotype_cols <- grep(paste0("^n_g.*",condition,"$"), names(stats_all), value = TRUE)
  cols_keep <- c("snp", "gene", "dist", qtl_cols, n_genotype_cols,"qc_genotypes")
  stats_relevant <- stats_all[, ..cols_keep]
  colnames(stats_relevant)[4:6] <- c("beta","p","se")
  colnames(stats_relevant)[7:9] <- c("n_geno0","n_geno1","n_geno2")
  
  CPM_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
  CPM <- data.table::fread(CPM_FILE) %>%
    dplyr::select(-name, -gene) %>%
    dplyr::filter(final_gene != 'ARMCX5-GPRASP2') %>%
    tibble::column_to_rownames("final_gene")
  
  pattern <- paste0(ct,"_", condition,"_")
  if (condition == "IFNG") {
    pattern <- paste0(ct, "_(", "IFN|IFNG", ")_")
  }
  selected_columns <- grep(pattern, names(CPM), value = TRUE)
  CPM <- CPM[ , selected_columns]
  
  cpm_thresh <- 1
  min_frac   <- 0.25
  n_samples <- ncol(CPM)
  min_n     <- ceiling(n_samples * min_frac)
  
  expressed_genes <- rownames(CPM)[rowSums(CPM >= cpm_thresh, na.rm = TRUE) >= min_n]
}

#### reQTL: select paired donors --> stats_relevant variable ####
if (QTLtype == "reQTL") {
  # for reQTLs: count number of donors that have both paired samples in PBS and stimulation
  chunk_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
  GENO_FILE <- paste0(dir,"/chunks/genotype_pairs_chunk_",chunk_id,".tsv") # columns: CHR POS ID REF ALT sample1 sample2 ...
  VST_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
  META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
  conds <- c("PBS","IFNG","IFNB","TNF")
  
  # Loading pairs from chunk
  pairs <- fread(chunk_file, header = TRUE)
  colnames(pairs)[4] <- "gene"
  colnames(pairs)[11] <- "snp"
  
  # Loading VST
  expr_wide <- fread(VST_FILE) %>%
    dplyr::filter(gene %in% pairs$gene) %>%
    dplyr::select(c("gene", contains(ct)))
  
  # Loading sample metadata
  selected_cond <- unique(c("PBS",condition))
  meta <- fread(META_FILE) %>% 
    dplyr::filter(celltype==ct) %>%
    dplyr::filter(condition %in% selected_cond)
  # Keep VST columns that have metadata (for this celltype)
  expr_samples <- intersect(colnames(expr_wide)[-1], meta$sample)
  expr_wide <- as.data.table(expr_wide)[, c("gene", expr_samples), with = FALSE]
  meta <- meta[sample %in% expr_samples]
  meta$condition <- factor(as.character(meta$condition), levels = conds)
  
  donor_cond_counts <- as.data.frame(table(meta$donor))
  paired_donors <- donor_cond_counts %>%
    dplyr::filter(Freq == 2) %>%
    pull(Var1) %>%
    as.character()
  
  # Melt VST to long (gene_id, sample, VST) and merge metadata
  expr_long <- melt(as.data.table(expr_wide), id.vars = "gene",
                    variable.name = "sample", value.name = "vst")
  expr_long <- merge(expr_long, meta, by = "sample", all.x = TRUE, all.y = FALSE)
  rm(expr_wide); invisible(gc())
  
  # Loading genotype table
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
  
  # Keep only genotype donors that appear in expression data
  geno_keep_samples <- intersect(geno_sample_cols, unique(expr_long$donor))
  geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples), with = FALSE]
  
  # for reQTL: Keep only genotype donors that have paired data
  geno_keep_samples_paired <- intersect(paired_donors, geno_keep_samples)
  geno_dt <- geno_dt[, c(fixed_cols, geno_keep_samples_paired), with = FALSE]
  
  # Read QTL statistics result table
  stats_all <- fread(model_stats_file, header = TRUE) %>%
    select(-qc_genotypes, -qc_n_per_condition_by_genotype)
  
  pattern <- paste0("^", QTLtype, "_.*_", condition, "$")
  qtl_cols <- grep(pattern, names(stats_all), value = TRUE)
  cols_keep <- c("snp", "gene", "dist", qtl_cols)
  
  stats_relevant <- stats_all[, ..cols_keep]
  colnames(stats_relevant)[4:6] <- c("beta","p","se")
  
  stats_relevant <- stats_relevant %>% left_join(geno_dt, by = c("snp" = "ID"))
  
  donor_cols <- intersect(paired_donors, colnames(stats_relevant))
  
  stats_relevant <- stats_relevant %>%
    mutate(
      n_geno0 = rowSums(across(all_of(donor_cols), ~ .x == 0), na.rm = TRUE),
      n_geno1 = rowSums(across(all_of(donor_cols), ~ .x == 1), na.rm = TRUE),
      n_geno2 = rowSums(across(all_of(donor_cols), ~ .x == 2), na.rm = TRUE),
      qc_genotypes = (n_geno0 > 0) + (n_geno1 > 0) + (n_geno2 > 0)
    ) %>%
    dplyr::select(c(snp,gene,dist,beta,p,se,n_geno0,n_geno1,n_geno2,qc_genotypes))
  
  # get expressed genes
  kept.samples <- meta %>% dplyr::filter(donor %in% donor_cols) %>% pull(sample)
  
  CPM_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
  CPM <- data.table::fread(CPM_FILE) %>%
    dplyr::select(-name, -gene) %>%
    dplyr::filter(final_gene != 'ARMCX5-GPRASP2') %>%
    tibble::column_to_rownames("final_gene") %>%
    dplyr::select(all_of(kept.samples))
  
  CPM.pbs <- CPM[ , grep(paste0(ct,"_PBS_"), names(CPM), value=TRUE)]
  
  pattern <- paste0(ct,"_", condition,"_")
  if (condition == "IFNG") {
    pattern <- paste0(ct, "_(", "IFN|IFNG", ")_")
  }
  selected_columns <- grep(pattern, names(CPM), value = TRUE)
  CPM.stim <- CPM[ , selected_columns]
  
  cpm_thresh <- 1
  min_frac   <- 0.25
  n_samples <- ncol(CPM.pbs)
  min_n     <- ceiling(n_samples * min_frac)
  
  expressed_genes.pbs <- rownames(CPM.pbs)[rowSums(CPM.pbs >= cpm_thresh, na.rm = TRUE) >= min_n]
  expressed_genes.stim <- rownames(CPM.stim)[rowSums(CPM.stim >= cpm_thresh, na.rm = TRUE) >= min_n]
  expressed_genes <- union(expressed_genes.pbs, expressed_genes.stim)
}

#### Apply QC filters on stats_relevant variable ####
# Rule A: Require all three genotypes globally #
stats_relevant.QCa  <- stats_relevant[qc_genotypes == 3]

# Rule B: Minimum donors per genotype is 3. if looking at reQTL, donor is required to be paired
# and Total samples per condition ≥ 15
stats_relevant.QCab <- stats_relevant.QCa[pmin(n_geno0, n_geno1, n_geno2, na.rm=TRUE) >= 3]
stats_relevant.QCab <- stats_relevant.QCab[(n_geno0 + n_geno1 + n_geno2) >= 15]

# Rule C: gene is expressed
stats_relevant.QCabc <- stats_relevant.QCab %>% dplyr::filter(gene %in% expressed_genes)

# # Rule D: z-score triage
# 2) What Z-score threshold should you use?
#   Since you’re using this as triage (not inference), pick something that:
#   shrinks compute meaningfull
# doesn’t throw away moderate-but-real effects
# I’d recommend starting with:
#   |Z| ≥ 2.5 (moderate screen; keeps more)
# To reduce computational burden, we applied a screening step based on the maximum absolute Z-score across cis-SNPs for each gene. 
# Genes with max |Z| ≥ 2.5 were retained for permutation testing. 
# This screening step was used solely to allocate permutation resources and does not constitute inference. 
# Empirical FDR was therefore estimated conditional on the screened gene set.

# Also: compute Z for the relevant statistic:
#   eQTL screen: z = beta / se per condition
# reQTL screen: z = dbeta / se per condition (interaction term)
# Then per gene:
#   max_abs_z = max(|z|) across SNPs in the cis window.
selected_genes <- stats_relevant.QCabc %>%
  dplyr::mutate(z = beta/se) %>%
  dplyr::filter( abs(z) >= 2.5) %>%
  pull(gene) %>%
  unique()

stats_relevant.QCabcd <- stats_relevant.QCabc %>% 
  dplyr::filter(gene %in% selected_genes) %>%
  dplyr::mutate(z = beta/se)

#### write output to file ####
outfile <- paste0(dir,"/results_QC/modeling_stats_postQC_",ct,"_",condition,"_",QTLtype,"_",chunk_id,".txt")
data.table::fwrite(stats_relevant.QCabcd, 
                   file = outfile,
                   quote = F, sep = "\t", row.names = F, col.names = T)

message(paste0("Done QC filtering for: ",outfile))