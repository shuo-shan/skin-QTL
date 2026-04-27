#!/usr/bin/env Rscript
# script for reclassifying QTL genes based on QTL dynamics across conditions.

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
ct         <- args[[1]]
this_condition  <- args[[2]]
this_gene <- args[[3]]

# # toy example
# ct <- "FRB"
# this_condition <- "IFNG"
# this_gene <- "ERAP2"

message(sprintf(
  "ct=%s gene=%s cytokine=%s",
  ct, this_gene, this_condition
))

# ----- Global paths ------
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)

chunk_id_lookup <- data.table::fread(paste0(dir, "/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[chunk_id_lookup$gene == this_gene, ]$chunk)
chunk_id <- sprintf("%03d", chunk_id)

pair_file <- paste0(dir, "/chunks/pairs_chunk_", chunk_id, ".tsv")
geno_file <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
vst_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
modelstats_file <- paste0(dir, "/results/result_", chunk_id, ".tsv")
meta_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"


qval_PBS_eQTL_file <- paste0(dir,"/eigenMT/results/",ct,"_","PBS","_","eQTL",".eigenMT.txt")
qval_stim_eQTL_file <- paste0(dir,"/eigenMT/results/",ct,"_",this_condition,"_","eQTL",".eigenMT.txt")
qval_stim_reQTL_file <- paste0(dir,"/eigenMT/results/",ct,"_",this_condition,"_","reQTL",".eigenMT.txt")
rm(chunk_id_lookup)

# ----- 1. Load pairs / genotype / meta / modelstats -------
pairs <- fread(pair_file, header = TRUE) %>%
  dplyr::filter(gene_name == this_gene) %>%
  dplyr::mutate(
    SNP = stringr::str_to_lower(SNP_ID),
    key = paste0(gene_name, "_", SNP_ID)
  )

genotype_all <- fread(geno_file) %>%
  dplyr::mutate(ID = stringr::str_to_lower(ID)) %>%
  dplyr::filter(ID %in% pairs$SNP_ID)

meta_all <- readr::read_tsv(meta_file, show_col_types = FALSE)

modelstats <- fread(modelstats_file) %>%
  dplyr::filter(gene == this_gene) %>%
  dplyr::filter(snp %in% unique(pairs$SNP_ID))

modelstats.PBS.eQTL <- modelstats %>%
  dplyr::select(snp, gene, contains("PBS")) %>%
  set_colnames(c("snp","gene","baseline","beta","p","se"))

modelstats.stim.eQTL <- modelstats %>%
  dplyr::select(snp, gene, contains(this_condition) & starts_with("eQTL")) %>%
  set_colnames(c("snp","gene","baseline","beta","p","se"))

modelstats.stim.reQTL <- modelstats %>%
  dplyr::select(snp, gene, contains(this_condition) & starts_with("reQTL")) %>%
  set_colnames(c("snp","gene","beta","p","se"))

# ----- 2. pick lead SNPs -----
# PBS eQTL
qval_table <- fread(qval_PBS_eQTL_file) %>%
  dplyr::filter(gene == this_gene)

if (nrow(qval_table) > 0) {
  # get eigenMT gene values
  res <- data.frame(celltype = qval_table$celltype[1],
              condition = qval_table$condition[1],
              QTLtype = qval_table$QTLtype[1],
              sig      = qval_table$q_gene[1] < 0.05,
              lead_snp = qval_table$lead_snp[1],
              lead_snp_p     = qval_table$pmin[1],
              lead_snp_q        = qval_table$q_gene[1],
              lead_snp_source = "PBS_eQTL"
              )
  
} else {
  # check modelstats SNP gene pair file
  res <- data.frame(celltype = ct,
              condition = "PBS",
              QTLtype = "eQTL",
              sig      = FALSE,
              lead_snp = modelstats.PBS.eQTL$snp[1],
              lead_snp_p     = modelstats.PBS.eQTL$p[1],
              lead_snp_q        = NA,
              lead_snp_source = "PBS_eQTL")
}
data_PBS_eQTL <- res

# stim eQTL 
qval_table <- fread(qval_stim_eQTL_file) %>%
  dplyr::filter(gene == this_gene)

if (nrow(qval_table) > 0) {
  # get eigenMT gene values
  res <- data.frame(celltype = qval_table$celltype[1],
              condition = qval_table$condition[1],
              QTLtype = qval_table$QTLtype[1],
              sig      = qval_table$q_gene[1] < 0.05,
              lead_snp = qval_table$lead_snp[1],
              lead_snp_p     = qval_table$pmin[1],
              lead_snp_q        = qval_table$q_gene[1],
              lead_snp_source = paste0(this_condition,"_eQTL"))
  
} else {
  # check modelstats SNP gene pair file
  res <- data.frame(celltype = ct,
              condition = this_condition,
              QTLtype = "eQTL",
              sig      = FALSE,
              lead_snp = modelstats.stim.eQTL$snp[1],
              lead_snp_p     = modelstats.stim.eQTL$p[1],
              lead_snp_q        = NA,
              lead_snp_source = paste0(this_condition,"_eQTL"))
}
data_stim_eQTL <- res

# stim reQTL 
qval_table <- fread(qval_stim_reQTL_file) %>%
  dplyr::filter(gene == this_gene)

if (nrow(qval_table) > 0) {
  # get eigenMT gene values
  res <- data.frame(celltype = qval_table$celltype[1],
                    condition = qval_table$condition[1],
                    QTLtype = "reQTL",
                    sig      = qval_table$q_gene[1] < 0.05,
                    lead_snp = qval_table$lead_snp[1],
                    lead_snp_p     = qval_table$pmin[1],
                    lead_snp_q        = qval_table$q_gene[1],
                    lead_snp_source = paste0(this_condition,"_reQTL"))
  
} else {
  # check modelstats SNP gene pair file
  res <- data.frame(celltype = ct,
                    condition = this_condition,
                    QTLtype = "reQTL",
                    sig      = FALSE,
                    lead_snp = modelstats.stim.reQTL$snp[1],
                    lead_snp_p     = modelstats.stim.reQTL$p[1],
                    lead_snp_q        = NA,
                    lead_snp_source = paste0(this_condition,"_reQTL"))
}
data_stim_reQTL <- res
rm(qval_table, res)

# summarize lead SNPs 
data_lead_snp <- rbind(data_PBS_eQTL,rbind(data_stim_eQTL, data_stim_reQTL))
data_all <- left_join(data_lead_snp, modelstats, by=c("lead_snp"="snp"))

# ----- 3. collapse lead SNPs if high LD -----
# Compute LD r^2 between two SNPs from a genotype table
ld_r2 <- function(genotype_all, snp1, snp2) {
  
  fixed_cols <- intersect(c("CHROM", "POS", "ID", "REF", "ALT"), colnames(genotype_all))
  donor_cols <- setdiff(colnames(genotype_all), fixed_cols)
  genotype_donors <- genotype_all %>% 
    column_to_rownames("ID") %>%
    dplyr::select(matches(donor_cols))
  
  # Extract donor genotypes (as numeric vectors)
  g1 <- as.numeric(genotype_donors[snp1, ])
  g2 <- as.numeric(genotype_donors[snp2, ])
  
  # Pairwise complete observations
  ok <- is.finite(g1) & is.finite(g2)
  g1 <- g1[ok]
  g2 <- g2[ok]
  n_used <- length(g1)
  
  if (stats::sd(g1) == 0 || stats::sd(g2) == 0) {
    return(list(r2 = NA_real_, r = NA_real_, n_used = n_used))
  }
  
  r <- stats::cor(g1, g2)
  list(r2 = r^2, n_used = n_used, snp1 = snp1, snp2 = snp2)
}

# Pick one lead SNP to represent the PBS or stimulated condition
# if they are all in the same LD with high r2 (r2>=0.6), only pick one lead snp and annotate the anchor sources and collapse sources with | delimiter
# for any pair of lead SNPs with low LD, then treat them separately and annotate the anchor sources respectively
# this lead snp only matters when the gene is sig in both PBSeQTL and stimeQTL (constitutive), and
# the lead snp's modelstats could tell us the dynamics of change, stable/amplified/attenuated/switched
# in the rare case where top PBS eQTL and top stim eQTL are in different LD, then the dynamics of change can also be interesting to know.
lead_snps <- list(
  PBS_eQTL = data_all[which(data_all$lead_snp_source=="PBS_eQTL"),]$lead_snp,
  stim_eQTL = data_all[which(data_all$lead_snp_source==paste0(this_condition,"_eQTL")),]$lead_snp
)
ld_pair1 <- ld_r2(genotype_all, snp1 = lead_snps$PBS_eQTL, snp2 = lead_snps$stim_eQTL)

# 1) helper: pick representative SNP in a set by smallest p
pick_rep <- function(snps) {
  d <- data_all[data_all$lead_snp %in% snps, c("lead_snp", "lead_snp_p"), drop = FALSE]
  d$lead_snp_p2 <- ifelse(is.na(d$lead_snp_p), Inf, d$lead_snp_p)
  d <- d[order(d$lead_snp_p2), , drop = FALSE]
  d$lead_snp[1]
}

# 2) define LD components
ld_thr = 0.6
snp1=data_all[1,]$lead_snp
snp2=data_all[2,]$lead_snp

components <- list()
if (ld_pair1$r2 >= ld_thr) {
  components <- list(LD1 = c(snp1, snp2))
} else {
  components <- list(LD1 = c(snp1), LD2 = c(snp2))
}

# 3) representative SNP per component
reps <- lapply(components, pick_rep)

# ----- 4. Assign primary class per gene × cytokine -----
PBS_sig <- data_all[which(data_all$lead_snp_source=="PBS_eQTL"),]$sig
cytokine_eQTL_sig <- data_all[which(data_all$lead_snp_source==paste0(this_condition,"_eQTL")),]$sig

primary_class = case_when(
  PBS_sig & cytokine_eQTL_sig  ~ "constitutive",
  !PBS_sig & cytokine_eQTL_sig ~ "emergent",
  PBS_sig & !cytokine_eQTL_sig ~ "vanishing",
  TRUE                          ~ NA  # neither
)

rep_snps <- unname(unlist(reps))
n_LDcomponent <- length(rep_snps)

# gene-level presence flags (based on your data_all rows)
has_PBS_eQTL      <- any(data_all$lead_snp_source == "PBS_eQTL" &
                           isTRUE(data_all$sig[data_all$lead_snp_source == "PBS_eQTL"]), na.rm = TRUE)
has_cyt_eQTL      <- any(data_all$lead_snp_source == paste0(this_condition, "_eQTL") &
                           isTRUE(data_all$sig[data_all$lead_snp_source == paste0(this_condition, "_eQTL")]), na.rm = TRUE)
has_cyt_reQTL     <- any(data_all$lead_snp_source == paste0(this_condition, "_reQTL") &
                           isTRUE(data_all$sig[data_all$lead_snp_source == paste0(this_condition, "_reQTL")]), na.rm = TRUE)

# ----- 5. Assign sub-class for constitutive only for each LD component
# ----- 6. ALSO Compile final classification table ----
rep_snps <- unname(unlist(reps))
out_rows <- lapply(rep_snps, function(this_snp) {
  
  # default
  sub_class <- NA_character_

  cytokine_beta <- modelstats.stim.eQTL[which(modelstats.stim.eQTL$snp==this_snp),]$beta
  cytokine_SE <- modelstats.stim.eQTL[which(modelstats.stim.eQTL$snp==this_snp),]$se
  cytokine_p <- modelstats.stim.eQTL[which(modelstats.stim.eQTL$snp==this_snp),]$p
  cytokine_baseline <- modelstats.stim.eQTL[which(modelstats.stim.eQTL$snp==this_snp),]$baseline
  PBS_beta <- modelstats.PBS.eQTL[which(modelstats.PBS.eQTL$snp==this_snp),]$beta
  PBS_SE <- modelstats.PBS.eQTL[which(modelstats.PBS.eQTL$snp==this_snp),]$se
  PBS_p <- modelstats.PBS.eQTL[which(modelstats.PBS.eQTL$snp==this_snp),]$p
  PBS_baseline <- modelstats.PBS.eQTL[which(modelstats.PBS.eQTL$snp==this_snp),]$baseline
  
  ratio = cytokine_beta / PBS_beta   # same-SNP betas
  Z     = (cytokine_beta - PBS_beta) / sqrt(cytokine_SE^2 + PBS_SE^2)
  
  sub_class = case_when(
    primary_class != "constitutive"    ~ "",
    abs(Z) < 1.96 ~ "_stable",   # change not stat. supported → stable
    # Among statistically supported changes, classify by direction and magnitude
    ratio < 0          ~ "_switched",    # direction flip
    ratio > 1    ~ "_amplified",
    ratio < 1    ~ "_attenuated",
    TRUE               ~ "_stable"        # supported but small change
  )
  
  if (!is.na(primary_class)) {
  primary_class <- paste0(primary_class, sub_class)
  }
  this_res <- data.frame(
    celltype = ct,
    cytokine = this_condition,
    gene = this_gene,
    gene_class = primary_class,
    has_PBS_eQTL = ifelse(has_PBS_eQTL, "yes", "no"),
    has_cytokine_eQTL = ifelse(has_cyt_eQTL, "yes", "no"),
    has_cytokine_reQTL = ifelse(has_cyt_reQTL, "yes", "no"),
    anchorSNP_betaComparison_ratio = ratio,
    anchorSNP_betaComparison_z = Z,
    anchorSNP_cytokine_beta = cytokine_beta,
    anchorSNP_cytokine_SE = cytokine_SE,
    anchorSNP_cytokine_p = cytokine_p,
    anchorSNP_cytokine_baseline = cytokine_baseline,
    anchorSNP_PBS_beta = PBS_beta,
    anchorSNP_PBS_SE = PBS_SE,
    anchorSNP_PBS_p = PBS_p,
    anchorSNP_PBS_baseline = PBS_baseline,
    anchorSNP_n_independent = n_LDcomponent,
    anchorSNP = this_snp,
    stringsAsFactors = FALSE
  )
})
out_df <- dplyr::bind_rows(out_rows)

# if this gene isn't interesting (not QTL anywhere), output row with NA
if (is.na(primary_class)){
  out_df <- data.frame(
    celltype = ct,
    cytokine = this_condition,
    gene = this_gene,
    gene_class = NA,
    has_PBS_eQTL = ifelse(has_PBS_eQTL, "yes", "no"),
    has_cytokine_eQTL = ifelse(has_cyt_eQTL, "yes", "no"),
    has_cytokine_reQTL = ifelse(has_cyt_reQTL, "yes", "no"),
    anchorSNP_betaComparison_ratio = NA,
    anchorSNP_betaComparison_z = NA,
    anchorSNP_cytokine_beta = NA,
    anchorSNP_cytokine_SE = NA,
    anchorSNP_cytokine_p = NA,
    anchorSNP_cytokine_baseline = NA,
    anchorSNP_PBS_beta = NA,
    anchorSNP_PBS_SE = NA,
    anchorSNP_PBS_p = NA,
    anchorSNP_PBS_baseline = NA,
    anchorSNP_n_independent = NA,
    anchorSNP = NA,
    stringsAsFactors = FALSE
  )
} 


# ----- 7. Write to file ----
outdir <- paste0(dir,"/reclassified/",this_condition)
outfname <- paste0("reclassified_",this_gene,".txt")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fwrite(out_df, file.path(outdir,outfname), quote=F, sep="\t")

message(sprintf(
  "ct=%s gene=%s cytokine=%s file=%s",
  ct, this_gene, this_condition, file.path(outdir,outfname)
))
