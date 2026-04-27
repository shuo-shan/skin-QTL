#!/usr/bin/env Rscript
# QTL reclassification per gene
# Usage: Rscript QTL_reclassify_perGene.R <celltype> <condition> <gene>
# Example: Rscript QTL_reclassify_perGene.R KRT IFNG ERAP2
#
# Logic:
#   1. Load lead SNPs from eigenMT tables for PBS eQTL, stim eQTL, stim reQTL
#   2. Compute pairwise LD (r2) among lead SNPs using your own genotype data
#   3. Cluster lead SNPs into LD groups (r2 >= 0.6 = same signal)
#   4. For each unique LD group, pull beta/SE from modelstats
#   5. Assign primary class, sub-class, canonical reQTL flag
#   6. Return one row per LD group (anchor) for this gene x condition

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
})

# ── CLI args ──────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript QTL_reclassify_perGene.R <celltype> <condition> <gene>\n")
  quit(save = "no", status = 1)
}
ct             <- args[[1]]
this_condition <- args[[2]]
this_gene      <- args[[3]]

# toy example of common LD SNP
# sig QTLs: ERAP2, not sig QTLs: IRF3
ct <- "KRT"
this_condition <- "IFNG"
this_gene <- "ERAP2"


# ── Config ────────────────────────────────────────────────────
LD_THRESHOLD     <- 0.6    # r2 >= this -> same signal
FDR_CUTOFF       <- 0.05
RATIO_AMPLIFIED  <- 1.5    # beta_cytokine / beta_PBS > this -> amplified
RATIO_ATTENUATED <- 0.67   # beta_cytokine / beta_PBS < this -> attenuated
Z_NOMINAL        <- 1.96   # |Z| for confirming direction change (switched)

# ── Paths ─────────────────────────────────────────────────────
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)

chunk_id_lookup <- fread(paste0(dir, "/data/gene_chunk_dict.txt"))
chunk_id        <- unique(chunk_id_lookup[gene == this_gene, chunk])
chunk_id        <- sprintf("%03d", chunk_id)

pair_file        <- paste0(dir, "/chunks/pairs_chunk_",    chunk_id, ".tsv")
geno_file        <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
modelstats_file  <- paste0(dir, "/results/result_",        chunk_id, ".tsv")
meta_file        <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"

eigenmt_dir      <- paste0(dir, "/eigenMT/results/")
qval_PBS_eQTL_file    <- paste0(eigenmt_dir, ct, "_PBS_eQTL.eigenMT.txt")
qval_stim_eQTL_file   <- paste0(eigenmt_dir, ct, "_", this_condition, "_eQTL.eigenMT.txt")
qval_stim_reQTL_file  <- paste0(eigenmt_dir, ct, "_", this_condition, "_reQTL.eigenMT.txt")

out_dir          <- paste0(dir, "/reclassified/",this_condition,"/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file         <- paste0(out_dir, this_gene, "_classified.tsv")

# ── Load data ─────────────────────────────────────────────────
message(sprintf("[ %s | %s | %s ] Loading data...", ct, this_condition, this_gene))

pairs <- fread(pair_file, header = TRUE) %>%
  filter(gene_name == this_gene) %>%
  mutate(SNP = str_to_lower(SNP_ID),
         key = paste0(gene_name, "_", SNP_ID))

snp_ids_for_gene <- unique(pairs$SNP_ID)

genotype_all <- fread(geno_file) %>%
  mutate(ID = str_to_lower(ID)) %>%
  filter(ID %in% snp_ids_for_gene)

modelstats <- fread(modelstats_file) %>%
  filter(gene == this_gene, snp %in% snp_ids_for_gene) %>%
  dplyr::select(snp, gene, dist, contains("_PBS"), contains(this_condition))

# ── Load eigenMT q-values ─────────────────────────────────────
load_qval <- function(f, gene) {
  if (!file.exists(f)) return(list(sig = FALSE, lead_snp = NA_character_,
                                   q = NA_real_,  pmin = NA_real_))
  dt <- fread(f) %>% filter(gene == !!gene)
  if (nrow(dt) == 0) return(list(sig = FALSE, lead_snp = NA_character_,
                                 q = NA_real_,  pmin = NA_real_))
  list(sig      = dt$q_gene[1] < FDR_CUTOFF,
       lead_snp = dt$lead_snp[1],
       q        = dt$q_gene[1],
       pmin     = dt$pmin[1])
}

qval_PBS_eQTL   <- load_qval(qval_PBS_eQTL_file,   this_gene)
qval_stim_eQTL  <- load_qval(qval_stim_eQTL_file,  this_gene)
qval_stim_reQTL <- load_qval(qval_stim_reQTL_file, this_gene)

is_PBS_eQTL  <- qval_PBS_eQTL$sig
is_stim_eQTL <- qval_stim_eQTL$sig
is_stim_reQTL<- qval_stim_reQTL$sig

message(sprintf("  PBS_eQTL=%s  %s_eQTL=%s  %s_reQTL=%s",
                is_PBS_eQTL, this_condition, is_stim_eQTL,
                this_condition, is_stim_reQTL))

# ── Collect lead SNPs ─────────────────────────────────────────
# Named list of lead SNPs; drop NAs (conditions not significant / missing)
lead_snps_raw <- list(
  PBS_eQTL          = qval_PBS_eQTL$lead_snp,
  stim_eQTL         = qval_stim_eQTL$lead_snp,
  stim_reQTL        = qval_stim_reQTL$lead_snp
)
# Rename stim keys to include the actual condition name for clarity
names(lead_snps_raw) <- c("PBS_eQTL",
                          paste0(this_condition, "_eQTL"),
                          paste0(this_condition, "_reQTL"))

lead_snps <- Filter(Negate(is.na), lead_snps_raw)

# ── LD function ───────────────────────────────────────────────
ld_r2 <- function(genotype_all, snp1, snp2) {
  if (is.na(snp1) || is.na(snp2)) return(list(r2=NA_real_, n_used=NA_integer_))
  if (snp1 == snp2)                return(list(r2=1.0,      n_used=NA_integer_))
  
  fixed_cols  <- intersect(c("CHROM","POS","ID","REF","ALT"), colnames(genotype_all))
  donor_cols  <- setdiff(colnames(genotype_all), fixed_cols)
  
  geno_matrix <- genotype_all %>%
    column_to_rownames("ID") %>%
    select(all_of(donor_cols))
  
  if (!snp1 %in% rownames(geno_matrix) || !snp2 %in% rownames(geno_matrix)) {
    message(sprintf("    WARNING: SNP not in genotype matrix: %s or %s", snp1, snp2))
    return(list(r2=NA_real_, n_used=NA_integer_))
  }
  
  g1 <- as.numeric(geno_matrix[snp1, ])
  g2 <- as.numeric(geno_matrix[snp2, ])
  ok <- is.finite(g1) & is.finite(g2)
  g1 <- g1[ok]; g2 <- g2[ok]
  n  <- length(g1)
  
  if (n < 10 || sd(g1) == 0 || sd(g2) == 0)
    return(list(r2=NA_real_, n_used=n))
  
  r <- cor(g1, g2)
  list(r2=r^2, n_used=n, snp1=snp1, snp2=snp2)
}

# ── Compute pairwise LD among all lead SNPs ───────────────────
unique_snps <- unique(unlist(lead_snps))  # deduplicated SNP IDs

message(sprintf("  Lead SNPs: %s", paste(names(lead_snps), unlist(lead_snps),
                                         sep="=", collapse="  ")))

# Build LD matrix for all unique lead SNPs
n_snps <- length(unique_snps)
ld_mat <- matrix(NA_real_, nrow=n_snps, ncol=n_snps,
                 dimnames=list(unique_snps, unique_snps))
diag(ld_mat) <- 1.0

if (n_snps > 1) {
  pairs_to_check <- combn(unique_snps, 2, simplify=FALSE)
  for (pr in pairs_to_check) {
    res <- ld_r2(genotype_all, pr[[1]], pr[[2]])
    ld_mat[pr[[1]], pr[[2]]] <- res$r2
    ld_mat[pr[[2]], pr[[1]]] <- res$r2
    message(sprintf("    LD r2(%s, %s) = %.3f",
                    pr[[1]], pr[[2]],
                    ifelse(is.na(res$r2), NA, round(res$r2, 3))))
  }
}

# ── Cluster lead SNPs into LD groups ─────────────────────────
# Two lead SNPs are in the same LD group if r2 >= LD_THRESHOLD.
# We use simple single-linkage clustering: if A~B and B~C, all three
# are grouped together. This is intentionally conservative — we only
# collapse SNPs we are confident are tagging the same signal.

get_ld_groups <- function(snp_names, ld_matrix, threshold) {
  n      <- length(snp_names)
  groups <- as.list(seq_len(n))          # start: each SNP its own group
  names(groups) <- snp_names
  
  # Build adjacency: which pairs exceed threshold?
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      r2_ij <- ld_matrix[snp_names[i], snp_names[j]]
      if (!is.na(r2_ij) && r2_ij >= threshold) {
        # Merge group of j into group of i
        gi <- groups[[snp_names[i]]]
        gj <- groups[[snp_names[j]]]
        merged <- union(gi, gj)
        for (k in snp_names) {
          if (groups[[k]] %in% c(gi, gj)) groups[[k]] <- merged
        }
      }
    }
  }
  
  # Convert to named factor: each SNP -> group_id (representative = first SNP)
  group_ids <- lapply(snp_names, function(s) {
    members <- snp_names[sapply(snp_names,
                                function(x) identical(groups[[x]], groups[[s]]))]
    members[1]   # representative = first alphabetically / by order
  })
  names(group_ids) <- snp_names
  group_ids
}

ld_groups <- get_ld_groups(unique_snps, ld_mat, LD_THRESHOLD)

# Map each lead role to its LD group representative
lead_snp_groups <- lapply(lead_snps, function(s) ld_groups[[s]])

# Invert: for each unique LD group representative, which roles does it cover?
unique_groups <- unique(unlist(lead_snp_groups))

anchor_roles <- lapply(unique_groups, function(rep_snp) {
  roles <- names(lead_snp_groups)[sapply(lead_snp_groups,
                                         function(g) identical(g, rep_snp))]
  list(anchor_snp   = rep_snp,
       anchor_sources = paste(sort(roles), collapse = "|"))
})

message(sprintf("  %d unique lead SNP(s) -> %d LD group(s)",
                length(unique_snps), length(unique_groups)))

# ── Pull modelstats for each anchor SNP ──────────────────────
get_snp_stats <- function(snp, modelstats, condition) {
  row <- modelstats[modelstats$snp == snp, ]
  if (nrow(row) == 0) return(NULL)
  row <- row[1, ]
  
  list(
    beta_PBS        = row$eQTL_beta_PBS,
    se_PBS          = row$eQTL_se_PBS,
    p_PBS           = row$eQTL_p_PBS,
    
    beta_cyt        = row[[paste0("eQTL_beta_", condition)]],
    se_cyt          = row[[paste0("eQTL_se_",   condition)]],
    p_cyt           = row[[paste0("eQTL_p_",    condition)]],
    
    dbeta_reQTL     = row[[paste0("reQTL_dbeta_", condition)]],
    se_reQTL        = row[[paste0("reQTL_se_",    condition)]],
    p_reQTL         = row[[paste0("reQTL_p_",     condition)]]
  )
}

# ── Classify each LD group ────────────────────────────────────
classify_anchor <- function(anchor_info, stats,
                            is_PBS_eQTL, is_stim_eQTL, is_stim_reQTL,
                            condition) {
  if (is.null(stats)) return(NULL)
  
  beta_pbs <- stats$beta_PBS
  se_pbs   <- stats$se_PBS
  beta_cyt <- stats$beta_cyt
  se_cyt   <- stats$se_cyt
  
  # Primary class — based on gene-level significance flags
  primary_class <- dplyr::case_when(
    is_PBS_eQTL &  is_stim_eQTL ~ "constitutive_eQTL",
    !is_PBS_eQTL &  is_stim_eQTL ~ "emergent_eQTL",
    is_PBS_eQTL & !is_stim_eQTL ~ "vanishing_eQTL",
    TRUE                          ~ NA_character_
  )
  
  # Beta ratio and Z-score (only meaningful for constitutive)
  beta_ratio    <- NA_real_
  Z_beta_change <- NA_real_
  
  if (!is.na(beta_pbs) && is.finite(beta_pbs) && beta_pbs != 0 &&
      !is.na(beta_cyt) && is.finite(beta_cyt) &&
      !is.na(se_pbs)   && !is.na(se_cyt)) {
    beta_ratio    <- beta_cyt / beta_pbs
    Z_beta_change <- (beta_cyt - beta_pbs) / sqrt(se_cyt^2 + se_pbs^2)
  }
  
  # Sub-class (constitutive only)
  sub_class <- NA_character_
  if (!is.na(primary_class) && primary_class == "constitutive_eQTL" &&
      !is.na(beta_ratio)) {
    sub_class <- dplyr::case_when(
      beta_ratio < 0  & abs(Z_beta_change) >= Z_NOMINAL ~ "switched",
      beta_ratio > RATIO_AMPLIFIED                       ~ "amplified",
      beta_ratio >= 0 & beta_ratio < RATIO_ATTENUATED    ~ "attenuated",
      beta_ratio >= RATIO_ATTENUATED &
        beta_ratio <= RATIO_AMPLIFIED                    ~ "stable",
      TRUE                                               ~ NA_character_
    )
  }
  
  full_class <- dplyr::case_when(
    !is.na(sub_class)    ~ paste0("constitutive_", sub_class),
    !is.na(primary_class)~ primary_class,
    TRUE                 ~ NA_character_
  )
  
  data.table(
    gene                 = this_gene,
    celltype             = ct,
    cytokine             = condition,
    anchor_snp           = anchor_info$anchor_snp,
    anchor_sources       = anchor_info$anchor_sources,
    
    # Gene-level significance
    PBS_eQTL_sig         = is_PBS_eQTL,
    cytokine_eQTL_sig    = is_stim_eQTL,
    canonical_reQTL_sig  = is_stim_reQTL,
    
    # Q-values
    PBS_eQTL_q           = qval_PBS_eQTL$q,
    cytokine_eQTL_q      = qval_stim_eQTL$q,
    canonical_reQTL_q    = qval_stim_reQTL$q,
    
    # Classification
    primary_class        = primary_class,
    sub_class            = sub_class,
    full_class           = full_class,
    
    # Effect sizes at anchor SNP
    beta_PBS             = beta_pbs,
    se_PBS               = se_pbs,
    p_PBS                = stats$p_PBS,
    beta_cytokine        = beta_cyt,
    se_cytokine          = se_cyt,
    p_cytokine           = stats$p_cyt,
    beta_ratio           = beta_ratio,
    Z_beta_change        = Z_beta_change,
    
    # Interaction term at anchor SNP
    dbeta_reQTL          = stats$dbeta_reQTL,
    se_reQTL             = stats$se_reQTL,
    p_reQTL              = stats$p_reQTL,
    
    # LD metadata
    LD_threshold_used    = LD_THRESHOLD,
    n_unique_lead_snps   = length(unique_snps),
    n_ld_groups          = length(unique_groups)
  )
}

results <- rbindlist(lapply(anchor_roles, function(anc) {
  stats <- get_snp_stats(anc$anchor_snp, modelstats, this_condition)
  classify_anchor(anc, stats,
                  is_PBS_eQTL, is_stim_eQTL, is_stim_reQTL,
                  this_condition)
}), fill = TRUE)

# ── Add LD pairwise summary columns ──────────────────────────
# Attach the raw r2 values between the three lead SNPs for transparency
pbs_snp  <- lead_snps[[paste0("PBS_eQTL")]]
cyt_snp  <- lead_snps[[paste0(this_condition, "_eQTL")]]
rqtl_snp <- lead_snps[[paste0(this_condition, "_reQTL")]]

safe_r2 <- function(s1, s2) {
  if (is.null(s1) || is.null(s2) || is.na(s1) || is.na(s2)) return(NA_real_)
  ld_mat[s1, s2]
}

results[, `:=`(
  PBS_lead_snp           = ifelse(is.null(pbs_snp),  NA_character_, pbs_snp),
  cytokine_eQTL_lead_snp = ifelse(is.null(cyt_snp),  NA_character_, cyt_snp),
  cytokine_reQTL_lead_snp= ifelse(is.null(rqtl_snp), NA_character_, rqtl_snp),
  r2_PBS_vs_cytEQTL      = safe_r2(pbs_snp,  cyt_snp),
  r2_PBS_vs_reQTL        = safe_r2(pbs_snp,  rqtl_snp),
  r2_cytEQTL_vs_reQTL    = safe_r2(cyt_snp,  rqtl_snp)
)]

# ── Write output ──────────────────────────────────────────────
fwrite(results, out_file, sep = "\t")
message(sprintf("  Written: %s  (%d row(s))", out_file, nrow(results)))
print(results[, .(anchor_snp, anchor_sources, full_class,
                  canonical_reQTL_sig, beta_ratio,
                  r2_PBS_vs_cytEQTL, r2_PBS_vs_reQTL, r2_cytEQTL_vs_reQTL)])