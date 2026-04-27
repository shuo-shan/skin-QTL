#!/usr/bin/env Rscript
# script for running adaptive permutation for reQTLs per gene

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
  library(coloc) # coloc 5.2.3
  library(susieR)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

#### Input ####
# ------------------------------#
# Main: parse CLI args ####
# ------------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("
Usage:
  Rscript step7_susie_perGene.R <ct> <gene>

Example:
  Rscript step7_susie_perGene.R MEL ERAP2
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
this_gene       <- args[[2]] # ERAP2

# # toy example for debugging
# ct         <- "MEL"
# this_gene <- "FAM210A"

# set-up global variables
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
chunk_id_lookup <- data.table::fread(paste0(dir,"/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[which(chunk_id_lookup$gene==this_gene),]$chunk)
chunk_id <- sprintf("%03d", chunk_id)
pair_file <- paste0(dir,"/data/chunk/pair_chunk_",chunk_id,".txt")
geno_file  <- paste0(dir,"/data/chunk/genotype_pair_chunk_",chunk_id,".txt")
modelstats_file <- paste0(dir,"/permutation/model_stats/model_stats_",chunk_id,".txt")
meta_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
rm(chunk_id_lookup)

# load SNP:gene pairs (global)
pairs <- fread(pair_file, header=F) %>% 
  magrittr::set_colnames(c("SNP","gene")) %>%
  dplyr::filter(gene == this_gene)
pairs$key <- paste0(pairs$gene,"_",pairs$SNP)

# load genotype of all SNPs in chunk (global)
genotype_all <- fread(geno_file) %>% dplyr::filter(ID %in% pairs$SNP)

# load metadata for all samples (global)
meta_all <- readr::read_tsv(meta_file, show_col_types = FALSE)

# ---- initialize master long table with correct columns/types ----
susie_master_long <- data.frame(
    snp        = character(),
    gene       = character(),
    dist       = numeric(),
    condition  = character(),
    QTLtype    = character(),
    beta       = numeric(),
    se         = numeric(),
    pval       = numeric(),
    z          = numeric(),
    pip        = numeric(),
    cs90_index = integer(),
    cs90_size  = integer(),
    is_in_cs90 = logical(),
    cs95_index = integer(),
    cs95_size  = integer(),
    is_in_cs95 = logical(),
    stringsAsFactors = FALSE
  )

#### FUNCTION: run SuSiE for 1 trait
run_susie_one_trait <- function(ct, this_gene, condition, QTLtype,
                                dir, modelstats_file, meta_all,
                                pairs, genotype_all,
                                L = 10, coverage_vec = c(0.90, 0.95),
                                min_abs_corr_cs = 0.60,
                                min_abs_corr_fit = 0.10) {
  
  message("=== running susie for ", ct, " ", condition, " ", QTLtype, ", gene: ", this_gene, " ===")
  
  # ---- Load modeling stats for the QTL of interest of chunk ----
  modelstats_all <- data.table::fread(modelstats_file)
  
  fixed_col <- c("snp","gene","dist")
  stats_col <- setdiff(colnames(modelstats_all), fixed_col)
  
  pattern <- paste0("^", QTLtype, ".*", condition, "$")
  selected_stats_col <- grep(pattern, stats_col, value = TRUE)
  
  if (length(selected_stats_col) == 0) {
    warning("No model stats columns matched pattern: ", pattern)
    return(NULL)
  }
  
  modelstats_subset <- modelstats_all[ , c(fixed_col, selected_stats_col), with = FALSE] %>%
    dplyr::mutate(key = paste0(gene,"_",snp)) %>%
    dplyr::filter(key %in% pairs$key)
  
  # Fix p=0 for display / stability of -log10(p) later (does not affect z, unless you use p to recompute)
  pval_col <- grep("_p_", colnames(modelstats_subset), value = TRUE)
  if (length(pval_col) == 1) {
    p <- modelstats_subset[[pval_col]]
    min_nonzero <- min(p[p > 0], na.rm = TRUE)
    if (is.finite(min_nonzero)) {
      modelstats_subset[[pval_col]] <- ifelse(p == 0, min_nonzero / 2, p)
    }
  }
  rm(modelstats_all)
  
  # ---- Load metadata to define LD donor set ----
  if (QTLtype == "eQTL") {
    kept_condition <- condition
  } else {
    kept_condition <- unique(c("PBS", condition))
  }
  
  meta <- meta_all %>%
    dplyr::filter(celltype == ct, condition %in% kept_condition) %>%
    dplyr::select(sample, donor, condition)
  
  donors_ld <- unique(meta$donor)
  
  # ---- build input data for susie ----
  beta_col <- grep("beta", colnames(modelstats_subset), value = TRUE)
  se_col   <- grep("_se_", colnames(modelstats_subset), value = TRUE)
  pval_col <- grep("_p_",  colnames(modelstats_subset), value = TRUE)
  
  if (length(beta_col) != 1 || length(se_col) != 1 || length(pval_col) != 1) {
    warning("Could not uniquely identify beta/se/p columns for: ", QTLtype, " ", condition)
    return(NULL)
  }
  
  qtl_df <- modelstats_subset %>%
    dplyr::transmute(
      snp = snp,
      gene = gene,
      dist = dist,
      beta = .data[[beta_col]],
      se   = .data[[se_col]],
      pval = .data[[pval_col]]
    ) %>%
    dplyr::mutate(z = beta / se) %>%
    dplyr::filter(is.finite(z))
  
  # Keep SNPs present in both tables
  common_snps <- intersect(qtl_df$snp, genotype_all$ID)
  if (length(common_snps) < 10) {
    warning("Too few SNPs after intersecting genotype and modelstats: ", length(common_snps))
    return(NULL)
  }
  
  qtl_df <- qtl_df %>% dplyr::filter(snp %in% common_snps) %>% dplyr::arrange(match(snp, common_snps))
  geno_df <- genotype_all %>% dplyr::filter(ID %in% common_snps) %>% dplyr::arrange(match(ID, common_snps))
  stopifnot(all(qtl_df$snp == geno_df$ID))
  
  # Subset genotype to LD donors and compute LD matrix
  donors_in_geno <- intersect(donors_ld, colnames(geno_df))
  if (length(donors_in_geno) < 5) stop("Too few donors present in genotype table for LD: ", length(donors_in_geno))
  
  G <- as.matrix(geno_df[, donors_in_geno, with = FALSE])
  storage.mode(G) <- "numeric"
  
  # drop SNPs with no variance
  snp_sd <- apply(G, 1, sd, na.rm = TRUE)
  keep_var <- is.finite(snp_sd) & snp_sd > 0
  if (sum(keep_var) < 10) {
    warning("Too few variable SNPs after donor subsetting: ", sum(keep_var))
    return(NULL)
  }
  
  G <- G[keep_var, , drop = FALSE]
  qtl_df <- qtl_df[keep_var, , drop = FALSE]
  
  # LD across SNPs
  R <- stats::cor(t(G), use = "pairwise.complete.obs")
  eig_min <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(eig_min) || eig_min < -1e-8) {
    warning("LD matrix has negative eigenvalues (min eig = ", signif(eig_min, 3), "). making PD.")
    # Even small negative eigenvalues from numerical noise can cause issues.
    R <- corpcor::make.positive.definite(R)
  }
  
  n_eff <- ncol(G)
  z <- qtl_df$z
  
  if (!all(is.finite(z))) stop("Non-finite z-scores present.")
  if (any(diag(R) < 0.99)) warning("LD diagonal not ~1; check genotype coding/missingness.")
  
  # ---- run SuSiE ----
  start_time <- Sys.time()
  fit <- susieR::susie_rss(
    z = z,
    R = R,
    n = n_eff,
    L = L,
    estimate_residual_variance = FALSE,
    min_abs_corr = min_abs_corr_fit
  )
  runtime_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  message(sprintf("Runtime: %.2f seconds", runtime_sec))
  
  pip_df <- data.frame(
    snp  = qtl_df$snp,
    gene = this_gene,
    dist = qtl_df$dist,
    condition = condition,
    QTLtype = QTLtype,
    beta = qtl_df$beta,
    se   = qtl_df$se,
    pval = qtl_df$pval,
    z    = qtl_df$z,
    pip  = fit$pip,
    stringsAsFactors = FALSE
  )
  
  # ---- Credible sets helper ----
  build_cs_membership <- function(fit, R, snp_ids, coverage, min_abs_corr = 0.6) {
    cs_obj <- susieR::susie_get_cs(
      fit,
      Xcorr = R,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      dedup = TRUE
    )
    n_cs <- length(cs_obj$cs)
    
    if (n_cs == 0) {
      return(tibble::tibble(snp = character(0), cs_index = integer(0), cs_size = integer(0)))
    }
    
    cs_summary <- tibble::tibble(
      cs_index = seq_len(n_cs),
      cs_size  = vapply(cs_obj$cs, length, integer(1)),
      component_lbf = if (!is.null(fit$lbf)) fit$lbf[seq_len(n_cs)] else NA_real_
    )
    
    mem_long <- tibble::tibble(
      cs_index = rep(seq_len(n_cs), times = vapply(cs_obj$cs, length, integer(1))),
      snp_idx  = unlist(cs_obj$cs, use.names = FALSE)
    ) %>%
      dplyr::mutate(snp = snp_ids[snp_idx]) %>%
      dplyr::left_join(cs_summary, by = "cs_index")
    
    # tie-break: higher lbf, then smaller cs, then earlier cs
    mem_best <- mem_long %>%
      dplyr::arrange(snp, dplyr::desc(component_lbf), cs_size, cs_index) %>%
      dplyr::group_by(snp) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(snp, cs_index, cs_size)
    
    mem_best
  }
  
  # ---- build CS90 + CS95 and join ----
  cs90_mem <- build_cs_membership(fit, R, qtl_df$snp, coverage = 0.90, min_abs_corr = min_abs_corr_cs) %>%
    dplyr::rename(cs90_index = cs_index, cs90_size = cs_size) %>%
    dplyr::mutate(is_in_cs90 = TRUE)
  
  cs95_mem <- build_cs_membership(fit, R, qtl_df$snp, coverage = 0.95, min_abs_corr = min_abs_corr_cs) %>%
    dplyr::rename(cs95_index = cs_index, cs95_size = cs_size) %>%
    dplyr::mutate(is_in_cs95 = TRUE)
  
  # ---- compile result table ----
  pip_df_aug <- pip_df %>%
    dplyr::left_join(cs90_mem, by = "snp") %>%
    dplyr::left_join(cs95_mem, by = "snp") %>%
    dplyr::mutate(
      is_in_cs90 = ifelse(is.na(is_in_cs90), FALSE, is_in_cs90),
      is_in_cs95 = ifelse(is.na(is_in_cs95), FALSE, is_in_cs95),
      cs90_index = ifelse(is.na(cs90_index), NA_integer_, as.integer(cs90_index)),
      cs95_index = ifelse(is.na(cs95_index), NA_integer_, as.integer(cs95_index)),
      cs90_size  = ifelse(is.na(cs90_size),  NA_integer_, as.integer(cs90_size)),
      cs95_size  = ifelse(is.na(cs95_size),  NA_integer_, as.integer(cs95_size))
    )
  
  # Ensure expected columns exist (in case cs tables are empty)
  needed <- c("snp","gene","dist","condition","QTLtype","beta","se","pval","z","pip",
              "cs90_index","cs90_size","is_in_cs90","cs95_index","cs95_size","is_in_cs95")
  for (nm in setdiff(needed, colnames(pip_df_aug))) pip_df_aug[[nm]] <- NA
  
  pip_df_aug <- pip_df_aug[, needed]
  
  pip_df_aug
}

#### MAIN: loop over different traits: condition X QTLtype ####
trait_grid <- tibble::tribble(
  ~QTLtype, ~condition,
  "eQTL",  "PBS",
  "eQTL",  "IFNG",
  "eQTL",  "IFNB",
  "eQTL",  "TNF",
  "reQTL", "IFNG",
  "reQTL", "IFNB",
  "reQTL", "TNF"
)

for (i in seq_len(nrow(trait_grid))) {
  QTLtype_i   <- trait_grid$QTLtype[[i]]
  condition_i <- trait_grid$condition[[i]]
  
  res <- run_susie_one_trait(
    ct = ct,
    this_gene = this_gene,
    condition = condition_i,
    QTLtype = QTLtype_i,
    dir = dir,
    modelstats_file = modelstats_file,
    meta_all = meta_all,
    pairs = pairs,
    genotype_all = genotype_all,
    L = 10
  ) %>% arrange(desc(is_in_cs90), desc(is_in_cs95), desc(pip))
  
  if (!is.null(res)) {
    output_file <- paste0("susie_",this_gene,".txt")
    data.table::fwrite(res, file = file.path(dir,"susie",condition_i,QTLtype_i,output_file),
                       quote=F, sep = "\t")
    message("Wrote: ", file.path(dir,"susie",condition_i,QTLtype_i,output_file))
    
    susie_master_long <- rbind(susie_master_long, res)
  }
}

output_file <- paste0("susie_",this_gene,".txt")
data.table::fwrite(susie_master_long, file = file.path(dir,"susie","output_long",output_file),
                   quote=F, sep="\t")
message("Wrote: ", file.path(dir,"susie","output_long",output_file))

##### QC plots #####
#---------------------------
# Helpers
#---------------------------

add_trait_label <- function(df) {
  df %>%
    mutate(trait = paste0(QTLtype, ":", condition))
}

safe_log10p <- function(p) {
  # p should already have p==0 fixed upstream, but keep this safe
  p <- ifelse(is.na(p), NA_real_, p)
  p <- pmin(pmax(p, 1e-300), 1)  # clamp
  -log10(p)
}

#---------------------------
# 1) PIP vs -log10(p) per trait
plot_pip_vs_p <- function(df_gene, pip_floor = 1e-6) {
  df_gene %>%
    add_trait_label() %>%
    mutate(
      log10p = safe_log10p(pval),
      pip_plot = pmax(pip, pip_floor)
    ) %>%
    ggplot(aes(x = log10p, y = pip_plot)) +
    geom_point(aes(shape = is_in_cs90, alpha = is_in_cs95), size = 1.6) +
    scale_y_continuous(trans = "log10", labels = label_number()) +
    facet_wrap(~ trait, scales = "free_x") +
    labs(
      title = unique(df_gene$gene),
      subtitle = "PIP vs -log10(p) (shape = CS90 membership, alpha = CS95 membership)",
      x = "-log10(p)",
      y = "PIP (log scale)"
    ) +
    theme_bw(base_size = 12)
}

#---------------------------
# 2) CS size QC: distribution of CS sizes per trait (for SNPs in CS)
#   This shows how big your CS are. If you see a ton of huge CS sizes,
#   it’s a sign many components are “diffuse / uninformative”.
plot_cs_size <- function(df_gene, which_cs = c("cs90", "cs95")) {
  which_cs <- match.arg(which_cs)
  
  if (which_cs == "cs90") {
    tmp <- df_gene %>%
      add_trait_label() %>%
      filter(is_in_cs90) %>%
      distinct(trait, cs90_index, cs90_size) %>%
      rename(cs_index = cs90_index, cs_size = cs90_size)
    ttl <- "CS90 size distribution (unique CS per trait)"
  } else {
    tmp <- df_gene %>%
      add_trait_label() %>%
      filter(is_in_cs95) %>%
      distinct(trait, cs95_index, cs95_size) %>%
      rename(cs_index = cs95_index, cs_size = cs95_size)
    ttl <- "CS95 size distribution (unique CS per trait)"
  }
  
  ggplot(tmp, aes(x = trait, y = cs_size)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.6, size = 1.6) +
    scale_y_continuous(trans = "log10") +
    coord_flip() +
    labs(
      title = unique(df_gene$gene),
      subtitle = ttl,
      x = NULL,
      y = "CS size (log10 scale)"
    ) +
    theme_bw(base_size = 12)
}

#---------------------------
# 3) PIP concordance between two traits (scatter)
#   Great for: PBS eQTL vs IFNG eQTL; IFNG eQTL vs IFNG reQTL; etc.
plot_pip_concordance <- function(df_gene, trait_a, trait_b, pip_floor = 1e-6) {
  wide <- df_gene %>%
    add_trait_label() %>%
    filter(trait %in% c(trait_a, trait_b)) %>%
    select(snp, trait, pip, is_in_cs95) %>%
    mutate(pip = pmax(pip, pip_floor)) %>%
    tidyr::pivot_wider(
      names_from = trait,
      values_from = c(pip, is_in_cs95),
      values_fill = list(pip = 0, is_in_cs95 = FALSE)
    )
  
  pip_a  <- paste0("pip_", trait_a)
  pip_b  <- paste0("pip_", trait_b)
  cs95_a <- paste0("is_in_cs95_", trait_a)
  cs95_b <- paste0("is_in_cs95_", trait_b)
  
  wide <- wide %>%
    mutate(
      cs_tier = dplyr::case_when(
        .data[[cs95_a]] & .data[[cs95_b]] ~ "CS95 both",
        .data[[cs95_a]] & !.data[[cs95_b]] ~ paste0("CS95 ", trait_a, " only"),
        !.data[[cs95_a]] & .data[[cs95_b]] ~ paste0("CS95 ", trait_b, " only"),
        TRUE ~ "Neither"
      ),
      cs_tier = factor(cs_tier, levels = c(
        "Neither",
        paste0("CS95 ", trait_a, " only"),
        paste0("CS95 ", trait_b, " only"),
        "CS95 both"
      ))
    )
  
  ggplot(wide, aes(x = .data[[pip_a]], y = .data[[pip_b]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(aes(color = cs_tier), alpha = 0.75, size = 1.8) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(
      values = setNames(
        c("grey80", "orange2", "orchid3", "deepskyblue2"),
        c(
          "Neither",
          paste0("CS95 ", trait_a, " only"),
          paste0("CS95 ", trait_b, " only"),
          "CS95 both"
        )
      )
    ) +
    labs(
      title = unique(df_gene$gene),
      subtitle = paste0("PIP concordance: ", trait_a, " vs ", trait_b, " (color = CS95 tier)"),
      x = paste0("PIP (", trait_a, ")"),
      y = paste0("PIP (", trait_b, ")"),
      color = "CS95 tier"
    ) +
    theme_bw(base_size = 12)
}



#---------------------------
# 4) Optional locus track: PIP along distance (per trait)
#   Helpful to see if fine-mapping is concentrated in a small LD block.
plot_locus_track <- function(df_gene, this_trait,
                             col_none = "grey75",
                             col_cs95 = "steelblue3",
                             col_cs90 = "deepskyblue2") {
  
  wide <- df_gene %>%
    add_trait_label() %>%
    filter(trait==this_trait) %>%
    select(snp, gene, dist, pval, trait, pip, is_in_cs90, is_in_cs95) %>%
    mutate(
      # tier priority: CS90 > CS95 > none
      tier = dplyr::case_when(
        is_in_cs90 ~ "CS90",
        is_in_cs95 ~ "CS95",
        TRUE       ~ "none"
      ),
      tier = factor(tier, levels = c("none", "CS95", "CS90")),
      log10p = -log10(pval)
    )
  
  out_prefix <- paste(unique(wide$gene), this_trait)
  
  # choose lead SNP (smallest p); if ties, label all ties
  lead <- wide %>% filter(pval == min(pval, na.rm = TRUE))
  
  # shared theme bits
  tier_scale <- scale_color_manual(
    values = c("none" = col_none, "CS95" = col_cs95, "CS90" = col_cs90),
    drop = FALSE
  )
    
  plot_model_pval <- ggplot(wide, aes(x = dist, y = log10p, color = tier)) +
    geom_point(alpha = 0.85, size = 1.6) +
    geom_label_repel(
      data = lead,
      aes(label = snp),
      size = 2.6,
      max.overlaps = 40,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0
    ) +
    tier_scale +
    labs(
      title = paste0("Model p-values: ", out_prefix),
      x = "SNP distance to gene TSS (bp)",
      y = "-log10(p)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  
  plot_model_pip <- ggplot(wide, aes(x = dist, y = pip, color = tier)) +
    geom_point(alpha = 0.85, size = 1.6) +
    tier_scale +
    labs(
      title = paste0("SuSiE PIP: ", out_prefix),
      x = "SNP distance to gene TSS (bp)",
      y = "PIP"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  # 2x1 stacked layout
  plot_model_pval / plot_model_pip + plot_layout(heights = c(1, 1))
  
  
}

df_gene <- susie_master_long %>% filter(gene == this_gene)

plots <- list(
  p1 = plot_pip_vs_p(df_gene),
  p2 = plot_cs_size(df_gene, which_cs = "cs90"),
  p3 = plot_cs_size(df_gene, which_cs = "cs95"),
  p4 = plot_pip_concordance(df_gene, "eQTL:PBS", "eQTL:IFNG"),
  p5 = plot_pip_concordance(df_gene, "eQTL:PBS", "eQTL:IFNB"),
  p6 = plot_pip_concordance(df_gene, "eQTL:PBS", "eQTL:TNF"),
  p7 = plot_pip_concordance(df_gene, "eQTL:IFNG", "reQTL:IFNG"),
  p8 = plot_pip_concordance(df_gene, "eQTL:IFNB", "reQTL:IFNB"),
  p9 = plot_pip_concordance(df_gene, "eQTL:TNF", "reQTL:TNF"),
  p10 = plot_locus_track(df_gene, "eQTL:PBS"),
  p11 = plot_locus_track(df_gene, "eQTL:IFNG"),
  p12 = plot_locus_track(df_gene, "eQTL:IFNB"),
  p13 = plot_locus_track(df_gene, "eQTL:TNF"),
  p14 = plot_locus_track(df_gene, "reQTL:IFNG"),
  p15 = plot_locus_track(df_gene, "reQTL:IFNB"),
  p16 = plot_locus_track(df_gene, "reQTL:TNF")
)

# ---- output path ----
output_file <- paste0("susie_", this_gene, ".pdf")
pdf_file <- file.path(dir, "susie", "QC", output_file)

# ---- write multipage pdf ----
pdf(pdf_file, width = 9.3, height = 8, onefile = TRUE)

for (nm in names(plots)) {
  p <- plots[[nm]]
  if (!is.null(p)) print(p)
}

dev.off()
message("Wrote: ", pdf_file)

