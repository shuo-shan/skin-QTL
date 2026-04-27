#!/usr/bin/env Rscript
# script for running coloc on vitiligo per gene

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
  Rscript step8_coloc_perGene_vitiligo.R <ct> <gene>

Example:
  Rscript step8_coloc_perGene_vitiligo.R MEL ERAP2
", "\n")
  quit(save = "no", status = 1)
}

ct         <- args[[1]]
this_gene       <- args[[2]] # ERAP2

# # toy example for debugging
# ct         <- "MEL"
# this_gene <- "AC069234.2"

#### ---- set-up ---- ####
# set-up global variables of my data
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
chunk_id_lookup <- data.table::fread(paste0(dir,"/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[which(chunk_id_lookup$gene==this_gene),]$chunk)
chunk_id <- sprintf("%03d", chunk_id)
pair_file <- paste0(dir,"/data/chunk/pair_chunk_",chunk_id,".txt")
geno_file  <- paste0(dir,"/data/chunk/genotype_pair_chunk_",chunk_id,".txt")
vst_file   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
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

# get chromosome
chr <- unique(genotype_all$CHROM)

#### ---- function for locus plot ---- ####
add_trait_label <- function(df) {
  df %>%
    mutate(trait = paste0(QTLtype, ":", condition))
}

plot_locus_track <- function(harmonized, this_trait) {
  
  wide <- harmonized %>%
    dplyr::select(SNP, GENE, DIST, P, P_QTL) %>%
    dplyr::mutate(
      log10p.gwas = -log10(P),
      log10p.qtl  = -log10(P_QTL)
    )
  
  # ---- guards ----
  if (nrow(wide) == 0) {
    message("[plot_locus_track] No rows for trait = ", this_trait)
    return(NULL)
  }
  
  if (all(is.na(wide$P_QTL))) {
    message("[plot_locus_track] P_QTL is all NA for trait = ", this_trait,
            " (n=", nrow(wide), "). Skipping lead label.")
    lead <- wide[0, ]
  } else {
    lead_p <- min(wide$P_QTL, na.rm = TRUE)
    lead <- wide %>% dplyr::filter(P_QTL == lead_p)
  }
  
  out_prefix <- paste(unique(wide$GENE), this_trait)
  
  plot_pval_qtl <- ggplot(wide, aes(x = DIST, y = log10p.qtl)) +
    geom_point(alpha = 0.85, size = 1.6) +
    geom_label_repel(
      data = lead,
      aes(label = SNP),
      size = 2.6,
      max.overlaps = 40,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0
    ) +
    labs(
      title = paste0("QTL model p-values: ", out_prefix),
      x = "SNP distance to gene TSS (bp)",
      y = "-log10(p.qtl)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  
  plot_pval_gwas <- ggplot(wide, aes(x = DIST, y = log10p.gwas)) +
    geom_point(alpha = 0.85, size = 1.6) +
    labs(
      title = paste0("Vitiligo GWAS locus: ", unique(wide$GENE)),
      x = "SNP distance to gene TSS (bp)",
      y = "-log10(p.gwas)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  plot_pval_qtl / plot_pval_gwas + plot_layout(heights = c(1, 1))
}


#### ---- function for running coloc ---- ####
run_coloc_one_trait <- function(
    ct,
    this_gene,
    QTLtype,
    condition,
    chr,
    dir,
    pairs,
    genotype_all,
    meta_all,
    modelstats_file,
    vst_file,
    # GWAS inputs
    gwas_dir = "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/GWAS_vitiligo/JinY_27723757_GCST004785",
    n_cases = 2853,
    n_controls = 37405
) {
  #### ---- load QTL metadata and VST file ---- ####
  # Load metadata
  if (QTLtype == "eQTL") {
    kept_condition <- condition
  } else {
    kept_condition <- unique(c("PBS", condition))
  }
  
  meta <- meta_all %>%
    dplyr::filter(celltype == ct, condition %in% kept_condition) %>%
    dplyr::select(sample, donor, condition)
  
  # compute sdY for the exact phenotype used in the QTL model
  y <- data.table::fread(vst_file) %>%
    dplyr::filter(gene == this_gene) %>%
    dplyr::select(all_of(meta$sample)) %>%
    t() %>% as.data.frame() %>% dplyr::pull(V1)
  
  sdY_qtl <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sdY_qtl) || sdY_qtl <= 0) {
    message(sprintf(
      "[SKIP] %s %s %s: constant VST (sdY = %s)",
      ct, condition, this_gene, as.character(sdY_qtl)
    ))
    return(NULL)
  }
 
  #### ---- load GWAS data ---- ####
  gwas_file <- file.path(gwas_dir, paste0("GWAS123", chr, "cmh.txt.gz"))
  if (!file.exists(gwas_file)) stop("GWAS file not found: ", gwas_file)
  
  # GWAS: MUST fill these in from the GWAS paper / GCST metadata
  N_gwas <- n_cases + n_controls
  s_gwas <- n_cases / N_gwas
  
  gwas <- readr::read_tsv(
    gwas_file,
    show_col_types = FALSE,
    col_types = cols(
      CHR   = col_integer(),
      SNP   = col_character(),
      BP    = col_integer(),
      A1    = col_character(),
      MAF   = col_double(),
      A2    = col_character(),
      CHISQ = col_double(),
      P     = col_double(),
      ORX   = col_double(),
      SE    = col_double(),
      L95   = col_double(),
      U95   = col_double()
    )
  ) %>%
    dplyr::mutate(
      SNP = stringr::str_to_lower(SNP),  # "RS..." -> "rs..."
      A1  = stringr::str_to_upper(A1),
      A2  = stringr::str_to_upper(A2)
    )
  
  #### ---- load genotype data ---- ####
  # calculate MAF of my data
  fixed_col <- c("CHROM","POS","ID","REF","ALT")
  donor_col <- setdiff(colnames(genotype_all), fixed_col)
  
  geno_mat <- genotype_all %>%
    dplyr::select(all_of(donor_col)) %>%
    dplyr::mutate(dplyr::across(everything(), as.numeric)) %>%
    as.matrix()
  rownames(geno_mat) <- genotype_all$ID %>% stringr::str_to_lower()
  
  # ALT allele frequency = mean(dosage) / 2
  alt_af <- rowMeans(geno_mat, na.rm = TRUE) / 2
  
  # MAF = min(AF, 1-AF)
  maf <- pmin(alt_af, 1 - alt_af)
  
  # Also useful: number of non-missing genotypes per SNP
  n_called <- rowSums(!is.na(geno_mat))
  
  maf_tbl <- tibble::tibble(
    ID = names(alt_af) %>% stringr::str_to_lower(),
    AF_ALT = unname(alt_af),
    MAF_QTL = unname(maf),
    N_CALLED = unname(n_called)
  )
  
  # Append to variant_tbl (make sure IDs match casing)
  # genotype_all has: CHROM POS ID REF ALT and dosage columns (0/1/2)
  variant_tbl <- genotype_all %>%
    dplyr::select(POS, ID, REF, ALT) %>%
    dplyr::mutate(ID = stringr::str_to_lower(ID)) %>%   # <-- move lowercase BEFORE join (important)
    dplyr::left_join(maf_tbl, by = "ID")
  
  #### ---- load QTL stats and metadata ---- ####
  modelstats_all <- data.table::fread(modelstats_file)
  
  fixed_col <- c("snp","gene","dist")
  stats_col <- setdiff(colnames(modelstats_all), fixed_col)
  
  pattern <- paste0("^", QTLtype, ".*", condition, "$")
  selected_stats_col <- grep(pattern, stats_col, value = TRUE)
  
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
  
  # Organize QTL table more neatly
  beta_col <- grep("beta", colnames(modelstats_subset), value = TRUE)
  se_col   <- grep("_se_", colnames(modelstats_subset), value = TRUE)
  pval_col <- grep("_p_",  colnames(modelstats_subset), value = TRUE)
  
  qtl_df <- modelstats_subset %>%
    dplyr::transmute(
      ID = snp,
      GENE = gene,
      DIST = dist,
      sdY_QTL = sdY_qtl,
      BETA_QTL = .data[[beta_col]],
      SE_QTL   = .data[[se_col]],
      P_QTL    = .data[[pval_col]]
    ) %>%
    dplyr::mutate(
      ID = stringr::str_to_lower(ID),
      Z_QTL = BETA_QTL / SE_QTL
    )
  
  # number of unique donors in the samples used for this QTL
  N_qtl <- meta %>% dplyr::pull(donor) %>% dplyr::n_distinct()
  
  #### ---- harmonize QTL and GWAS alleles ---- ####
  # 1) restrict to SNPs for this gene
  snp_list_QTL <- variant_tbl$ID %>% unique()
  snp_list_GWAS <- unique(gwas$SNP)
  snp_list <- intersect(snp_list_QTL, snp_list_GWAS)
  
  gwas_sub <- gwas %>%
    dplyr::filter(SNP %in% snp_list)
  
  # 2) Join GWAS to variant table (for REF/ALT/POS), only keep SNPs present in both GWAS and QTL
  merged <- gwas_sub %>%
    dplyr::inner_join(variant_tbl, by = c("SNP" = "ID")) %>%
    dplyr::inner_join(qtl_df, by = c("SNP" = "ID")) %>%
    dplyr::mutate(
      REF = stringr::str_to_upper(REF),
      ALT = stringr::str_to_upper(ALT),
      A1  = stringr::str_to_upper(A1),
      A2  = stringr::str_to_upper(A2)
    )
  
  # 3) Harmonize alleles
  # QTL effect allele = ALT
  # GWAS reports effect for A1 (assumed effect allele, consistent with ORX)
  is_ambiguous_pair <- function(ref, alt) {
    (ref == "A" & alt == "T") |
      (ref == "T" & alt == "A") |
      (ref == "C" & alt == "G") |
      (ref == "G" & alt == "C")
  }
  
  harmonized_all <- merged %>%
    dplyr::mutate(
      # GWAS beta on log-OR scale
      beta_gwas_raw    = log(ORX),
      varbeta_gwas     = SE^2,
      beta_qtl         = BETA_QTL,
      varbeta_qtl      = SE_QTL^2,
      ref = REF,
      alt = ALT,
      a1  = A1,
      a2  = A2,
      
      allele_status = dplyr::case_when(
        a1 == alt & a2 == ref ~ "aligned",
        a1 == ref & a2 == alt ~ "flip_gwas",
        TRUE                  ~ "mismatch"
      ),
      ambiguous = is_ambiguous_pair(ref, alt),
      
      beta_gwas = dplyr::if_else(allele_status == "flip_gwas", -beta_gwas_raw, beta_gwas_raw)
    )
  
  # 4) Filter to keep usable SNPs
  harmonized <- harmonized_all %>%
    dplyr::filter(allele_status %in% c("aligned", "flip_gwas")) %>%
    dplyr::filter(!ambiguous) %>%   # recommended
    dplyr::distinct(SNP, .keep_all = TRUE)
  
  message("[", QTLtype, " ", condition, "] Harmonized SNPs kept: ", nrow(harmonized))
  
  #### ---- make locus plot ---- ####
  this_trait <- paste0(QTLtype,":",condition)
  p <- plot_locus_track(harmonized, this_trait)
  #### ---- run coloc.abf ---- ####
  # ---- 4.1 Build coloc input tables (one row per SNP) ----
  coloc_tbl <- harmonized %>%
    dplyr::transmute(
      snp = SNP,
      
      # QTL dataset (quantitative)
      beta_qtl    = beta_qtl,
      varbeta_qtl = SE_QTL^2,
      maf_qtl     = MAF_QTL,
      
      # GWAS dataset (case-control; beta = log(OR) aligned to ALT already)
      beta_gwas    = beta_gwas,
      varbeta_gwas = varbeta_gwas,
      maf_gwas     = MAF
    ) %>%
    dplyr::filter(
      !is.na(beta_qtl), !is.na(varbeta_qtl), !is.na(maf_qtl),
      !is.na(beta_gwas), !is.na(varbeta_gwas), !is.na(maf_gwas)
    ) %>%
    dplyr::distinct(snp, .keep_all = TRUE)
  
  # Basic sanity checks
  if (nrow(coloc_tbl) < 50) {
    warning("[", QTLtype, " ", condition, "] Only ", nrow(coloc_tbl),
            " SNPs available for coloc. This may be underpowered/unstable.")
    return(NULL)
  }
  
  # coloc requires MAF in (0,1). Avoid exact 0 or 1
  eps <- 1e-6
  coloc_tbl <- coloc_tbl %>%
    dplyr::mutate(
      maf_qtl  = pmin(pmax(maf_qtl,  eps), 0.5 - eps),  # <-- safer upper bound
      maf_gwas = pmin(pmax(maf_gwas, eps), 0.5 - eps)
    )
  
  # ---- 4.2 Construct coloc datasets ----
  dataset_qtl <- list(
    snp     = coloc_tbl$snp,
    beta    = coloc_tbl$beta_qtl,
    varbeta = coloc_tbl$varbeta_qtl,
    maf     = coloc_tbl$maf_qtl,
    type    = "quant",
    N       = N_qtl,
    sdY     = sdY_qtl
  )
  
  dataset_gwas <- list(
    snp     = coloc_tbl$snp,
    beta    = coloc_tbl$beta_gwas,
    varbeta = coloc_tbl$varbeta_gwas,
    maf     = coloc_tbl$maf_gwas,
    type    = "cc",
    N       = N_gwas,
    s       = s_gwas
  )
  
  # ---- 4.3 Run coloc ----
  coloc_res <- coloc::coloc.abf(dataset1 = dataset_qtl, dataset2 = dataset_gwas)
  
  # ---- 4.4 Summaries ----
  summary_tbl <- tibble::tibble(
    ct   = ct,
    gene = this_gene,
    QTLtype = QTLtype,
    condition = condition,
    nsnps = coloc_res$summary[["nsnps"]],
    PP.H0 = coloc_res$summary[["PP.H0.abf"]],
    PP.H1 = coloc_res$summary[["PP.H1.abf"]],
    PP.H2 = coloc_res$summary[["PP.H2.abf"]],
    PP.H3 = coloc_res$summary[["PP.H3.abf"]],
    PP.H4 = coloc_res$summary[["PP.H4.abf"]]
  )
  
  # SNP-level table
  snp_res_tbl <- coloc_res$results %>%
    as_tibble()
  
  if ("snp" %in% names(snp_res_tbl)) {
    # already fine
  } else if ("SNP" %in% names(snp_res_tbl)) {
    snp_res_tbl <- snp_res_tbl %>% dplyr::rename(snp = SNP)
  } else {
    stop("coloc_res$results has no snp/SNP column. Columns are: ",
         paste(names(snp_res_tbl), collapse = ", "))
  }
  
  snp_res_tbl <- snp_res_tbl %>%
    dplyr::left_join(coloc_tbl, by = c("snp" = "snp")) %>%
    dplyr::mutate(ct = ct, gene = this_gene, QTLtype = QTLtype, condition = condition)
  
  # ---- 4.5 Output ----
  return(list(
    summary_tbl = summary_tbl,
    snp_res_tbl = snp_res_tbl,
    coloc_res   = coloc_res,
    plot        = p
  ))
}

#### ---- loop over QTL traits ---- ####
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

coloc_master_summary <- tibble::tibble()
coloc_master_snps <- tibble::tibble()
plot_list <- list()

for (i in seq_len(nrow(trait_grid))) {
  QTLtype_i   <- trait_grid$QTLtype[[i]]
  condition_i <- trait_grid$condition[[i]]
  
  res <- run_coloc_one_trait(
    ct = ct,
    this_gene = this_gene,
    QTLtype = QTLtype_i,
    condition = condition_i,
    chr = chr,
    dir = dir,
    pairs = pairs,
    genotype_all = genotype_all,
    meta_all = meta_all,
    modelstats_file = modelstats_file,
    vst_file = vst_file
  )
  
  if (!is.null(res)) {
    # write outputs per trait
    out_dir <- file.path(dir, "coloc", "vitiligo", condition_i, QTLtype_i)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    summary_file <- paste0("coloc_summary_", this_gene, ".tsv")
    snps_file    <- paste0("coloc_snps_", this_gene, ".tsv")
    
    readr::write_tsv(res$summary_tbl, file.path(out_dir, summary_file))
    readr::write_tsv(res$snp_res_tbl, file.path(out_dir, snps_file))
    
    message("Wrote: ", out_dir)
    
    coloc_master_summary <- dplyr::bind_rows(coloc_master_summary, res$summary_tbl)
    coloc_master_snps    <- dplyr::bind_rows(coloc_master_snps, res$snp_res_tbl)
    
    # make plots
    trait_key <- paste0(QTLtype_i, ":", condition_i)
    
    if (!is.null(res$plot)) {
      plot_list[[trait_key]] <- res$plot
    }
  }
}


#### ---- Plot ---- ####
pdf_dir <- file.path(dir, "coloc", "vitiligo", "plots")
pdf_file <- file.path(pdf_dir, paste0(ct, "_", this_gene, ".locus_tracks.pdf"))
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

if (length(plot_list) == 0) {
  message("No plots to write for ", ct, " ", this_gene)
} else {
  grDevices::pdf(pdf_file, width = 9.3, height = 8)
  for (nm in names(plot_list)) {
    print(plot_list[[nm]])
  }
  grDevices::dev.off()
  message("Wrote locus track PDF: ", pdf_file)
}

w <- warnings()
if (!is.null(w)) {
  message("---- warnings() ----")
  print(w)
}
