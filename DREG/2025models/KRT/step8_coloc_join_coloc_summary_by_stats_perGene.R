#!/usr/bin/env Rscript
# step8_coloc_join_coloc_summary_by_stats_perGene.R
# Join coloc outputs with QTL and GWAS stats for later discovery

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
  library(coloc)   # coloc 5.2.3
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(gridExtra)
})

# ------------------------------
# CLI args
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("
Usage:
  Rscript step8_coloc_join_coloc_summary_by_stats_perGene.R <ct> <gene> <GWAS_trait> 

Examples:
  Rscript step8_coloc_join_coloc_summary_by_stats_perGene.R MEL PIGX psoriasis 
", "\n")
  quit(save = "no", status = 1)
}

ct <- args[[1]]
this_gene <- args[[2]]
GWAS_trait <- args[[3]]

# # for testing only
# ct <- "KRT"
# this_gene <- "LCE3A"
# GWAS_trait <- "psoriasis"

# ------------------------------
# Global paths 
# ------------------------------
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)

chunk_id_lookup <- data.table::fread(paste0(dir, "/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[chunk_id_lookup$gene == this_gene, ]$chunk)
chunk_id <- sprintf("%03d", chunk_id)

pair_file <- paste0(dir, "/chunks/pairs_chunk_", chunk_id, ".tsv")
geno_file <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
modelstats_file <- paste0(dir, "/results/result_", chunk_id, ".tsv")
meta_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
rm(chunk_id_lookup)

# ------------------------------
# Load pairs / genotype / meta
# ------------------------------
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

# ------------------------------
# Fetch genome locus range
# ------------------------------
chr <- unique(pairs$gene_chr)
ciswindow_left <- unique(pairs$gene_start) - 500000
ciswindow_right <- unique(pairs$gene_start) + 500000
this_range <- paste(chr, ciswindow_left, ciswindow_right, sep="_")

# ------------------------------
# Standardized GWAS path 
# ------------------------------
gwas_stats_dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/", GWAS_trait)
gwas_stats_file <- paste0("standardized_", GWAS_trait, ".", chr, ".tsv.gz")
standardized_gwas_path <- file.path(gwas_stats_dir, gwas_stats_file)

# ------------------------------
# Small utils
# ------------------------------
stop2 <- function(...) stop(paste0(...), call. = FALSE)

is_finite01 <- function(x) is.finite(x) & x > 0 & x < 1

clamp_maf <- function(x, eps = 1e-6) {
  x <- as.numeric(x)
  pmin(pmax(x, eps), 0.5 - eps)
}

fix_p_zeros <- function(p) {
  p <- as.numeric(p)
  if (!any(is.finite(p) & p > 0, na.rm = TRUE)) return(p)
  min_nonzero <- min(p[p > 0], na.rm = TRUE)
  if (!is.finite(min_nonzero)) return(p)
  p[p == 0] <- min_nonzero / 2
  p
}

is_ambiguous_pair <- function(ref, alt) {
  ref <- toupper(ref); alt <- toupper(alt)
  (ref == "A" & alt == "T") |
    (ref == "T" & alt == "A") |
    (ref == "C" & alt == "G") |
    (ref == "G" & alt == "C")
}

# 1-row -> key/value table
transpose_for_table <- function(x) {
  stopifnot(nrow(x) == 1)
  data.frame(
    field = names(x),
    value = as.character(x[1, ]),
    stringsAsFactors = FALSE
  )
}

slice_rows <- function(df, from, to) {
  n <- nrow(df)
  if (from > n) return(df[0, , drop = FALSE])
  df[from:min(to, n), , drop = FALSE]
}

# Wrap any grob with a border (outline)
with_border <- function(g, lwd = 1) {
  grobTree(
    rectGrob(gp = gpar(fill = NA, col = "black", lwd = lwd)),
    g
  )
}

make_tbl_grob <- function(df, title) {
  if (nrow(df) == 0) {
    g <- arrangeGrob(
      textGrob(paste0(title, "\n(no rows)"),
               x = 0, just = "left",
               gp = gpar(fontface = "bold", fontsize = 10)),
      ncol = 1
    )
    return(with_border(g))
  }
  
  g <- arrangeGrob(
    textGrob(title, x = 0, just = "left",
             gp = gpar(fontface = "bold", fontsize = 10)),
    tableGrob(
      df, rows = NULL,
      theme = ttheme_minimal(
        core = list(fg_params = list(hjust = 0, x = 0.02, fontsize = 8)),
        colhead = list(fg_params = list(fontface = "bold", fontsize = 8))
      )
    ),
    ncol = 1,
    heights = unit.c(unit(0.25, "in"), unit(1, "null"))
  )
  
  with_border(g)
}

# New layout: col1=t1, col2=t2/t3, col3=t4/t5
make_summary_page <- function(out_summary) {
  kv <- transpose_for_table(out_summary)
  
  t1 <- make_tbl_grob(slice_rows(kv,  1, 26), "Table 1 (main)")
  t2 <- make_tbl_grob(slice_rows(kv, 27, 34), "Table 2 (best QTL)")
  t3 <- make_tbl_grob(slice_rows(kv, 35, 42), "Table 3 (best harmonized QTL)")
  t4 <- make_tbl_grob(slice_rows(kv, 43, 49), "Table 4 (best harmonized GWAS)")
  t5 <- make_tbl_grob(slice_rows(kv, 50, 57), "Table 5 (best GWAS)")
  
  # Layout matrix with 2 rows x 3 cols:
  # Row1: t1 | t2 | t4
  # Row2: t1 | t3 | t5
  layout <- rbind(
    c(1, 2, 3),
    c(1, 5, 4)
  )
  
  arrangeGrob(
    t1, t2, t3, t4, t5,
    layout_matrix = layout,
    widths  = c(0.40, 0.30, 0.30),
    heights = c(0.50, 0.50)
  )
}

# ------------------------------
# LD helper
# ------------------------------
fast_r2_to_lead <- function(geno_sub, lead_vec) {
  X <- as.matrix(geno_sub)
  y <- as.numeric(lead_vec)
  
  row_means <- rowMeans(X, na.rm = TRUE)
  Xc <- X - row_means
  
  yc <- y - mean(y, na.rm = TRUE)
  
  okY <- !is.na(yc)
  Xc[, !okY] <- NA
  
  okXY <- !is.na(Xc)
  n_pair <- rowSums(okXY)
  
  X0 <- Xc
  X0[!okXY] <- 0
  
  cov_num <- as.vector(X0 %*% yc)
  ssx <- rowSums(X0^2)
  
  y2 <- yc^2
  ssy <- as.vector(okXY %*% y2)
  
  denom <- sqrt(ssx * ssy)
  cor <- cov_num / denom
  cor[denom == 0 | n_pair < 2] <- NA_real_
  
  cor^2
}

# ------------------------------
# GWAS ingestion (STANDARDIZED, already per-chr file)
# ------------------------------
read_standardized_gwas <- function(path, range, verbose = TRUE) {
  
  # filter to SNPs within cis-window range
  this_chr <- strsplit(range, "_")[[1]][1]
  this_ciswindow_left <- as.numeric(strsplit(range, "_")[[1]][2])
  this_ciswindow_right <- as.numeric(strsplit(range, "_")[[1]][3])
  
  # read GWAS data table for specified range
  dt <- data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    data.table = TRUE,
    showProgress = FALSE,
    na.strings = c("NA", "", "NaN")
  ) %>%
    dplyr::filter(pos >= this_ciswindow_left & pos <= this_ciswindow_right)
  
  need <- c("trait","type","snp","rsid","variant_id","chr","pos","ea","oa",
            "beta","se","varbeta","p","MAF","N","num_cases","num_controls","s","sdY")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop2("Standardized GWAS missing required columns: ", paste(miss, collapse = ", "),
          "\nFound: ", paste(names(dt), collapse = ", "),
          "\nFile: ", path)
  }
  
  dt[, trait := as.character(trait)]
  dt[, type := as.character(type)]
  dt[, snp := stringr::str_to_lower(as.character(snp))]
  dt[, rsid := stringr::str_to_lower(as.character(rsid))]
  dt[, variant_id := stringr::str_to_lower(as.character(variant_id))]
  dt[, chr := as.integer(chr)]
  dt[, pos := as.integer(pos)]
  dt[, ea := stringr::str_to_upper(as.character(ea))]
  dt[, oa := stringr::str_to_upper(as.character(oa))]
  dt[, beta := as.numeric(beta)]
  dt[, se := as.numeric(se)]
  dt[, varbeta := as.numeric(varbeta)]
  dt[, p := as.numeric(p)]
  dt[, MAF := as.numeric(MAF)]
  dt[, N := as.numeric(N)]
  dt[, num_cases := as.numeric(num_cases)]
  dt[, num_controls := as.numeric(num_controls)]
  dt[, s := as.numeric(s)]
  dt[, sdY := as.numeric(sdY)]
  
  # clean essentials
  dt <- dt[!is.na(snp) & snp != ""]
  dt[(!is.finite(varbeta) | varbeta <= 0) & is.finite(se) & se > 0, varbeta := se^2]
  dt[, p := fix_p_zeros(p)]
  
  if (verbose) {
    message(sprintf(
      "[GWAS stats loaded for locus] %s | rows=%d | chr=%s | type=%s | snp missing=%d",
      basename(path), nrow(dt),
      paste(unique(dt$chr), collapse = ","),
      paste(unique(dt$type), collapse = ","),
      sum(is.na(dt$snp) | dt$snp == "")
    ))
  }
  dt[]
}

derive_gwas_globals <- function(gwas_dt, cases_total = NA_real_, controls_total = NA_real_,
                                N_override = NA_real_, sdY_override = NA_real_, verbose = TRUE) {
  gwas_type <- unique(gwas_dt$type)
  if (length(gwas_type) != 1) {
    message("[WARN] GWAS has multiple 'type' values: ", paste(gwas_type, collapse = ", "),
            " | Using first: ", gwas_type[1])
    gwas_type <- gwas_type[1]
  }
  if (!(gwas_type %in% c("cc","quant"))) stop2("GWAS standardized 'type' must be 'cc' or 'quant'. Got: ", gwas_type)
  
  # N
  N_ok <- gwas_dt$N[is.finite(gwas_dt$N) & gwas_dt$N > 0]
  if (length(N_ok) > 0) {
    N_global <- stats::median(N_ok)
    N_source <- "per_row_median"
    if (verbose) message(sprintf("[GWAS] N median(per-row) = %.0f", N_global))
  } else if (is.finite(N_override) && N_override > 0) {
    N_global <- N_override
    N_source <- "N_override"
    if (verbose) message(sprintf("[GWAS] N_override = %.0f", N_global))
  } else {
    stop2("[GWAS] No usable N found in standardized file and no N_override provided.")
  }
  
  s_global <- NA_real_
  s_source <- NA_character_
  sdY_global <- NA_real_
  sdY_source <- NA_character_
  
  if (gwas_type == "cc") {
    s_ok <- gwas_dt$s[is_finite01(gwas_dt$s)]
    if (length(s_ok) > 0) {
      s_global <- stats::median(s_ok)
      s_source <- "per_row_median"
      if (verbose) message(sprintf("[GWAS cc] s median(per-row) = %.5f", s_global))
    } else {
      ok <- is.finite(gwas_dt$num_cases) & gwas_dt$num_cases > 0 &
        is.finite(gwas_dt$num_controls) & gwas_dt$num_controls > 0
      if (any(ok)) {
        s_vec <- gwas_dt$num_cases[ok] / (gwas_dt$num_cases[ok] + gwas_dt$num_controls[ok])
        s_vec <- s_vec[is_finite01(s_vec)]
        if (length(s_vec) > 0) {
          s_global <- stats::median(s_vec)
          s_source <- "num_cases_controls_median"
          if (verbose) message(sprintf("[GWAS cc] s median(from num_cases/num_controls) = %.5f", s_global))
        }
      }
    }
    
    if (!is_finite01(s_global)) {
      if (is.finite(cases_total) && is.finite(controls_total) && cases_total > 0 && controls_total > 0) {
        s_global <- cases_total / (cases_total + controls_total)
        s_source <- "cases_total_controls_total"
        if (verbose) message(sprintf("[GWAS cc] s from provided totals = %.5f (cases=%.0f controls=%.0f)",
                                     s_global, cases_total, controls_total))
      } else {
        stop2("[GWAS cc] No usable s found; provide cases_total and controls_total.")
      }
    }
  } else {
    sdY_ok <- gwas_dt$sdY[is.finite(gwas_dt$sdY) & gwas_dt$sdY > 0]
    if (length(sdY_ok) > 0) {
      sdY_global <- stats::median(sdY_ok)
      sdY_source <- "per_row_median"
      if (verbose) message(sprintf("[GWAS quant] sdY median(per-row) = %.5f", sdY_global))
    } else if (is.finite(sdY_override) && sdY_override > 0) {
      sdY_global <- sdY_override
      sdY_source <- "sdY_override"
      if (verbose) message(sprintf("[GWAS quant] sdY_override = %.5f", sdY_global))
    } else {
      stop2("[GWAS quant] No usable sdY found; provide sdY_override.")
    }
  }
  
  list(
    type = gwas_type,
    N = N_global, N_source = N_source,
    s = s_global, s_source = s_source,
    sdY = sdY_global, sdY_source = sdY_source
  )
}

# ------------------------------
# Load GWAS ONCE (before loop)
# ------------------------------
message(sprintf("[GWAS] Loading once: %s", standardized_gwas_path))
gwas_dt <- read_standardized_gwas(standardized_gwas_path, this_range, verbose = TRUE) %>%
  arrange(p)
GWAS_best_rsid <- gwas_dt[1, ]$rsid
GWAS_best_pos <- gwas_dt[1, ]$pos
GWAS_best_ea <- gwas_dt[1, ]$ea
GWAS_best_oa <- gwas_dt[1, ]$oa
GWAS_best_beta <- gwas_dt[1, ]$beta
GWAS_best_se <- gwas_dt[1, ]$se
GWAS_best_p <- gwas_dt[1, ]$p
GWAS_best_maf <- gwas_dt[1, ]$MAF

if (GWAS_trait == "sunburn") {
  opt <- list(
    N_override = as.numeric(unique(gwas_dt$N)),
    cases_total = NA_real_,
    controls_total = NA_real_,
    sdY_override = 0.483
  )
}

if (GWAS_trait == "skin_pigmentation") {
  opt <- list(
    N_override = 415018,
    cases_total = NA_real_,
    controls_total = NA_real_,
    sdY_override = 0.556
  )
}

if (GWAS_trait == "height") {
  opt <- list(
    N_override = 424305,
    cases_total = NA_real_,
    controls_total = NA_real_,
    sdY_override = 1.289
  )
}

gwas_globals <- derive_gwas_globals(
  gwas_dt,
  cases_total = opt$cases_total,
  controls_total = opt$controls_total,
  N_override = opt$N_override,
  sdY_override = opt$sdY_override,
  verbose = TRUE
)

# crosswalk (for completeness; prefer_snp fixed to rsid, so mostly unused)
gwas_xwalk <- unique(gwas_dt[, .(snp, rsid, variant_id)])
gwas_xwalk[, rsid := as.character(rsid)]
gwas_xwalk[, variant_id := as.character(variant_id)]

# ------------------------------
# Core per-(QTLtype, condition) coloc (uses preloaded gwas_dt + globals)
# ------------------------------
merge_coloc_result_with_QTL_and_GWAS_stats <- function(
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
    gwas_dt,
    gwas_globals,
    nsnps_min = 50,
    verbose = TRUE
) {
  # ----------------------------
  # QTL meta selection
  # ----------------------------
  kept_condition <- if (QTLtype == "eQTL") condition else unique(c("PBS", condition))
  
  meta <- meta_all %>%
    dplyr::filter(celltype == ct, condition %in% kept_condition) %>%
    dplyr::select(sample, donor, condition)
  
  # ----------------------------
  # Genotype-derived REF/ALT + MAF_QTL
  # ----------------------------
  fixed_col <- c("CHROM","POS","ID","REF","ALT")
  donor_col <- setdiff(colnames(genotype_all), fixed_col)
  
  geno_mat <- genotype_all %>%
    dplyr::select(all_of(donor_col)) %>%
    dplyr::mutate(dplyr::across(everything(), as.numeric)) %>%
    as.matrix()
  rownames(geno_mat) <- stringr::str_to_lower(genotype_all$ID)
  
  alt_af <- rowMeans(geno_mat, na.rm = TRUE) / 2
  maf_qtl <- pmin(alt_af, 1 - alt_af)
  n_called <- rowSums(!is.na(geno_mat))
  
  maf_tbl <- tibble::tibble(
    rsid = names(alt_af) %>% stringr::str_to_lower(),
    AF_ALT_QTL = unname(alt_af),
    MAF_QTL = unname(maf_qtl),
    N_CALLED = unname(n_called)
  )
  
  variant_tbl <- genotype_all %>%
    dplyr::select(POS, ID, REF, ALT) %>%
    dplyr::mutate(
      rsid = stringr::str_to_lower(ID),
      REF = stringr::str_to_upper(REF),
      ALT = stringr::str_to_upper(ALT)
    ) %>%
    dplyr::left_join(maf_tbl, by = "rsid")
  
  # ----------------------------
  # QTL stats slice from model_stats
  # ----------------------------
  modelstats_all <- data.table::fread(modelstats_file, showProgress = FALSE)
  
  fixed_cols <- c("snp","gene","dist")
  stats_cols <- setdiff(colnames(modelstats_all), fixed_cols)
  pattern <- paste0("^", QTLtype, ".*", condition, "$")
  selected_stats_cols <- grep(pattern, stats_cols, value = TRUE)
  
  modelstats_subset <- modelstats_all[, c(fixed_cols, selected_stats_cols), with = FALSE] %>%
    dplyr::mutate(
      snp = stringr::str_to_lower(snp),
      key = paste0(gene, "_", snp)
    ) %>%
    dplyr::filter(key %in% pairs$key)
  
  # identify beta/se/p columns for this slice
  beta_col <- if (QTLtype == "eQTL") {
    grep("^eQTL_beta_", colnames(modelstats_subset), value = TRUE)
  } else {
    grep("^reQTL_dbeta_", colnames(modelstats_subset), value = TRUE)
  }
  se_col <- grep("_se_", colnames(modelstats_subset), value = TRUE)
  p_col  <- grep("_p_", colnames(modelstats_subset), value = TRUE)
  
  qtl_df <- modelstats_subset %>%
    dplyr::transmute(
      rsid = snp,
      gene = gene,
      dist = as.numeric(dist),
      beta_qtl = as.numeric(.data[[beta_col]]),
      se_qtl = as.numeric(.data[[se_col]]),
      varbeta_qtl = as.numeric(.data[[se_col]])^2,
      p_qtl = fix_p_zeros(as.numeric(.data[[p_col]]))
    ) %>%
    left_join( . , variant_tbl, by="rsid") %>%
    arrange(p_qtl)
  
  # ----------------------------
  # fetch best QTL result
  # ----------------------------
  QTL_best_rsid <- qtl_df[1, ]$rsid
  QTL_best_pos <- qtl_df[1, ]$POS
  QTL_best_REF <- qtl_df[1, ]$REF
  QTL_best_ALT <- qtl_df[1, ]$ALT
  QTL_best_beta <- qtl_df[1, ]$beta_qtl
  QTL_best_se <- qtl_df[1, ]$se_qtl
  QTL_best_p <- qtl_df[1, ]$p_qtl
  QTL_best_AF_ALT <- qtl_df[1, ]$AF_ALT_QTL

  # ----------------------------
  # GWAS subset
  # ----------------------------
  gwas_df <- as.data.frame(gwas_dt) %>%
    dplyr::transmute(
      rsid = rsid,                       # explicit
      snp_key = snp,                     # should match rsid per your standardization choice
      ea = ea,
      oa = oa,
      beta_gwas_raw = beta,
      varbeta_gwas = varbeta,
      p_gwas = p,
      MAF_gwas_raw = MAF
    )
  
  # ----------------------------
  # Join QTL (rsid) <-> GWAS (snp_key/rsid) <-> variant_tbl
  # Primary: qtl$rsid == gwas$snp_key; rescue: qtl$rsid == gwas$rsid
  # ----------------------------
  merged_primary <- qtl_df %>%
    dplyr::inner_join(gwas_df, by = c("rsid" = "snp_key")) 
  
  merged <- merged_primary
  
  message(sprintf("[JOIN primary rsid] matched=%d (qtl=%d gwas=%d)", nrow(merged_primary), nrow(qtl_df), nrow(gwas_df)))
  
  # ----------------------------
  # Harmonize alleles to ALT (QTL effect allele = ALT dosage)
  # ----------------------------
  harmonized_all <- merged %>%
    dplyr::mutate(
      ref = stringr::str_to_upper(REF),
      alt = stringr::str_to_upper(ALT),
      ea  = stringr::str_to_upper(ea),
      oa  = stringr::str_to_upper(oa),
      ambiguous = is_ambiguous_pair(ref, alt),
      allele_status = dplyr::case_when(
        ea == alt & oa == ref ~ "aligned",
        ea == ref & oa == alt ~ "flip_gwas",
        TRUE ~ "mismatch"
      ),
      beta_gwas = dplyr::if_else(allele_status == "flip_gwas", -beta_gwas_raw, beta_gwas_raw),
      maf_qtl = MAF_QTL,
      maf_gwas = dplyr::case_when(
        is.finite(MAF_gwas_raw) & MAF_gwas_raw > 0 & MAF_gwas_raw < 0.5 ~ MAF_gwas_raw,
        TRUE ~ NA_real_
      ),
      maf_source = dplyr::case_when(
        is.finite(MAF_gwas_raw) & MAF_gwas_raw > 0 & MAF_gwas_raw < 0.5 ~ "gwas_MAF",
        TRUE ~ "missing"
      )
    ) %>%
    dplyr::mutate(
      maf_gwas = dplyr::if_else(is.na(maf_gwas) & is.finite(maf_qtl), maf_qtl, maf_gwas),
      maf_source = dplyr::if_else(maf_source == "missing" & !is.na(maf_gwas), "fallback_qtl_MAF", maf_source)
    )
  
  harmonized <- harmonized_all %>%
    dplyr::filter(allele_status %in% c("aligned", "flip_gwas")) %>%
    dplyr::filter(!ambiguous) %>%
    dplyr::distinct(rsid, .keep_all = TRUE)
  

  message(sprintf(
    "[HARM] %s %s:%s | kept=%d | aligned=%d flip=%d | mismatch=%d | ambig_dropped=%d | maf_fallback=%d",
    QTLtype, condition, this_gene,
    nrow(harmonized),
    sum(harmonized$allele_status == "aligned"),
    sum(harmonized$allele_status == "flip_gwas"),
    sum(harmonized_all$allele_status == "mismatch", na.rm = TRUE),
    sum(harmonized_all$ambiguous, na.rm = TRUE),
    sum(harmonized$maf_source == "fallback_qtl_MAF", na.rm = TRUE)
  ))
  
  
  coloc_tbl <- harmonized %>%
    dplyr::filter(
      is.finite(beta_qtl),
      is.finite(varbeta_qtl) & varbeta_qtl > 0,
      is.finite(maf_qtl) & maf_qtl > 0 & maf_qtl < 0.5,
      is.finite(beta_gwas),
      is.finite(varbeta_gwas) & varbeta_gwas > 0,
      is.finite(maf_gwas) & maf_gwas > 0 & maf_gwas < 0.5,
      !is.na(rsid) & rsid != ""
    ) %>%
    dplyr::distinct(rsid, .keep_all = TRUE)
  
  if (nrow(coloc_tbl) < nsnps_min) {
    message(sprintf("[SKIP] %s %s %s:%s nsnps=%d < %d",
                    ct, this_gene, QTLtype, condition, nrow(coloc_tbl), nsnps_min))
    return(NULL)
  }
  
  coloc_tbl <- coloc_tbl %>%
    dplyr::mutate(
      maf_qtl = clamp_maf(maf_qtl),
      maf_gwas = clamp_maf(maf_gwas)
    )
  
  # ----------------------------
  # fetch best QTL result in harmonized (GWAS x QTL) dataset
  # ----------------------------
  coloc_tbl <- coloc_tbl %>% arrange(p_qtl)
  QTL_bestHarmonized_rsid <- coloc_tbl[1, ]$rsid
  QTL_bestHarmonized_pos <- coloc_tbl[1, ]$POS
  QTL_bestHarmonized_REF <- coloc_tbl[1, ]$REF
  QTL_bestHarmonized_ALT <- coloc_tbl[1, ]$ALT
  QTL_bestHarmonized_beta <- coloc_tbl[1, ]$beta_qtl
  QTL_bestHarmonized_se <- coloc_tbl[1, ]$se_qtl
  QTL_bestHarmonized_p <- coloc_tbl[1, ]$p_qtl
  QTL_bestHarmonized_AF_ALT <- coloc_tbl[1, ]$AF_ALT_QTL
  
  # ----------------------------
  # fetch best GWAS result in harmonized (GWAS x QTL) dataset
  # ----------------------------
  coloc_tbl <- coloc_tbl %>% arrange(p_gwas)
  GWAS_bestHarmonized_rsid <- coloc_tbl[1, ]$rsid
  GWAS_bestHarmonized_pos <- coloc_tbl[1, ]$POS
  GWAS_bestHarmonized_ea <- coloc_tbl[1, ]$ea
  GWAS_bestHarmonized_oa <- coloc_tbl[1, ]$oa
  GWAS_bestHarmonized_beta <- coloc_tbl[1, ]$beta_gwas_raw
  GWAS_bestHarmonized_p <- coloc_tbl[1, ]$p_gwas
  GWAS_bestHarmonized_maf <- coloc_tbl[1, ]$maf_gwas

  
  
  # return result
  # ----------------------------
  # Build a single-row summary output with fixed columns
  # (easy to rbind across runs and fwrite with stable header)
  # ----------------------------
  stats_summary <- data.table::data.table(
    gene = this_gene,
    
    nsnps_qtl = nrow(qtl_df),
    nsnps_gwas = nrow(gwas_df),
    nsnps_harmonized = nrow(harmonized),
    
    # Best QTL in raw QTL table
    QTL_best_rsid = QTL_best_rsid,
    QTL_best_pos = QTL_best_pos,
    QTL_best_REF = QTL_best_REF,
    QTL_best_ALT = QTL_best_ALT,
    QTL_best_beta = QTL_best_beta,
    QTL_best_se = QTL_best_se,
    QTL_best_p = QTL_best_p,
    QTL_best_AF_ALT = QTL_best_AF_ALT,
    
    # Best QTL in harmonized/coloc table
    QTL_bestHarmonized_rsid = QTL_bestHarmonized_rsid,
    QTL_bestHarmonized_pos = QTL_bestHarmonized_pos,
    QTL_bestHarmonized_REF = QTL_bestHarmonized_REF,
    QTL_bestHarmonized_ALT = QTL_bestHarmonized_ALT,
    QTL_bestHarmonized_beta = QTL_bestHarmonized_beta,
    QTL_bestHarmonized_se = QTL_bestHarmonized_se,
    QTL_bestHarmonized_p = QTL_bestHarmonized_p,
    QTL_bestHarmonized_AF_ALT = QTL_bestHarmonized_AF_ALT,
    
    # Best GWAS in harmonized/coloc table
    GWAS_bestHarmonized_rsid = GWAS_bestHarmonized_rsid,
    GWAS_bestHarmonized_pos = GWAS_bestHarmonized_pos,
    GWAS_bestHarmonized_ea = GWAS_bestHarmonized_ea,
    GWAS_bestHarmonized_oa = GWAS_bestHarmonized_oa,
    GWAS_bestHarmonized_beta = GWAS_bestHarmonized_beta,
    GWAS_bestHarmonized_p = GWAS_bestHarmonized_p,
    GWAS_bestHarmonized_maf = GWAS_bestHarmonized_maf
  )
  # ----------------------------
  # read coloc summary file
  # ----------------------------
  coloc_file=paste0(dir,"/coloc/summary/coloc_",GWAS_trait,"_",condition,"_",QTLtype,".txt")
  coloc_summary <- fread(coloc_file) %>% 
    dplyr::filter(gene==this_gene) %>%
    dplyr::mutate(gene_chr=unique(pairs$gene_chr)) %>%
    dplyr::mutate(gene_TSS=unique(pairs$gene_start))
  
  # get out_summary file ----
  out_summary <- left_join(coloc_summary, stats_summary, by="gene")
  
  return(out_summary)
} 


# ------------------------------
# Trait grid
# ------------------------------
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

# ------------------------------
# Main loop
# ------------------------------
coloc_master_summary <- tibble::tibble()
plot_list <- list()

message(sprintf(
  "[START] ct=%s gene=%s GWAS_trait=%s chr=%s gwas_file=%s",
  ct, this_gene, GWAS_trait, as.character(chr), standardized_gwas_path
))

plot_dir <- paste0(dir,"/coloc/",GWAS_trait,"/plots_table/")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

pdf(paste0(plot_dir,"/",ct,"_",this_gene,".table.pdf"), width = 11, height = 8.5)  # landscape letter
for (i in seq_len(nrow(trait_grid))) {
  QTLtype_i <- trait_grid$QTLtype[[i]]
  condition_i <- trait_grid$condition[[i]]
  
  res <- merge_coloc_result_with_QTL_and_GWAS_stats(
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
    gwas_dt = gwas_dt,
    gwas_globals = gwas_globals,
    nsnps_min = 50,
    verbose = TRUE
  )
  
  GWAS_best_summary <- data.table::data.table(
    gene = this_gene,
    GWAS_best_rsid = GWAS_best_rsid,
    GWAS_best_pos = GWAS_best_pos,
    GWAS_best_ea = GWAS_best_ea,
    GWAS_best_oa = GWAS_best_oa,
    GWAS_best_beta = GWAS_best_beta,
    GWAS_best_se = GWAS_best_se,
    GWAS_best_p = GWAS_best_p,
    GWAS_best_maf = GWAS_best_maf)
  
  res <- left_join(res, GWAS_best_summary, by="gene")
  
  if (!is.null(res)) {
    out_dir <- file.path(dir, "coloc", GWAS_trait, condition_i, QTLtype_i)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    summary_file <- paste0("joined_coloc_summary_", this_gene, ".tsv")
    
    readr::write_tsv(res, file.path(out_dir, summary_file))
    
    message("Wrote: ", out_dir)
    
    coloc_master_summary <- dplyr::bind_rows(coloc_master_summary, res)
    
    trait_key <- paste0(QTLtype_i, ":", condition_i)
    if (!is.null(res$plot)) plot_list[[trait_key]] <- res$plot
  }
  
  # ------------------------
  # make table into plot
  # ------------------------
  # New page + draw table
  grid.newpage()
  grid.draw(make_summary_page(res))
}
dev.off()

