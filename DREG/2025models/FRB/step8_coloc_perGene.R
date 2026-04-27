#!/usr/bin/env Rscript
# step8_coloc_perGene.R
# Run coloc per gene between QTL summary stats (model_stats) and *standardized* GWAS SSF.
# Outputs:
#  1) per-trait coloc summary table (PP.H0..PP.H4, nsnps) with ct/gene/QTLtype/condition
#  2) per-SNP coloc results table joined with key fields
#  3) locus plots (QTL p, GWAS p, coloc SNP.PP.H4) with LD coloring from genotype matrix
#  4) writes outputs into the same directory structure as old script
#
# Notes:
# - GWAS input must already be standardized to fixed schema (tab-separated) with literal "NA" for blanks.
# - "snp" column is the intended harmonization key; it depends on rsid
# - We align alleles to QTL ALT dosage (REF/ALT from genotype table), flipping GWAS beta if needed.
# - Strand-ambiguous SNPs (A/T, C/G) are dropped by default.

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
})

# ------------------------------
# CLI args
# ------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("
Usage:
  Rscript step8_coloc_perGene.R <ct> <gene> <GWAS_trait> 

Optional (in any order as key=value):
  cases_total=<num>          # cc only fallback if GWAS has no usable s/num_cases/num_controls
  controls_total=<num>
  N_override=<num>           # fallback if GWAS has no usable N
  sdY_override=<num>         # quant only fallback if GWAS has no usable sdY

Examples:
  Rscript step8_coloc_perGene.R MEL PIGX psoriasis 
  Rscript step8_coloc_perGene.R FRB ERAP2 atopic_dermatitis cases_total=42963 controls_total=408472 N_override=451435
", "\n")
  quit(save = "no", status = 1)
}

ct <- args[[1]]
this_gene <- args[[2]]
GWAS_trait <- args[[3]]

## for testing only
#ct <- "FRB"
#this_gene <- "RCC2-AS1"
#GWAS_trait <- "atopic_dermatitis"

# ------------------------------
# Global paths 
# ------------------------------
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)

chunk_id_lookup <- data.table::fread(paste0(dir, "/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[chunk_id_lookup$gene == this_gene, ]$chunk)
chunk_id <- sprintf("%03d", chunk_id)

pair_file <- paste0(dir, "/chunks/pairs_chunk_", chunk_id, ".tsv")
geno_file <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
vst_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
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
# Locus plotting (expects columns: snp, gene, dist, P_QTL, P_GWAS)
# ------------------------------
plot_locus_track <- function(
    harmonized,
    this_trait,
    geno_mat = NULL,
    coloc_snp_tbl = NULL,
    coloc_summary = NULL,
    LDanchor = c("lead_qtl", "lead_gwas", "lead_coloc_pph4"),
    ld_colname = "r2_to_anchor"
) {
  LDanchor <- match.arg(LDanchor)
  
  pick_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  get_pp_from_summary <- function(coloc_summary, key) {
    if (is.null(coloc_summary)) return(NA_real_)
    if (is.list(coloc_summary)) {
      if (!is.null(coloc_summary[[key]])) return(as.numeric(coloc_summary[[key]]))
      if (!is.null(names(coloc_summary)) && key %in% names(coloc_summary)) return(as.numeric(coloc_summary[[key]]))
    }
    if (is.numeric(coloc_summary) && !is.null(names(coloc_summary)) && key %in% names(coloc_summary)) {
      return(as.numeric(coloc_summary[[key]]))
    }
    NA_real_
  }
  
  snp_col    <- pick_col(harmonized, c("snp"))
  gwas_p_col <- pick_col(harmonized, c("P_GWAS"))
  qtl_p_col  <- pick_col(harmonized, c("P_QTL"))
  gene_col   <- pick_col(harmonized, c("gene"))
  dist_col   <- pick_col(harmonized, c("dist"))
  
  needed <- c(snp_col, gwas_p_col, qtl_p_col, gene_col, dist_col)
  if (any(is.na(needed))) {
    stop2("[plot_locus_track] Missing required columns. Available: ", paste(names(harmonized), collapse = ", "))
  }
  
  wide <- harmonized %>%
    dplyr::transmute(
      snp   = stringr::str_to_lower(.data[[snp_col]]),
      GENE  = .data[[gene_col]],
      DIST  = as.numeric(.data[[dist_col]]),
      P_GWAS = as.numeric(.data[[gwas_p_col]]),
      P_QTL  = as.numeric(.data[[qtl_p_col]])
    ) %>%
    dplyr::mutate(
      P_GWAS = fix_p_zeros(P_GWAS),
      P_QTL  = fix_p_zeros(P_QTL),
      log10p.gwas = -log10(P_GWAS),
      log10p.qtl  = -log10(P_QTL)
    )
  
  if (nrow(wide) == 0) return(NULL)
  
  lead_qtl <- wide$snp[which.min(wide$P_QTL)]
  if (!is.finite(wide$P_QTL[match(lead_qtl, wide$snp)])) lead_qtl <- wide$snp[1]
  
  lead_gwas <- wide$snp[which.min(wide$P_GWAS)]
  if (!is.finite(wide$P_GWAS[match(lead_gwas, wide$snp)])) lead_gwas <- wide$snp[1]
  
  coloc_snp_tbl <- tibble::as_tibble(coloc_snp_tbl)
  if (nrow(coloc_snp_tbl) > 0) {
    if (!("snp" %in% names(coloc_snp_tbl)) && ("SNP" %in% names(coloc_snp_tbl))) {
      coloc_snp_tbl <- coloc_snp_tbl %>% dplyr::rename(snp = SNP)
    }
    coloc_snp_tbl <- coloc_snp_tbl %>% dplyr::mutate(snp = stringr::str_to_lower(.data$snp))
  }
  
  pph4_col <- pick_col(coloc_snp_tbl, c("SNP.PP.H4", "PP.H4", "PPH4"))
  pph3_col <- pick_col(coloc_snp_tbl, c("SNP.PP.H3", "PP.H3", "PPH3"))
  
  coloc_map <- if (nrow(coloc_snp_tbl) > 0) {
    coloc_snp_tbl %>%
      dplyr::transmute(
        snp = .data$snp,
        coloc_pph4 = if (!is.na(pph4_col)) as.numeric(.data[[pph4_col]]) else NA_real_,
        coloc_pph3 = if (!is.na(pph3_col)) as.numeric(.data[[pph3_col]]) else NA_real_
      )
  } else {
    tibble::tibble(snp = character(), coloc_pph4 = numeric(), coloc_pph3 = numeric())
  }
  
  wide <- wide %>% dplyr::left_join(coloc_map, by = "snp")
  
  lead_coloc_pph4 <- if (!all(is.na(wide$coloc_pph4))) wide$snp[which.max(wide$coloc_pph4)] else NA_character_
  lead_coloc_pph3 <- if (!all(is.na(wide$coloc_pph3))) wide$snp[which.max(wide$coloc_pph3)] else NA_character_
  
  PPH4_global <- get_pp_from_summary(coloc_summary, "PP.H4.abf")
  if (is.na(PPH4_global)) PPH4_global <- get_pp_from_summary(coloc_summary, "PP.H4")
  PPH3_global <- get_pp_from_summary(coloc_summary, "PP.H3.abf")
  if (is.na(PPH3_global)) PPH3_global <- get_pp_from_summary(coloc_summary, "PP.H3")
  
  ld_anchor_snp <- dplyr::case_when(
    LDanchor == "lead_qtl"        ~ lead_qtl,
    LDanchor == "lead_gwas"       ~ lead_gwas,
    LDanchor == "lead_coloc_pph4" ~ lead_coloc_pph4,
    TRUE                          ~ lead_qtl
  )
  ld_anchor_label <- dplyr::case_when(
    LDanchor == "lead_qtl"        ~ "lead QTL",
    LDanchor == "lead_gwas"       ~ "lead GWAS",
    LDanchor == "lead_coloc_pph4" ~ "lead coloc PPH4",
    TRUE                          ~ "lead QTL"
  )
  
  wide[[ld_colname]] <- NA_real_
  has_ld <- FALSE
  if (!is.null(geno_mat) && !is.null(rownames(geno_mat))) {
    rownames(geno_mat) <- stringr::str_to_lower(rownames(geno_mat))
    snps_here <- intersect(wide$snp, rownames(geno_mat))
    if (length(snps_here) >= 3 && !is.na(ld_anchor_snp) && ld_anchor_snp %in% rownames(geno_mat)) {
      geno_sub <- geno_mat[snps_here, , drop = FALSE]
      lead_vec <- geno_mat[ld_anchor_snp, ]
      r2 <- fast_r2_to_lead(geno_sub, lead_vec)
      wide[[ld_colname]][match(snps_here, wide$snp)] <- r2
      has_ld <- TRUE
      wide <- wide %>% dplyr::arrange(.data[[ld_colname]])
    }
  }
  
  title_txt <- paste0("QTL/GWAS/coloc locus: ", unique(wide$GENE), " | ", this_trait,
                      " | LD anchor SNP = ", ld_anchor_label)
  
  subtitle_txt <- paste0(
    "lead_qtl=", lead_qtl,
    " | lead_gwas=", lead_gwas,
    " | lead_coloc_PPH4=", lead_coloc_pph4,
    " | PPH4=", ifelse(is.na(PPH4_global), "NA", formatC(PPH4_global, digits = 4, format = "f")),
    " | lead_coloc_PPH3=", lead_coloc_pph3,
    " | PPH3=", ifelse(is.na(PPH3_global), "NA", formatC(PPH3_global, digits = 4, format = "f"))
  )
  
  color_aes <- if (has_ld) ggplot2::aes(color = .data[[ld_colname]]) else ggplot2::aes()
  ld_scale <- ggplot2::scale_color_viridis_c(
    option = "viridis",
    limits = c(0, 1),
    oob = scales::squish,
    na.value = "grey80",
    name = "r² to anchor"
  )
  
  p_qtl <- ggplot2::ggplot(wide, ggplot2::aes(x = DIST, y = log10p.qtl)) +
    ggplot2::geom_point(color_aes, alpha = 0.6, size = 1.5) +
    ggplot2::labs(title = title_txt, subtitle = subtitle_txt,
                  x = "SNP distance to gene TSS (bp)", y = "-log10(p.qtl)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(angle = 45, hjust = 1))
  if (has_ld) p_qtl <- p_qtl + ld_scale
  
  p_gwas <- ggplot2::ggplot(wide, ggplot2::aes(x = DIST, y = log10p.gwas)) +
    ggplot2::geom_point(color_aes, alpha = 0.6, size = 1.5) +
    ggplot2::labs(title = paste0("GWAS locus: ", unique(wide$GENE)),
                  x = "SNP distance to gene TSS (bp)", y = "-log10(p.gwas)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  if (has_ld) p_gwas <- p_gwas + ld_scale
  
  p_coloc <- NULL
  if (!all(is.na(wide$coloc_pph4))) {
    p_coloc <- ggplot2::ggplot(wide, ggplot2::aes(x = DIST, y = coloc_pph4)) +
      ggplot2::geom_point(color_aes, alpha = 0.6, size = 1.5) +
      ggplot2::labs(title = "coloc SNP posterior (SNP.PP.H4)",
                    x = "SNP distance to gene TSS (bp)", y = "Posterior") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "none")
    if (has_ld) p_coloc <- p_coloc + ld_scale
  }
  
  if (!is.null(p_coloc)) {
    p_qtl / p_gwas / p_coloc + patchwork::plot_layout(heights = c(1, 1, 1))
  } else {
    p_qtl / p_gwas + patchwork::plot_layout(heights = c(1, 1))
  }
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
gwas_dt <- read_standardized_gwas(standardized_gwas_path, this_range, verbose = TRUE)

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
run_coloc_one_trait_preloadedGWAS <- function(
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
  
  # compute sdY_qtl
  y <- data.table::fread(vst_file, showProgress = FALSE) %>%
    dplyr::filter(gene == this_gene) %>%
    dplyr::select(all_of(meta$sample)) %>%
    t() %>% as.data.frame() %>% dplyr::pull(V1)
  
  sdY_qtl <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sdY_qtl) || sdY_qtl <= 0) {
    message(sprintf("[SKIP] %s %s %s:%s constant VST (sdY=%s)", ct, this_gene, QTLtype, condition, as.character(sdY_qtl)))
    return(NULL)
  }
  
  N_qtl <- meta %>% dplyr::pull(donor) %>% dplyr::n_distinct()
  if (!is.finite(N_qtl) || N_qtl < 10) {
    message(sprintf("[SKIP] %s %s %s:%s N_qtl too small (%s)", ct, this_gene, QTLtype, condition, as.character(N_qtl)))
    return(NULL)
  }
  
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
      varbeta_qtl = as.numeric(.data[[se_col]])^2,
      p_qtl = fix_p_zeros(as.numeric(.data[[p_col]]))
    )
  
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
    dplyr::inner_join(gwas_df, by = c("rsid" = "snp_key")) %>%
    dplyr::inner_join(variant_tbl, by = c("rsid" = "rsid"))
  
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
    dplyr::transmute(
      snp = rsid,
      beta_qtl = beta_qtl,
      varbeta_qtl = varbeta_qtl,
      maf_qtl = maf_qtl,
      beta_gwas = beta_gwas,
      varbeta_gwas = varbeta_gwas,
      maf_gwas = maf_gwas,
      gene = gene,
      dist = dist,
      P_QTL = p_qtl,
      P_GWAS = p_gwas,
      allele_status = allele_status,
      maf_source = maf_source
    ) %>%
    dplyr::filter(
      is.finite(beta_qtl),
      is.finite(varbeta_qtl) & varbeta_qtl > 0,
      is.finite(maf_qtl) & maf_qtl > 0 & maf_qtl < 0.5,
      is.finite(beta_gwas),
      is.finite(varbeta_gwas) & varbeta_gwas > 0,
      is.finite(maf_gwas) & maf_gwas > 0 & maf_gwas < 0.5,
      !is.na(snp) & snp != ""
    ) %>%
    dplyr::distinct(snp, .keep_all = TRUE)
  
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
  
  dataset_qtl <- list(
    snp = coloc_tbl$snp,
    beta = coloc_tbl$beta_qtl,
    varbeta = coloc_tbl$varbeta_qtl,
    maf = coloc_tbl$maf_qtl,
    type = "quant",
    N = N_qtl,
    sdY = sdY_qtl
  )
  
  dataset_gwas <- list(
    snp = coloc_tbl$snp,
    beta = coloc_tbl$beta_gwas,
    varbeta = coloc_tbl$varbeta_gwas,
    maf = coloc_tbl$maf_gwas,
    type = gwas_globals$type,
    N = gwas_globals$N
  )
  if (gwas_globals$type == "cc") {
    dataset_gwas$s <- gwas_globals$s
  } else {
    dataset_gwas$sdY <- gwas_globals$sdY
  }
  
  coloc_res <- tryCatch(
    coloc::coloc.abf(dataset1 = dataset_qtl, dataset2 = dataset_gwas),
    error = function(e) {
      message(sprintf("[ERROR coloc] %s %s %s:%s | %s", ct, this_gene, QTLtype, condition, e$message))
      NULL
    }
  )
  if (is.null(coloc_res)) return(NULL)
  
  summary_tbl <- tibble::tibble(
    ct = ct,
    gene = this_gene,
    QTLtype = QTLtype,
    condition = condition,
    GWAS_trait = GWAS_trait,
    prefer_snp = "rsid",
    
    nsnps = coloc_res$summary[["nsnps"]],
    PP.H0 = coloc_res$summary[["PP.H0.abf"]],
    PP.H1 = coloc_res$summary[["PP.H1.abf"]],
    PP.H2 = coloc_res$summary[["PP.H2.abf"]],
    PP.H3 = coloc_res$summary[["PP.H3.abf"]],
    PP.H4 = coloc_res$summary[["PP.H4.abf"]],
    
    N_qtl = N_qtl,
    sdY_qtl = sdY_qtl,
    gwas_type = gwas_globals$type,
    N_gwas = gwas_globals$N,
    s_gwas = gwas_globals$s,
    sdY_gwas = gwas_globals$sdY,
    N_gwas_source = gwas_globals$N_source,
    s_gwas_source = gwas_globals$s_source,
    sdY_gwas_source = gwas_globals$sdY_source
  )
  
  snp_res_tbl <- tibble::as_tibble(coloc_res$results)
  if (!("snp" %in% names(snp_res_tbl)) && ("SNP" %in% names(snp_res_tbl))) {
    snp_res_tbl <- snp_res_tbl %>% dplyr::rename(snp = SNP)
  }
  if (!("snp" %in% names(snp_res_tbl))) {
    stop2("coloc_res$results has no snp/SNP column. Columns: ", paste(names(snp_res_tbl), collapse = ", "))
  }
  
  snp_res_tbl <- snp_res_tbl %>%
    dplyr::left_join(
      coloc_tbl %>%
        dplyr::select(
          snp, gene, dist, P_QTL, P_GWAS,
          beta_qtl, varbeta_qtl, maf_qtl,
          beta_gwas, varbeta_gwas, maf_gwas,
          allele_status, maf_source
        ),
      by = "snp"
    ) %>%
    dplyr::mutate(
      ct = ct,
      gene = this_gene,
      QTLtype = QTLtype,
      condition = condition,
      GWAS_trait = GWAS_trait,
      prefer_snp = "rsid"
    )
  
  harmonized_for_plot <- coloc_tbl %>%
    dplyr::transmute(
      snp = snp,
      gene = gene,
      dist = dist,
      P_QTL = P_QTL,
      P_GWAS = P_GWAS
    )
  
  this_trait <- paste0(GWAS_trait, " | ", QTLtype, ":", condition)
  
  p <- plot_locus_track(
    harmonized = harmonized_for_plot,
    this_trait = this_trait,
    geno_mat = geno_mat,
    coloc_snp_tbl = coloc_res$results,
    coloc_summary = coloc_res$summary,
    LDanchor = "lead_qtl"
  )
  
  list(
    summary_tbl = summary_tbl,
    snp_res_tbl = snp_res_tbl,
    coloc_res = coloc_res,
    plot = p
  )
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
coloc_master_snps <- tibble::tibble()
plot_list <- list()

message(sprintf(
  "[START] ct=%s gene=%s GWAS_trait=%s chr=%s gwas_file=%s",
  ct, this_gene, GWAS_trait, as.character(chr), standardized_gwas_path
))

for (i in seq_len(nrow(trait_grid))) {
  QTLtype_i <- trait_grid$QTLtype[[i]]
  condition_i <- trait_grid$condition[[i]]
  
  res <- run_coloc_one_trait_preloadedGWAS(
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
    vst_file = vst_file,
    gwas_dt = gwas_dt,
    gwas_globals = gwas_globals,
    nsnps_min = 50,
    verbose = TRUE
  )
  
  if (!is.null(res)) {
    out_dir <- file.path(dir, "coloc", GWAS_trait, condition_i, QTLtype_i)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    summary_file <- paste0("coloc_summary_", this_gene, ".tsv")
    snps_file <- paste0("coloc_snps_", this_gene, ".tsv")
    
    readr::write_tsv(res$summary_tbl, file.path(out_dir, summary_file))
    readr::write_tsv(res$snp_res_tbl, file.path(out_dir, snps_file))
    message("Wrote: ", out_dir)
    
    coloc_master_summary <- dplyr::bind_rows(coloc_master_summary, res$summary_tbl)
    coloc_master_snps <- dplyr::bind_rows(coloc_master_snps, res$snp_res_tbl)
    
    trait_key <- paste0(QTLtype_i, ":", condition_i)
    if (!is.null(res$plot)) plot_list[[trait_key]] <- res$plot
  }
}

# ------------------------------
# Write plots PDF (same path pattern)
# ------------------------------
pdf_dir <- file.path(dir, "coloc", GWAS_trait, "plots")
pdf_file <- file.path(pdf_dir, paste0(ct, "_", this_gene, ".locus_tracks.pdf"))
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

if (length(plot_list) == 0) {
  message("No plots to write for ", ct, " ", this_gene)
} else {
  grDevices::pdf(pdf_file, width = 14, height = 12)
  for (nm in names(plot_list)) print(plot_list[[nm]])
  grDevices::dev.off()
  message("Wrote locus track PDF: ", pdf_file)
}

w <- warnings()
if (!is.null(w)) {
  message("---- warnings() ----")
  print(w)
}

message("[DONE]")
