#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(magrittr)
  library(readr)
  library(patchwork)
  library(tibble)
})

# ---------------------- Argument Parsing ----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("
Usage:
  Rscript step6_make_response_plot_for_any_pair.R <gene> <snp>

Example:
  Rscript step6_make_response_plot_for_any_pair.R ITGA1 rs2548496 \\
", "\n")
  quit(save = "no", status = 1)
}

basedir <- args[[1]] # "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots"
g <- args[[2]] # ITGA1
snp <- args[[3]] # rs2548496

# # toy example
# basedir <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots"
# g <- "IL1B"
# snp <- "rs4790797"


# ---------------------- Set-up --------------------- ####
letter <- substr(g, 1, 1)
dir <- paste0(basedir, "/",letter, "/", g, "/", snp)
GENO_FILE <- paste0(dir,"/genotype.txt")
STATS_FILE <- NULL
OUT_PDF <- paste0(basedir,"/temp_output/plot_response_",g,"_",snp,".pdf")
CPM_FILE   <- paste0(dir, "/cpm.txt")
META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"

message(paste0("Starting to plot for ", g, " : ", snp))
# ---------------------- Load Expression ---------------------- ####
message("[1] Loading CPM expression...")
cpm <- fread(CPM_FILE) %>%
  select(-c("gene","name")) %>%
  distinct(final_gene, .keep_all = TRUE) %>%
  column_to_rownames("final_gene") %>%
  rownames_to_column("gene_id") 

expr_long <- as.data.table(cpm) %>%
  melt(id.vars = "gene_id", variable.name = "sample", value.name = "CPM")
rm(cpm); invisible(gc())

# ---------------------- Metadata ---------------------- ####
message("[2] Loading metadata...")
meta <- fread(META_FILE)
expr_long <- expr_long %>%
  inner_join(meta, by = "sample") %>%
  mutate(condition = as.character(condition),
         celltype  = as.character(celltype))

# ---------------------- Genotype Dosage ---------------------- ####
message("[3] Loading genotype dosage...")
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")

geno_samples <- setdiff(colnames(geno_dt), fixed_cols)
expr_donors <- unique(expr_long$donor)
keep <- intersect(geno_samples, expr_donors)
geno_dt <- geno_dt[, c(fixed_cols, keep), with=FALSE]
setkey(geno_dt, ID)

dosage_to_label <- function(d, ref, alt){
  factor(ifelse(d==0, paste0(ref,ref),
                ifelse(d==1, paste0(ref,alt),
                       ifelse(d==2, paste0(alt,alt), NA))),
         levels=c(paste0(ref,ref),paste0(ref,alt),paste0(alt,alt)))
}

get_genotype_df <- function(snp){
  row <- geno_dt[.(snp)]
  if (nrow(row) == 0) return(NULL)
  
  ref <- row$REF[1]
  alt <- row$ALT[1]
  
  # Extract numeric dosage vector correctly
  dos <- unlist(row[, ..keep], use.names = TRUE)
  dos <- as.numeric(dos)   # convert "0","1","2" to numeric
  names(dos) <- keep       # make sure donor names are preserved
  
  tibble(
    donor = names(dos),
    genotype = dosage_to_label(dos, ref, alt),
    dosage = dos
  )
}


# ---------------------- Stats ---------------------- ####
message("[4] Loading eQTL and reQTL modeling stats...")

# Helper: safe empty stats table with expected columns
make_empty_stats <- function() {
  data.table(
    snp = character(),
    gene = character(),
    celltype = character(),
    dist = numeric(),
    eQTL_beta_PBS = numeric(), eQTL_p_PBS = numeric(), eQTL_se_PBS = numeric(),
    eQTL_beta_IFNG = numeric(), eQTL_p_IFNG = numeric(), eQTL_se_IFNG = numeric(),
    eQTL_beta_IFNB = numeric(), eQTL_p_IFNB = numeric(), eQTL_se_IFNB = numeric(),
    eQTL_beta_TNF  = numeric(), eQTL_p_TNF  = numeric(), eQTL_se_TNF  = numeric(),
    reQTL_dbeta_IFNG = numeric(), reQTL_p_IFNG = numeric(), reQTL_se_IFNG = numeric(),
    reQTL_dbeta_IFNB = numeric(), reQTL_p_IFNB = numeric(), reQTL_se_IFNB = numeric(),
    reQTL_dbeta_TNF  = numeric(), reQTL_p_TNF  = numeric(), reQTL_se_TNF  = numeric()
  )
}

# Load stats if available; otherwise use empty placeholder table
stats <- NULL
if (!is.null(STATS_FILE) && file.exists(STATS_FILE)) {
  stats <- fread(STATS_FILE)
} else {
  stats <- make_empty_stats()
}

# Robust getters: always return beta/p/se (all NA if missing)
get_eQTL_stats <- function(s, g, ct, stim){
  if (is.null(stats) || !is.data.frame(stats) || nrow(stats) == 0) {
    return(list(beta = NA_real_, p = NA_real_, se = NA_real_))
  }
  row <- stats %>% filter(snp == s, gene == g, celltype == ct)
  if (nrow(row) == 0) return(list(beta = NA_real_, p = NA_real_, se = NA_real_))
  
  beta_col <- paste0("eQTL_beta_", stim)
  p_col    <- paste0("eQTL_p_", stim)
  se_col   <- paste0("eQTL_se_", stim)
  
  if (!all(c(beta_col, p_col, se_col) %in% names(row))) {
    return(list(beta = NA_real_, p = NA_real_, se = NA_real_))
  }
  
  list(
    beta = row[[beta_col]][1],
    p    = row[[p_col]][1],
    se   = row[[se_col]][1]
  )
}

get_reQTL_stats <- function(s, g, ct, stim){
  if (is.null(stats) || !is.data.frame(stats) || nrow(stats) == 0) {
    return(list(beta = NA_real_, p = NA_real_, se = NA_real_))
  }
  row <- stats %>% filter(snp == s, gene == g, celltype == ct)
  if (nrow(row) == 0) return(list(beta = NA_real_, p = NA_real_, se = NA_real_))
  
  beta_col <- paste0("reQTL_dbeta_", stim)
  p_col    <- paste0("reQTL_p_", stim)
  se_col   <- paste0("reQTL_se_", stim)
  
  if (!all(c(beta_col, p_col, se_col) %in% names(row))) {
    return(list(beta = NA_real_, p = NA_real_, se = NA_real_))
  }
  
  list(
    beta = row[[beta_col]][1],
    p    = row[[p_col]][1],
    se   = row[[se_col]][1]
  )
}

# Robust formatter: NA placeholder if missing
fmt_stat <- function(beta, p){
  if (any(is.na(c(beta, p)))) return("p=NA, β=NA")
  stars <- if (p < 1e-5) "***" else if (p < 1e-4) "**" else if (p < 1e-3) "*" else ""
  paste0("p=", formatC(p, format = "e", digits = 1),
         ", β=", sprintf("%.1f", beta),
         if (stars != "") paste0(", ", stars) else "")
}

# Robust distance getter (per pair)
get_dist <- function(s, g, ct){
  if (is.null(stats) || !is.data.frame(stats) || nrow(stats) == 0 || !"dist" %in% names(stats)) return(NA_real_)
  row <- stats %>% filter(snp == s, gene == g, celltype == ct)
  if (nrow(row) == 0) return(NA_real_)
  as.numeric(row$dist[1])
}

# ---------------------- Plot Function ----------------------
# ---------------------- Plot Helper (PBS vs condition panel) ----------------------
make_subpanel <- function(df, ct, stim, this_snp, this_gene){
  
  sub <- df %>% filter(celltype == ct, condition %in% c("PBS", stim))
  if(nrow(sub) == 0) return(NULL)
  
  # Count donors per genotype + condition
  counts <- sub %>%
    group_by(genotype, dosage, condition) %>%
    summarize(n = n_distinct(donor), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
    mutate(strip_label = paste0(genotype, "\nPBS=", PBS, " ", stim, "=", !!sym(stim)))
  
  # join counts into sub
  sub <- sub %>%
    left_join(counts, by = c("genotype", "dosage"))
  
  # ---- enforce facet order by numeric dosage ----
  sub <- sub %>%
    arrange(dosage) %>%   # ensures 0 → 1 → 2 before converting to factor
    mutate(
      genotype_label = factor(strip_label, levels = unique(strip_label))
    )
  
  # ensure PBS appears left on x-axis
  sub$condition <- factor(sub$condition, levels = c("PBS", stim))
  
  # ---- plotting ----
  st.reQTL <- get_reQTL_stats(this_snp, this_gene, ct, stim)
  st.eQTL.PBS <- get_eQTL_stats(this_snp, this_gene, ct, "PBS")
  st.eQTL.stim <- get_eQTL_stats(this_snp, this_gene, ct, stim)
  
  subtitle <- paste0(ct, " — PBS vs ", stim, "\n",
                     "(reQTL ", stim, ": ", fmt_stat(st.reQTL$beta, st.reQTL$p), ")","\n",
                     "(eQTL PBS", ": ", fmt_stat(st.eQTL.PBS$beta, st.eQTL.PBS$p), ")", "\n",
                     "(eQTL ",stim,": ", fmt_stat(st.eQTL.stim$beta, st.eQTL.stim$p), ")")
  
  
  ggplot(sub, aes(x = condition, y = CPM, group = donor)) +
    geom_line(alpha = 0.25, color = "grey70") +
    geom_point(aes(color = condition), size = 2.2, alpha = 0.9) +
    scale_color_manual(values = c(
      "PBS" = "gray35",
      "IFNG" = "#0072B2",
      "IFNB" = "#E69F00",
      "TNF" = "#D55E00"
    )) +
    facet_wrap(~ genotype_label, scales = "fixed", nrow = 1) +
    labs(title = subtitle, y = "CPM", x = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 3),
      strip.text = element_text(size = 8, lineheight = 0.8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.major.x = element_blank(),   # already removing vertical major grid
      panel.grid.major.y = element_blank(),   # ← removes horizontal major grid
      panel.grid.minor.y = element_blank()    # ← removes horizontal minor grid
    )
  
}

# ---------------------- Main Plot for a SNP : gene ----------------------
make_plot <- function(this_gene, this_snp){
  
  df <- expr_long %>% filter(gene_id == this_gene)
  geno <- get_genotype_df(this_snp)
  df <- df %>% inner_join(geno, by="donor")
  
  ct_order <- c("FRB","KRT","MEL")
  stim_order <- c("IFNG","IFNB","TNF")
  
  panels <- list()
  for(ct in ct_order){
    for(stim in stim_order){
      p <- make_subpanel(df, ct, stim, this_snp, this_gene)
      if(!is.null(p)) panels[[paste(ct, stim)]] <- p
    }
  }
  
  combined <- wrap_plots(panels, ncol=3)
  
  # Apply global title WITHOUT overriding sub-panel title sizes
  # SNP coordinate (safe)
  snp_row <- geno_dt[geno_dt$ID == this_snp, c("CHROM","POS"), drop = FALSE]
  snp_coord <- if (nrow(snp_row) == 1) paste0(snp_row$CHROM[1], ":", snp_row$POS[1]) else "NA:NA"
  
  # Distance (safe; uses ct passed into script, not the loop variable)
  dist_bp <- get_dist(this_snp, this_gene, ct)
  
  combined <- combined + plot_annotation(
    title = paste0(this_gene, ":", this_snp,
                   " [dist ", ifelse(is.na(dist_bp), "NA", as.integer(dist_bp)), " bp]",
                   " [SNP ", snp_coord, "]")
  )
  combined & theme(
    plot.title = element_text(size=7, face="bold", hjust=0.5)  # affects ONLY main title now
  )
}

# ---------------------- Output ----------------------
message("[5] Writing PDF: ", OUT_PDF)
pdf(OUT_PDF, width=16, height=12)
print(make_plot(g, snp))
dev.off()
message("Done")


