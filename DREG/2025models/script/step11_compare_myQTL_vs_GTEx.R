#!/usr/bin/env Rscript
# ============================================================
# compare_eQTLs.R
#
# Merges my eQTL results with GTEx fibroblast eQTLs.
#
# Strategy:
#   1. Load GTEx sig.eQTLs.bed  → rsID + gene_symbol pairs
#   2. Load GTEx signif_variant_gene_pairs.txt.gz → slope (beta) + pval
#      and cross-reference via chr_pos_ref_alt → rsID using the bed file
#   3. Load my genotype_pairs_chunk_*.tsv files → SNP info (rsID, REF, ALT)
#      NOTE: my chunk files only contain GENOTYPE dosage, NOT my eQTL
#      results. You must supply my own eQTL summary stats (see MY_EQTL_FILE
#      below). The script will still work for harmonization once you
#      provide that file.
#   4. Harmonize alleles: if my REF/ALT is the flip of GTEx, negate beta.
#   5. Output a wide comparison table.
#
# Required inputs (edit the paths section below):
#   - GTEX_DIR      : directory containing the three GTEx files
#   - MY_EQTL_FILE  : my eQTL summary stats (see expected columns below)
#   - CHUNKS_DIR    : directory with genotype_pairs_chunk_*.tsv (for SNP lookup)
#   - OUT_FILE      : output TSV path
#
# Expected columns in MY_EQTL_FILE:
#   snp         – rsID  (e.g. rs12563495)
#   gene        – gene symbol  (e.g. WASH7P)
#   beta        – effect size (same direction as slope)
#   pvalue      – nominal p-value
#   REF         – reference allele
#   ALT         – alternative allele
#   (any other columns are carried through)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(gtexr)
})

# ── PATHS ───────────────────────────────────────────────
GTEX_DIR     <- "/pi/manuel.garber-umw/human/skin/eQTLs/literature/GTEx/GTEx_Analysis_v8_eQTL"
GTEX_BEDFILE <- "Cells_Cultured_fibroblasts.v8.sig.eQTLs.bed"
GTEX_FILE <- "Cells_Cultured_fibroblasts.v8.signif_variant_gene_pairs.txt.gz"
GTEX_EGENES_FILE <- "Cells_Cultured_fibroblasts.v8.egenes.txt.gz"
CHUNKS_DIR   <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/chunks"
MY_EQTL_DIR <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/results_QC"
MY_EQTL_TYPE <- "FRB_IFNG_eQTL"
OUT_DIR <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data"
OUT_FILE     <- paste0(OUT_DIR,"/GTEx_fibroblasts_comparison_my_FRB_IFNG_eQTL.tsv")
OUT_PLOT     <- paste0(OUT_DIR,"/GTEx_fibroblasts_comparison_my_FRB_IFNG_eQTL.pdf")
myqval = fread("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/eigenMT/results/FRB_IFNG_eQTL.eigenMT.txt")
p_threshold = max(myqval[which(myqval$q_gene<0.05, )]$pmin)

# ───────Step 1: Loading GTEx sig.eQTLs.bed───────────────────────── ####
cat("=== Step 1: Loading GTEx sig.eQTLs.bed (rsID + gene symbol) ===\n")
bed_cols <- c("chr","start","end","name","rsid","gene_symbol",
              "snp_coord","gene_coord","tss_distance")
gtex_bed <- fread(
  file.path(GTEX_DIR, GTEX_BEDFILE),
  col.names = bed_cols, header = FALSE
)
# Keep only the rsID → gene_symbol mapping (one row per pair)
gtex_pairs <- gtex_bed[, .(rsid, gene_symbol, chr, start)] %>% unique()
cat(sprintf("  %d unique rsID-gene pairs in BED\n", nrow(gtex_pairs)))

# ──────────────────────────────────────────────────────────── ####
cat("\n=== Step 2: Loading GTEx signif_variant_gene_pairs (slope/beta) ===\n")
gtex_sig <- fread(
  cmd = paste("zcat",
              file.path(GTEX_DIR,GTEX_FILE)),
  header = TRUE
)
# Parse variant_id: chr_pos_ref_alt_b38
gtex_sig[, c("v_chr","v_pos","gtex_ref","gtex_alt","build") :=
           tstrsplit(variant_id, "_", fixed = TRUE)]
gtex_sig[, v_pos := as.integer(v_pos)]

# ──────────────────────────────────────────────────────────── ####
# Load the egenes file just to get a gene_id → gene_name map
cat("\n=== Step 3: gene_id → gene_symbol map from egenes ===\n")
egenes <- fread(
  cmd = paste("zcat",
              file.path(GTEX_DIR,GTEX_EGENES_FILE)),
  header = TRUE,
  select = c("gene_id","gene_name")
) %>% unique()
setnames(egenes, c("gene_id","gene_name"), c("gene_id","gene_symbol_eg"))

# Join gene symbol onto sig pairs
gtex_sig <- merge(gtex_sig, egenes,
                  by.x = "gene_id", by.y = "gene_id", all.x = TRUE)

# Now cross-reference rsID via the BED file
# BED has chr + start (0-based) → v_pos is 1-based in GTEx variant_id
gtex_bed2 <- gtex_bed[, .(rsid, chr, pos1 = start)]   # convert to 1-based
gtex_sig <- merge(gtex_sig,
                  gtex_bed2[, .(v_chr = chr, v_pos = pos1, rsid)] %>% unique(),
                  by = c("v_chr","v_pos"), all.x = TRUE)

# Use gene symbol from egenes (more complete) but fall back to bed
# Final GTEx table: one row per rsid-gene_symbol pair with slope + pval
gtex_final <- gtex_sig[, .(
  rsid,
  gtex_gene   = gene_symbol_eg,
  gtex_ref,
  gtex_alt,
  gtex_slope  = slope,
  gtex_slope_se = slope_se,
  gtex_pval   = pval_nominal,
  gtex_pval_beta = pval_beta,
  gtex_variant_id = variant_id
)] %>% unique()

cat(sprintf("  %d GTEx significant variant-gene pairs with rsID\n",
            nrow(gtex_final[!is.na(rsid)])))
rm(gtex_sig, gtex_bed, gtex_bed2, gtex_pairs)

# ──────────────────────────────────────────────────────────── ####
cat("\n=== Step 4: Build SNP coordinate table from genotype chunks ===\n")
# The chunk files give us: rsID (ID column) + REF + ALT
# We only need the first 5 columns; read them efficiently
chunk_files <- list.files(CHUNKS_DIR, pattern = "genotype_pairs_chunk_.*\\.tsv$",
                          full.names = TRUE)
cat(sprintf("  Found %d chunk files\n", length(chunk_files)))

snp_info <- rbindlist(lapply(chunk_files, function(f) {
  fread(f, header = TRUE, select = c("CHROM","POS","ID","REF","ALT"))
})) %>% unique(by = "ID")

setnames(snp_info, c("CHROM","POS","ID","REF","ALT"),
         c("chrom","pos","rsid","my_ref","my_alt"))
cat(sprintf("  %d unique SNPs in genotype chunks\n", nrow(snp_info)))

# ──────────────────────────────────────────────────────────── ####
cat("\n=== Step 5: Load my eQTL summary stats ===\n")
myqtl_files <- list.files(MY_EQTL_DIR, pattern = paste0("modeling_stats_postQC_",MY_EQTL_TYPE,"_.*\\.txt$"),
                          full.names = TRUE)
cat(sprintf("  Found %d my qtl files\n", length(myqtl_files)))

my_eqtl <- rbindlist(lapply(myqtl_files, function(f) {
  fread(f, header = TRUE)
}))

# If REF/ALT not in my eQTL file, pull from chunk SNP info
if (!"ref" %in% names(my_eqtl) || !"alt" %in% names(my_eqtl)) {
  cat("  REF/ALT not found in my eQTL file – merging from genotype chunks\n")
  my_eqtl <- merge(my_eqtl,
                   snp_info[, .(rsid, my_ref, my_alt)],
                   by.x = "snp", by.y = "rsid", all.x = TRUE)
} else {
  setnames(my_eqtl, c("ref","alt"), c("my_ref","my_alt"))
}

cat(sprintf("  %d rows in my eQTL table\n", nrow(my_eqtl)))

# ──────────────────────────────────────────────────────────── ####
cat("\n=== Step 6: Merge on rsid + gene_symbol ===\n")
merged <- merge(
  my_eqtl,
  gtex_final[!is.na(rsid)],
  by.x = c("snp","gene"),
  by.y = c("rsid","gtex_gene"),
  all = FALSE    # inner join – keep only SNP-gene pairs in BOTH datasets
)
cat(sprintf("  %d overlapping SNP-gene pairs\n", nrow(merged)))

# ──────────────────────────────────────────────────────────── ####
cat("\n=== Step 7: Allele harmonization & beta flip ===\n")
# Possible situations (ignoring strand for now):
#   1. my REF == GTEx REF and my ALT == GTEx ALT  → concordant, no flip
#   2. my REF == GTEx ALT and my ALT == GTEx REF  → swapped, flip beta sign
#   3. anything else                               → ambiguous, flag

merged[, allele_status := fcase(
  my_ref == gtex_ref & my_alt == gtex_alt, "concordant",
  my_ref == gtex_alt & my_alt == gtex_ref, "flipped",
  is.na(gtex_alt), "gtex_not_found",
  default = "ambiguous"
)]

# For ambiguous cases, also try complement matching (A↔T, C↔G)
complement <- function(x) chartr("ACGT","TGCA", x)
merged[allele_status == "ambiguous",
       allele_status := fcase(
         complement(my_ref) == gtex_ref & complement(my_alt) == gtex_alt, "complement",
         complement(my_ref) == gtex_alt & complement(my_alt) == gtex_ref, "complement_flipped",
         default = "unresolved"
       )]

# Apply beta flip where needed
merged[, my_beta_harmonized := fcase(
  allele_status %in% c("concordant","complement"),         beta,
  allele_status %in% c("flipped","complement_flipped"),   -beta,
  default = NA_real_   # unresolved – exclude from downstream
)]

table(merged$allele_status)

# ---------- add another column of whether the p.nominal in my QTL is significant --------- ####
merged[, my_signif := fcase(
  p < p_threshold,         "sig",
  p >= p_threshold,         "not_sig"
)]

temp <- merged %>%
  dplyr::filter(my_signif=="sig" & allele_status != "unresolved") 
  
table(temp$allele_status)


# ──────────────────────────────────────────────────────────── ####
cat("\n=== Step 8: Write output ===\n")
# Reorder columns for clarity
out_cols <- c(
  "snp","gene",
  "my_ref","my_alt","gtex_ref","gtex_alt","allele_status","my_signif",
  "beta","my_beta_harmonized",
  "gtex_slope","gtex_slope_se",
  "p","gtex_pval","gtex_pval_beta",
  "gtex_variant_id"
)
# Add any extra columns from my data that aren't already listed
extra <- setdiff(names(merged), out_cols)
out <- merged[, c(out_cols[out_cols %in% names(merged)], extra), with = FALSE]

fwrite(out, OUT_FILE, sep = "\t", quote = FALSE)
cat(sprintf("  Written: %s  (%d rows, %d cols)\n",
            OUT_FILE, nrow(out), ncol(out)))

cat("\n=== Quick summary ===\n")
cat(sprintf("  Total harmonized pairs      : %d\n", nrow(out)))
cat(sprintf("  Concordant alleles          : %d\n", sum(out$allele_status == "concordant", na.rm=TRUE)))
cat(sprintf("  Flipped (beta negated)      : %d\n", sum(out$allele_status == "flipped",    na.rm=TRUE)))
cat(sprintf("  Complement-flipped          : %d\n", sum(out$allele_status %in% c("complement","complement_flipped"), na.rm=TRUE)))
cat(sprintf("  Unresolved (excluded)       : %d\n", sum(out$allele_status == "unresolved",  na.rm=TRUE)))
cat(sprintf("\n  Recommended comparison metric: my_beta_harmonized vs gtex_slope\n"))
cat(sprintf("  (Pearson correlation on significant pairs, optionally weighted by 1/slope_se)\n"))
cat(sprintf("  Not found in GTEx           : %d\n", sum(out$allele_status == "gtex_not_found", na.rm=TRUE)))

# ============================================================ ####
# METRIC 1: Spearman correlation (all vs sig only)
dt <- out[allele_status != "unresolved"]
sig <- dt[dt$my_signif == "sig"]

r_all  <- cor(dt$my_beta_harmonized,  dt$gtex_slope,  method = "spearman", use = "complete.obs")
r_sig  <- cor(sig$my_beta_harmonized, sig$gtex_slope, method = "spearman", use = "complete.obs")
r_pear_all <- cor(dt$my_beta_harmonized,  dt$gtex_slope,  method = "pearson", use = "complete.obs")
r_pear_sig <- cor(sig$my_beta_harmonized, sig$gtex_slope, method = "pearson", use = "complete.obs")

cat(sprintf("\n--- Correlation ---\n"))
cat(sprintf("  Pearson  r  (all pairs)  : %.3f\n", r_pear_all))
cat(sprintf("  Pearson  r  (sig only)   : %.3f\n", r_pear_sig))
cat(sprintf("  Spearman rho (all pairs) : %.3f\n", r_all))
cat(sprintf("  Spearman rho (sig only)  : %.3f\n", r_sig))

# ============================================================ ####
# METRIC 2: Sign concordance
sign_all <- mean(sign(dt$my_beta)  == sign(dt$gtex_slope),  na.rm = TRUE)
sign_sig <- mean(sign(sig$my_beta) == sign(sig$gtex_slope), na.rm = TRUE)

cat(sprintf("\n--- Sign concordance ---\n"))
cat(sprintf("  All pairs : %.1f%%\n", 100 * sign_all))
cat(sprintf("  Sig only  : %.1f%%\n", 100 * sign_sig))

# Binomial test vs 50% (chance)
binom_all <- binom.test(sum(sign(dt$my_beta)  == sign(dt$gtex_slope),  na.rm=TRUE), nrow(dt),  p=0.5)
binom_sig <- binom.test(sum(sign(sig$my_beta) == sign(sig$gtex_slope), na.rm=TRUE), nrow(sig), p=0.5)
cat(sprintf("  Binomial p (all): %.2e  (sig only): %.2e\n",
            binom_all$p.value, binom_sig$p.value))
# ============================================================ ####
# METRIC 3: pi1 – fraction of My sig eQTLs that replicate in GTEx
#   pi1 = 1 - pi0  where pi0 estimated from p-value histogram
#   Uses Storey's method via qvalue package if available,
#   otherwise uses a simple lambda-based estimate
cat(sprintf("\n--- pi1 replication statistic ---\n"))
pi1_gtex_pvals <- sig$gtex_pval   # GTEx p-values for My significant pairs

pi1_val <- tryCatch({
  if (!requireNamespace("qvalue", quietly = TRUE)) stop("no qvalue")
  qobj <- qvalue::qvalue(p = pi1_gtex_pvals)
  1 - qobj$pi0
}, error = function(e) {
  # Simple lambda = 0.5 estimate as fallback
  lambda <- 0.5
  pi0_hat <- min(1, mean(pi1_gtex_pvals > lambda) / (1 - lambda))
  1 - pi0_hat
})

cat(sprintf("  pi1 (proportion replicating in GTEx): %.3f\n", pi1_val))
cat(sprintf("  (Install Bioconductor 'qvalue' package for more accurate pi1)\n"))

# Nominal replication rate (GTEx p < 0.05)
nominal_rep <- mean(sig$gtex_pval < 0.05, na.rm = TRUE)
cat(sprintf("  Nominal replication rate (GTEx p<0.05): %.1f%%\n", 100 * nominal_rep))

# ============================================================ ####
# PLOTS
cat("\nGenerating plots...\n")

# Color by -log10(My pval), capped for visual clarity
dt[, neg_log10_p := pmin(-log10(p), 10)]
sig[, neg_log10_p := pmin(-log10(p), 10)]

# ── Panel A: All pairs, colored by significance ──────────────
p_all <- ggplot(dt, aes(gtex_slope, my_beta_harmonized, colour = neg_log10_p)) +
  geom_point(alpha = 0.35, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed", linewidth = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "navy", linewidth = 1) +
  scale_colour_gradient(low = "grey75", high = "firebrick",
                        name = "-log10(p)\n(my data)") +
  labs(
    title = sprintf("All harmonized pairs  (Spearman rho=%.3f)",
                    r_all),
    x = "GTEx slope", y = "My beta (harmonized)"
  ) +
  theme_bw(base_size = 11)

# ── Panel B: My significant pairs only ─────────────────────
p_sig <- ggplot(sig, aes(gtex_slope, my_beta_harmonized, colour = neg_log10_p)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed", linewidth = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "navy", linewidth = 1) +
  scale_colour_gradient(low = "steelblue2", high = "firebrick",
                        name = "-log10(p)\n(my data)") +
  labs(
    title = sprintf("My sig pairs (p≤%.3g, n=%d)  Spearman rho=%.3f",
                    p_threshold, nrow(sig), r_sig),
    x = "GTEx slope", y = "My beta (harmonized)"
  ) +
  theme_bw(base_size = 11)

# ── Panel C: GTEx p-value histogram for My sig pairs (pi1 visual) ──
p_hist <- ggplot(sig, aes(x = gtex_pval)) +
  geom_histogram(bins = 20, fill = "steelblue", colour = "white", boundary = 0) +
  geom_hline(yintercept = nrow(sig) * 0.05, colour = "red",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.75, y = nrow(sig) * 0.05,
           label = "H0 baseline (uniform)", colour = "red", vjust = -0.5, size = 3.2) +
  labs(
    title = sprintf("GTEx p-val distribution for My sig pairs\npi1 = %.3f  |  Nominal rep rate = %.1f%%",
                    pi1_val, 100 * nominal_rep),
    x = "GTEx p-value", y = "Count"
  ) +
  theme_bw(base_size = 11)

# ── Panel D: Sign concordance bar ────────────────────────────
conc_df <- data.frame(
  subset   = c("All pairs", sprintf("Sig (p≤%.3g)", p_threshold)),
  concordant = c(sign_all, sign_sig) * 100,
  n        = c(nrow(dt), nrow(sig))
)
conc_df$label <- sprintf("%.1f%%\n(n=%d)", conc_df$concordant, conc_df$n)

p_sign <- ggplot(conc_df, aes(x = subset, y = concordant, fill = subset)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
  geom_text(aes(label = label), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
  labs(
    title = "Sign concordance with GTEx",
    subtitle = sprintf("Binomial test vs 50%%: p(all)=%.2e, p(sig)=%.2e",
                       binom_all$p.value, binom_sig$p.value),
    x = NULL, y = "% same direction as GTEx"
  ) +
  theme_bw(base_size = 11)

# ── Combine with patchwork or cowplot if available ───────────
combined <- tryCatch({
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("no patchwork")
  library(patchwork)
  (p_all | p_sig) / (p_hist | p_sign)
}, error = function(e) {
  tryCatch({
    if (!requireNamespace("cowplot", quietly = TRUE)) stop("no cowplot")
    library(cowplot)
    plot_grid(p_all, p_sig, p_hist, p_sign, ncol = 2, labels = "AUTO")
  }, error = function(e2) {
    # fallback: save individually
    ggsave("eQTL_panel_A_all.png",       p_all,   width=7, height=5, dpi=150)
    ggsave("eQTL_panel_B_sig.png",       p_sig,   width=7, height=5, dpi=150)
    ggsave("eQTL_panel_C_hist.png",      p_hist,  width=7, height=5, dpi=150)
    ggsave("eQTL_panel_D_sign.png",      p_sign,  width=6, height=5, dpi=150)
    cat("Saved 4 individual panel PNGs (patchwork/cowplot not available)\n")
    NULL
  })
})

if (!is.null(combined)) {
  ggsave(OUT_PLOT, combined, width = 14, height = 10, dpi = 150)
  cat(sprintf("Saved combined figure: %s\n", OUT_PLOT))
}

# ── Print summary table ──────────────────────────────────────
cat("\n========== REPLICATION SUMMARY ==========\n")
cat(sprintf("  %-40s %s\n", "Metric", "All pairs  |  Sig only"))
cat(sprintf("  %-40s %.3f      |  %.3f\n",    "Pearson r",          r_pear_all, r_pear_sig))
cat(sprintf("  %-40s %.3f      |  %.3f\n",    "Spearman rho",       r_all,      r_sig))
cat(sprintf("  %-40s %.1f%%     |  %.1f%%\n", "Sign concordance",   100*sign_all, 100*sign_sig))
cat(sprintf("  %-40s —          |  %.3f\n",   "pi1 (Storey)",       pi1_val))
cat(sprintf("  %-40s —          |  %.1f%%\n", "Nominal rep (GTEx p<0.05)", 100*nominal_rep))
cat("==========================================\n")



