## ============================================================
##  ATAC-seq Preflight Script
##  Run this BEFORE atac_analysis.R
##
##  Outputs:
##    1. sample_metadata_with_libsize.tsv  <- feed into main script
##    2. PDF report with count distributions to pick CPM threshold
## ============================================================

suppressPackageStartupMessages({
  library(Rsamtools)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

## ============================================================
##  PARAMETERS
## ============================================================

BAM_DIR   <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/bam_dedupped"
METADATA  <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/metadata/sample_metadata.txt"
MULTICOV  <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks_multicov.txt"
OUT_DIR   <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/results/preflight"

# metadata must already have columns: sample, donor, celltype, condition
# (no total_reads yet — that's what we compute here)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ============================================================
##  1. COMPUTE LIBRARY SIZES FROM BAMs
##     Uses idxstats — counts mapped reads per chromosome,
##     sums to total mapped reads. Requires .bai index files.
## ============================================================

cat("=== Step 1: Computing library sizes from BAMs ===\n")

bam_files <- sort(list.files(BAM_DIR, pattern = "\\.nodup\\.bam$", full.names = TRUE))
cat("Found", length(bam_files), "BAM files\n\n")

lib_sizes <- vapply(bam_files, function(bam) {
  sample_name <- sub("\\.nodup\\.bam$", "", basename(bam))
  cat("  Reading:", sample_name, "... ")

  # Check index exists
  bai <- paste0(bam, ".bai")
  if (!file.exists(bai)) {
    cat("indexing... ")
    indexBam(bam)
  }

  # idxstatsBam returns a df with seqnames, seqlength, mapped, unmapped
  idx   <- idxstatsBam(bam)
  total <- sum(idx$mapped)
  cat(format(total, big.mark = ","), "mapped reads\n")
  total
}, numeric(1))

lib_df <- tibble(
  sample      = sub("\\.nodup\\.bam$", "", basename(bam_files)),
  total_reads = lib_sizes
)

cat("\nLibrary size summary:\n")
cat(sprintf("  Min:    %s\n", format(min(lib_sizes),    big.mark = ",")))
cat(sprintf("  Median: %s\n", format(median(lib_sizes), big.mark = ",")))
cat(sprintf("  Max:    %s\n", format(max(lib_sizes),    big.mark = ",")))

## ============================================================
##  2. MERGE WITH METADATA AND SAVE
## ============================================================

cat("\n=== Step 2: Merging with metadata ===\n")

meta_in <- read_tsv(METADATA, col_types = cols())

# Check required columns
required_cols <- c("sample", "donor", "celltype", "condition")
missing <- setdiff(required_cols, colnames(meta_in))
if (length(missing) > 0) stop("Metadata missing columns: ", paste(missing, collapse = ", "))

# Check all BAM samples are in metadata
missing_in_meta <- setdiff(lib_df$sample, meta_in$sample)
if (length(missing_in_meta) > 0) {
  warning("These BAM samples are NOT in metadata:\n  ",
          paste(missing_in_meta, collapse = "\n  "))
}

meta_out <- meta_in %>%
  left_join(lib_df, by = "sample") %>%
  mutate(
    condition = factor(condition, levels = c("PBS", "IFNG")),
    celltype  = factor(celltype,  levels = c("FRB", "KRT", "MEL")),
    donor     = factor(donor)
  )

out_meta <- file.path(OUT_DIR, "sample_metadata_with_libsize.tsv")
write_tsv(meta_out, out_meta)
cat("Saved:", out_meta, "\n")

## ============================================================
##  3. LOAD COUNT MATRIX AND COMPUTE CPM
## ============================================================

cat("\n=== Step 3: Loading count matrix and computing CPM ===\n")

multicov  <- read_tsv(MULTICOV, col_types = cols())
count_mat <- as.matrix(multicov[, 5:ncol(multicov)])
colnames(count_mat) <- colnames(multicov)[5:ncol(multicov)]
peak_names <- multicov[[4]]
rownames(count_mat) <- peak_names

# Align lib_sizes to count matrix column order
meta_aligned <- meta_out[match(colnames(count_mat), meta_out$sample), ]
stopifnot(all(meta_aligned$sample == colnames(count_mat)))

lib_sizes_aligned <- meta_aligned$total_reads
names(lib_sizes_aligned) <- meta_aligned$sample

cpm_mat <- sweep(count_mat, 2, lib_sizes_aligned / 1e6, FUN = "/")

cat(sprintf("Count matrix: %d peaks × %d samples\n",
            nrow(count_mat), ncol(count_mat)))

## ============================================================
##  4. EMPIRICAL DISTRIBUTION PLOTS
##     Goal: help you pick CPM threshold for open/closed calls
## ============================================================

cat("\n=== Step 4: Generating empirical distribution report ===\n")

CELLTYPES  <- levels(meta_aligned$celltype)
CONDITIONS <- levels(meta_aligned$condition)
ct_colors  <- c(FRB = "#E69F00", KRT = "#56B4E9", MEL = "#009E73")

pdf(file.path(OUT_DIR, "preflight_report.pdf"), width = 10, height = 7)

# ---- Plot 1: Library sizes by sample ----
p1 <- meta_out %>%
  mutate(sample = factor(sample, levels = sample[order(celltype, condition, donor)])) %>%
  ggplot(aes(x = sample, y = total_reads / 1e6, fill = celltype, alpha = condition)) +
    geom_col() +
    scale_fill_manual(values = ct_colors) +
    scale_alpha_manual(values = c(PBS = 0.55, IFNG = 1.0)) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Library sizes (mapped reads, deduplicated)",
         x = NULL, y = "Mapped reads (millions)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
print(p1)

# ---- Plot 2: Global CPM distribution (log10) ----
# Subsample peaks for speed (max 50k)
set.seed(42)
idx_sample <- sample(nrow(cpm_mat), min(50000, nrow(cpm_mat)))
cpm_long <- cpm_mat[idx_sample, ] %>%
  as.data.frame() %>%
  rownames_to_column("peak_name") %>%
  pivot_longer(-peak_name, names_to = "sample", values_to = "CPM") %>%
  left_join(meta_aligned %>% select(sample, celltype, condition, donor),
            by = "sample")

p2 <- ggplot(cpm_long, aes(x = log10(CPM + 0.01), color = celltype, group = sample)) +
  geom_density(alpha = 0.7, linewidth = 0.4) +
  scale_color_manual(values = ct_colors) +
  geom_vline(xintercept = log10(c(0.5, 1, 2)), lty = 2,
             color = c("grey50","black","grey50"), linewidth = 0.5) +
  annotate("text", x = log10(0.5), y = Inf, label = "CPM=0.5",
           hjust = -0.1, vjust = 1.5, size = 3, color = "grey50") +
  annotate("text", x = log10(1),   y = Inf, label = "CPM=1",
           hjust = -0.1, vjust = 1.5, size = 3, color = "black") +
  annotate("text", x = log10(2),   y = Inf, label = "CPM=2",
           hjust = -0.1, vjust = 1.5, size = 3, color = "grey50") +
  labs(title = "Global CPM distribution (50k peaks sampled)",
       subtitle = "Dashed lines = candidate thresholds. Aim to cut at the antimode (valley between two peaks).",
       x = "log10(CPM + 0.01)", y = "Density") +
  theme_bw(base_size = 11)
print(p2)

# ---- Plot 3: Per-celltype CPM density, PBS vs IFNG ----
for (ct in CELLTYPES) {
  p3 <- cpm_long %>%
    filter(celltype == ct) %>%
    ggplot(aes(x = log10(CPM + 0.01), color = condition, group = sample)) +
      geom_density(linewidth = 0.5) +
      scale_color_manual(values = c(PBS = "steelblue", IFNG = "firebrick")) +
      geom_vline(xintercept = log10(1), lty = 2, color = "black") +
      labs(title = paste0(ct, ": CPM distribution per sample"),
           subtitle = "Each line = one sample. Dashed = CPM 1.",
           x = "log10(CPM + 0.01)", y = "Density") +
      theme_bw(base_size = 11)
  print(p3)
}

# ---- Plot 4: % peaks open vs CPM threshold sweep (per celltype × condition) ----
thresholds <- c(0.25, 0.5, 1, 2, 3, 5)
min_donors_grid <- c(2, 3, 4, 5)

sweep_rows <- list()
for (ct in CELLTYPES) {
  for (cond in CONDITIONS) {
    samps <- meta_aligned$sample[meta_aligned$celltype == ct &
                                 meta_aligned$condition == cond]
    sub   <- cpm_mat[, samps, drop = FALSE]
    for (thr in thresholds) {
      donors_pass <- rowSums(sub > thr)
      for (md in min_donors_grid) {
        pct_open <- mean(donors_pass >= md) * 100
        sweep_rows[[length(sweep_rows) + 1]] <- tibble(
          celltype   = ct,
          condition  = cond,
          threshold  = thr,
          min_donors = md,
          pct_open   = pct_open
        )
      }
    }
  }
}
sweep_df <- bind_rows(sweep_rows) %>%
  mutate(label = paste0("≥", min_donors, " donors"))

p4 <- ggplot(sweep_df,
             aes(x = threshold, y = pct_open,
                 color = celltype, linetype = condition,
                 group = interaction(celltype, condition))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = ct_colors) +
  scale_x_continuous(breaks = thresholds) +
  facet_wrap(~ label, nrow = 2) +
  labs(title = "% peaks called 'open' vs CPM threshold",
       subtitle = "Facets = minimum donor requirement. Pick threshold where curves are stable.",
       x = "CPM threshold", y = "% peaks called open") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")
print(p4)

# ---- Plot 5: Fraction of peaks with exactly 0 counts ----
pct_zero <- colMeans(count_mat == 0) * 100
zero_df  <- tibble(sample = names(pct_zero), pct_zero = pct_zero) %>%
  left_join(meta_aligned %>% select(sample, celltype, condition, donor),
            by = "sample")

p5 <- ggplot(zero_df,
             aes(x = reorder(sample, pct_zero), y = pct_zero,
                 fill = celltype, alpha = condition)) +
  geom_col() +
  scale_fill_manual(values = ct_colors) +
  scale_alpha_manual(values = c(PBS = 0.55, IFNG = 1.0)) +
  labs(title = "% peaks with zero raw counts per sample",
       subtitle = "High % in one sample = possible QC failure",
       x = NULL, y = "% zero-count peaks") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
print(p5)

# ---- Plot 6: Donor-level consistency check ----
# For each peak, how variable are counts across donors within a celltype × condition?
# Use coefficient of variation (CV) — high CV = noisy peak
cv_rows <- list()
for (ct in CELLTYPES) {
  for (cond in CONDITIONS) {
    samps <- meta_aligned$sample[meta_aligned$celltype == ct &
                                 meta_aligned$condition == cond]
    sub   <- count_mat[, samps, drop = FALSE]
    # only on non-zero peaks
    m     <- rowMeans(sub)
    s     <- apply(sub, 1, sd)
    cv    <- s / (m + 1)   # +1 pseudocount to avoid div/0
    cv_rows[[length(cv_rows) + 1]] <- tibble(
      celltype  = ct,
      condition = cond,
      cv        = cv[is.finite(cv)]
    )
  }
}
cv_df <- bind_rows(cv_rows)

p6 <- ggplot(cv_df, aes(x = cv, color = celltype, linetype = condition)) +
  geom_density(linewidth = 0.7) +
  scale_color_manual(values = ct_colors) +
  coord_cartesian(xlim = c(0, 5)) +
  labs(title = "Within-group peak count variability (CV)",
       subtitle = "High CV = noisy peaks. Expect most peaks to have CV < 1.",
       x = "Coefficient of variation (raw counts)", y = "Density") +
  theme_bw(base_size = 11)
print(p6)

dev.off()

## ============================================================
##  5. PRINT DECISION GUIDE
## ============================================================

cat("\n=== Preflight complete ===\n\n")
cat("Files saved to:", OUT_DIR, "\n\n")

cat("------------------------------------------------------------\n")
cat("HOW TO CHOOSE YOUR CPM THRESHOLD\n")
cat("------------------------------------------------------------\n")
cat("Open preflight_report.pdf and look at:\n\n")
cat("  Plot 2 (global density): find the antimode — the valley\n")
cat("    between the closed-peak pile and the open-peak pile.\n")
cat("    That valley is your natural CPM cutoff. Typically 0.5–2.\n\n")
cat("  Plot 3 (per-celltype): check that all samples behave\n")
cat("    similarly. An outlier sample will have a shifted curve.\n\n")
cat("  Plot 4 (threshold sweep): pick a threshold where the\n")
cat("    '% open' curve flattens — meaning you're past the noise\n")
cat("    floor. Combine with a min_donors of 4 (out of 6) for\n")
cat("    reproducibility.\n\n")
cat("  Plot 5 (zero counts): any sample with >60% zero-count\n")
cat("    peaks is likely a QC failure — consider removing it.\n\n")
cat("Then set in atac_analysis.R:\n")
cat("  CPM_THRESHOLD   <- <your chosen value>\n")
cat("  MIN_DONORS_OPEN <- 4  # or adjust based on Plot 4\n")
cat("  METADATA        <- \"", out_meta, "\"\n", sep="")
cat("------------------------------------------------------------\n")
