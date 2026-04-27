## ============================================================
##  ATAC-seq Analysis Pipeline
##  Steps: load → normalize → open/closed calls →
##         differential chromatin accessibility → summary table
## ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(RColorBrewer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
})

## ============================================================
##  PARAMETERS — edit these
## ============================================================

MULTICOV      <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/master_peaks_multicov.txt"

# *** Use the file OUTPUT by atac_preflight.R (includes total_reads column) ***
METADATA      <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/results/preflight/sample_metadata_with_libsize.tsv"

OUT_DIR       <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/results"

# metadata must have columns: sample, donor, celltype, condition, total_reads
# condition values must be exactly: PBS, IFNG
CELLTYPES     <- c("FRB", "KRT", "MEL")
CONDITIONS    <- c("PBS", "IFNG")

# Open/closed calling
# *** Set CPM_THRESHOLD after reviewing preflight_report.pdf Plot 2 + Plot 4 ***
# Look for the antimode (valley) in the density plot — typically 0.5–2
CPM_THRESHOLD    <- 1       # REPLACE with your chosen value after preflight
MIN_DONORS_OPEN  <- 4       # peak must be open in >= this many donors (out of 6)

# DCA significance thresholds
PADJ_CUTOFF   <- 0.05
LFC_CUTOFF    <- 0.5        # |log2FoldChange| > this

# DESeq2 design: paired by donor
# formula: ~ donor + condition

## ============================================================
##  SETUP
## ============================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "DESeq2_objects"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "DCA_results"),    showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "QC_plots"),       showWarnings = FALSE)

cat("=== ATAC-seq Analysis Pipeline ===\n")
cat("Started:", format(Sys.time()), "\n\n")

## ============================================================
##  1. LOAD DATA
## ============================================================

cat("Loading multicov matrix...\n")
multicov <- read_tsv(MULTICOV, col_types = cols())

# First 4 columns are peak coordinates
peak_coords <- multicov[, 1:4]
colnames(peak_coords) <- c("chr", "start", "end", "peak_name")
count_mat   <- as.matrix(multicov[, 5:ncol(multicov)])
rownames(count_mat) <- peak_coords$peak_name

cat("  Peaks:", nrow(count_mat), "\n")
cat("  Samples:", ncol(count_mat), "\n")

cat("Loading metadata...\n")
meta <- read_tsv(METADATA, col_types = cols())
# Ensure sample column matches BAM-derived column names in multicov
stopifnot("sample"      %in% colnames(meta),
          "donor"       %in% colnames(meta),
          "celltype"    %in% colnames(meta),
          "condition"   %in% colnames(meta),
          "total_reads" %in% colnames(meta))

# Align metadata to count matrix column order
meta <- meta[match(colnames(count_mat), meta$sample), ]
stopifnot(all(meta$sample == colnames(count_mat)))  # hard check
rownames(meta) <- meta$sample

meta$condition <- factor(meta$condition, levels = c("PBS", "IFNG"))
meta$donor     <- factor(meta$donor)
meta$celltype  <- factor(meta$celltype, levels = CELLTYPES)

cat("  Donors:", nlevels(meta$donor), "\n")
cat("  Celltypes:", paste(levels(meta$celltype), collapse=", "), "\n\n")

## ============================================================
##  2. CPM NORMALIZATION (for open/closed calling only)
##     DESeq2 uses raw counts internally — do NOT pre-normalize
## ============================================================

cat("Computing CPM for open/closed calling...\n")
lib_sizes   <- meta$total_reads
names(lib_sizes) <- meta$sample

cpm_mat <- sweep(count_mat, 2, lib_sizes / 1e6, FUN = "/")

## ============================================================
##  3. OPEN / CLOSED CALLING PER CELLTYPE × CONDITION
## ============================================================

cat("Calling open/closed status per celltype × condition...\n")

# Returns "open" / "closed" for a set of samples (columns of cpm_mat)
# A peak is "open" if CPM > threshold in >= MIN_DONORS_OPEN donors
call_open_closed <- function(cpm_submat, threshold = CPM_THRESHOLD,
                             min_donors = MIN_DONORS_OPEN) {
  donors_passing <- rowSums(cpm_submat > threshold)
  ifelse(donors_passing >= min_donors, "open", "closed")
}

status_list <- list()

for (ct in CELLTYPES) {
  for (cond in CONDITIONS) {
    key     <- paste(ct, cond, sep = "_")
    samples <- meta$sample[meta$celltype == ct & meta$condition == cond]
    if (length(samples) == 0) {
      warning("No samples found for: ", key)
      next
    }
    status_list[[key]] <- call_open_closed(cpm_mat[, samples, drop = FALSE])
    n_open <- sum(status_list[[key]] == "open")
    cat(sprintf("  %s: %d open / %d closed (%.1f%% open)\n",
                key, n_open, nrow(cpm_mat) - n_open,
                100 * n_open / nrow(cpm_mat)))
  }
}

## ============================================================
##  4. QC PLOTS (all celltypes together)
## ============================================================

cat("\nGenerating QC plots...\n")

# 4a. Peak size distribution
peak_sizes <- peak_coords$end - peak_coords$start
pdf(file.path(OUT_DIR, "QC_plots", "peak_size_distribution.pdf"), width = 7, height = 5)
hist(peak_sizes, breaks = 100, main = "Peak size distribution",
     xlab = "Peak size (bp)", col = "steelblue", border = "white")
abline(v = median(peak_sizes), col = "red", lty = 2, lwd = 2)
legend("topright", legend = paste("Median:", median(peak_sizes), "bp"),
       col = "red", lty = 2, bty = "n")
dev.off()

# 4b. Library size bar plot
pdf(file.path(OUT_DIR, "QC_plots", "library_sizes.pdf"), width = 10, height = 5)
par(mar = c(8, 5, 3, 2))
cols_ct <- c(FRB = "#E69F00", KRT = "#56B4E9", MEL = "#009E73")
bar_cols <- cols_ct[as.character(meta$celltype)]
bp <- barplot(meta$total_reads / 1e6, names.arg = meta$sample,
              las = 2, col = bar_cols, cex.names = 0.6,
              ylab = "Total reads (millions)", main = "Library sizes")
legend("topright", legend = names(cols_ct), fill = cols_ct, bty = "n")
dev.off()

# 4c. PCA on VST-normalized counts (all samples)
cat("  Running PCA on all samples...\n")
# Build a quick DESeq2 object for vst
dds_all <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = ~ celltype + condition
)
# Filter low-count peaks before vst (keeps it fast)
keep_all <- rowSums(counts(dds_all) >= 10) >= 3
dds_all  <- dds_all[keep_all, ]
vsd_all  <- vst(dds_all, blind = TRUE)

pca_data <- plotPCA(vsd_all, intgroup = c("celltype", "condition"),
                    returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pdf(file.path(OUT_DIR, "QC_plots", "PCA_all_samples.pdf"), width = 7, height = 6)
print(
  ggplot(pca_data, aes(x = PC1, y = PC2,
                       color = celltype, shape = condition, label = name)) +
    geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20) +
    scale_color_manual(values = c(FRB = "#E69F00", KRT = "#56B4E9", MEL = "#009E73")) +
    labs(title = "PCA — all samples (VST)",
         x = paste0("PC1: ", pct_var[1], "% variance"),
         y = paste0("PC2: ", pct_var[2], "% variance")) +
    theme_bw(base_size = 12)
)
dev.off()

# 4d. Sample correlation heatmap
cor_mat <- cor(assay(vsd_all), method = "pearson")
ann_col <- data.frame(
  celltype  = meta$celltype,
  condition = meta$condition,
  row.names = meta$sample
)
ann_colors <- list(
  celltype  = c(FRB = "#E69F00", KRT = "#56B4E9", MEL = "#009E73"),
  condition = c(PBS = "grey70", IFNG = "firebrick")
)
pdf(file.path(OUT_DIR, "QC_plots", "sample_correlation_heatmap.pdf"), width = 12, height = 10)
pheatmap(cor_mat,
         annotation_col  = ann_col,
         annotation_row  = ann_col,
         annotation_colors = ann_colors,
         color           = colorRampPalette(c("white", "steelblue", "navy"))(100),
         main            = "Sample Pearson correlation (VST counts)",
         fontsize        = 7,
         show_rownames   = TRUE,
         show_colnames   = TRUE)
dev.off()

cat("  QC plots saved to", file.path(OUT_DIR, "QC_plots"), "\n\n")

## ============================================================
##  5. PEAK ANNOTATION (ChIPseeker)
## ============================================================

cat("Annotating peaks with ChIPseeker...\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_gr <- GRanges(
  seqnames = peak_coords$chr,
  ranges   = IRanges(start = peak_coords$start + 1,  # BED is 0-based
                     end   = peak_coords$end),
  peak_name = peak_coords$peak_name
)

anno <- annotatePeak(peak_gr, tssRegion = c(-2000, 200),
                     TxDb = txdb, annoDb = "org.Hs.eg.db",
                     verbose = FALSE)
anno_df <- as.data.frame(anno)
# Keep only the columns we care about
anno_df <- anno_df %>%
  select(peak_name, annotation, distanceToTSS,
         SYMBOL, GENENAME) %>%
  mutate(
    feature = case_when(
      grepl("Promoter",   annotation) ~ "promoter",
      grepl("Exon",       annotation) ~ "exon",
      grepl("Intron",     annotation) ~ "intron",
      grepl("Downstream", annotation) ~ "downstream",
      grepl("Intergenic", annotation) ~ "intergenic",
      TRUE                            ~ "other"
    )
  )

cat("  Annotation feature breakdown:\n")
print(table(anno_df$feature))
cat("\n")

## ============================================================
##  6. DIFFERENTIAL CHROMATIN ACCESSIBILITY (DESeq2)
##     One DESeq2 object per celltype; IFNG vs PBS, paired by donor
## ============================================================

cat("Running DESeq2 differential chromatin accessibility...\n\n")

dca_results  <- list()   # stores DESeq2 result tables
dds_objects  <- list()   # stores dds objects for saving

dynamics_list <- list()  # stores gainAcc / loseAcc / stable / NA per peak per celltype

for (ct in CELLTYPES) {
  cat(sprintf("--- %s ---\n", ct))

  # Subset to this celltype
  samples_ct <- meta$sample[meta$celltype == ct]
  meta_ct    <- droplevels(meta[samples_ct, ])
  counts_ct  <- count_mat[, samples_ct]

  # Build DESeq2 object — paired design
  dds <- DESeqDataSetFromMatrix(
    countData = counts_ct,
    colData   = meta_ct,
    design    = ~ donor + condition
  )
  dds$condition <- relevel(dds$condition, ref = "PBS")

  # Filter: keep peaks with >= 10 counts in >= half the samples
  keep <- rowSums(counts(dds) >= 10) >= (ncol(dds) / 2)
  dds  <- dds[keep, ]
  cat(sprintf("  Peaks passing filter: %d / %d\n", sum(keep), nrow(counts_ct)))

  # Run DESeq2
  dds  <- DESeq(dds, quiet = TRUE)

  # Extract IFNG vs PBS results
  res  <- results(dds,
                  contrast  = c("condition", "IFNG", "PBS"),
                  alpha     = PADJ_CUTOFF,
                  lfcThreshold = 0)

  # Shrink LFC (apeglm) for better effect size estimates
  res_shrunk <- lfcShrink(dds,
                          coef     = "condition_IFNG_vs_PBS",
                          type     = "apeglm",
                          quiet    = TRUE)

  # Summarize
  cat(sprintf("  DESeq2 summary (IFNG vs PBS):\n"))
  print(summary(res, alpha = PADJ_CUTOFF))

  # Dynamics calls (on shrunken LFC)
  res_df <- as.data.frame(res_shrunk) %>%
    rownames_to_column("peak_name") %>%
    mutate(
      dynamics = case_when(
        !is.na(padj) & padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "gainAcc",
        !is.na(padj) & padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "loseAcc",
        !is.na(padj)                                                       ~ "stable",
        TRUE                                                               ~ NA_character_
      )
    )

  cat(sprintf("  gainAcc: %d | loseAcc: %d | stable: %d | filtered/NA: %d\n",
              sum(res_df$dynamics == "gainAcc",  na.rm = TRUE),
              sum(res_df$dynamics == "loseAcc",  na.rm = TRUE),
              sum(res_df$dynamics == "stable",   na.rm = TRUE),
              sum(is.na(res_df$dynamics))))

  # Save full DCA table
  write_tsv(res_df,
            file.path(OUT_DIR, "DCA_results",
                      paste0(ct, "_IFNG_vs_PBS_DESeq2.tsv")))

  # MA plot
  pdf(file.path(OUT_DIR, "DCA_results", paste0(ct, "_MA_plot.pdf")),
      width = 7, height = 5)
  DESeq2::plotMA(res_shrunk, alpha = PADJ_CUTOFF,
                 main = paste0(ct, ": IFNG vs PBS (shrunken LFC)"),
                 ylim = c(-4, 4))
  dev.off()

  # Volcano plot
  volcano_df <- res_df %>%
    mutate(
      sig = case_when(
        dynamics == "gainAcc" ~ "gainAcc",
        dynamics == "loseAcc" ~ "loseAcc",
        TRUE                  ~ "ns"
      )
    )
  pdf(file.path(OUT_DIR, "DCA_results", paste0(ct, "_volcano.pdf")),
      width = 7, height = 6)
  print(
    ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
      geom_point(size = 0.6, alpha = 0.5) +
      scale_color_manual(values = c(gainAcc = "firebrick", loseAcc = "steelblue", ns = "grey70")) +
      geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), lty = 2, color = "black", linewidth = 0.4) +
      geom_hline(yintercept = -log10(PADJ_CUTOFF),        lty = 2, color = "black", linewidth = 0.4) +
      labs(title  = paste0(ct, ": IFNG vs PBS"),
           x      = "log2 Fold Change (shrunken)",
           y      = "-log10(padj)",
           color  = "Dynamics") +
      theme_bw(base_size = 12)
  )
  dev.off()

  # Store
  dca_results[[ct]]  <- res_df
  dds_objects[[ct]]  <- dds
  dynamics_list[[ct]] <- res_df %>% select(peak_name, dynamics)

  # Save individual DESeq2 object
  saveRDS(dds, file = file.path(OUT_DIR, "DESeq2_objects",
                                paste0(ct, ".dds.rds")))
  cat(sprintf("  Saved: %s.dds.rds\n\n", ct))
}

## ============================================================
##  7. SUMMARY TABLE
##     Rows = peaks, columns = coordinates + annotation +
##     status_{CT}_{COND} + dynamics_{CT}_IFNG + DCA stats
## ============================================================

cat("Building peak summary table...\n")

summary_tbl <- peak_coords  # chr, start, end, peak_name

# Add peak size
summary_tbl <- summary_tbl %>%
  mutate(peak_size_bp = end - start)

# Add annotation
summary_tbl <- summary_tbl %>%
  left_join(anno_df, by = "peak_name")

# Add open/closed status columns
for (ct in CELLTYPES) {
  for (cond in CONDITIONS) {
    key    <- paste(ct, cond, sep = "_")
    col_nm <- paste0("status_", key)
    if (!is.null(status_list[[key]])) {
      summary_tbl[[col_nm]] <- status_list[[key]][summary_tbl$peak_name]
    } else {
      summary_tbl[[col_nm]] <- NA_character_
    }
  }
}

# Add dynamics columns (IFNG vs PBS per celltype)
for (ct in CELLTYPES) {
  col_nm <- paste0("dynamics_", ct, "_IFNG")
  if (!is.null(dynamics_list[[ct]])) {
    summary_tbl <- summary_tbl %>%
      left_join(dynamics_list[[ct]] %>% rename(!!col_nm := dynamics),
                by = "peak_name")
  } else {
    summary_tbl[[col_nm]] <- NA_character_
  }
}

# Add DESeq2 stats columns per celltype
for (ct in CELLTYPES) {
  if (!is.null(dca_results[[ct]])) {
    stats_cols <- dca_results[[ct]] %>%
      select(peak_name,
             log2FC      = log2FoldChange,
             lfcSE,
             pvalue,
             padj) %>%
      rename_with(~ paste0(ct, "_", .), -peak_name)

    summary_tbl <- summary_tbl %>%
      left_join(stats_cols, by = "peak_name")
  }
}

# Order columns logically
coord_cols    <- c("peak_name", "chr", "start", "end", "peak_size_bp")
anno_cols     <- c("feature", "SYMBOL", "GENENAME", "distanceToTSS", "annotation")
status_cols   <- paste0("status_",   outer(CELLTYPES, CONDITIONS, paste, sep="_"))
dynamics_cols <- paste0("dynamics_", CELLTYPES, "_IFNG")
stats_cols_all <- unlist(lapply(CELLTYPES, function(ct)
  paste0(ct, c("_log2FC","_lfcSE","_pvalue","_padj"))))

col_order <- c(coord_cols, anno_cols, status_cols, dynamics_cols, stats_cols_all)
col_order <- col_order[col_order %in% colnames(summary_tbl)]
summary_tbl <- summary_tbl %>% select(all_of(col_order))

# Save
out_summary <- file.path(OUT_DIR, "peak_summary_table.tsv")
write_tsv(summary_tbl, out_summary)
cat(sprintf("  Saved: %s\n", out_summary))
cat(sprintf("  Dimensions: %d peaks × %d columns\n\n",
            nrow(summary_tbl), ncol(summary_tbl)))

## ============================================================
##  8. QUICK SUMMARY PRINTOUT
## ============================================================

cat("=== Analysis Complete ===\n")
cat(format(Sys.time()), "\n\n")

cat("Output files:\n")
cat("  DESeq2 objects:  ", file.path(OUT_DIR, "DESeq2_objects/"), "\n")
cat("  DCA results:     ", file.path(OUT_DIR, "DCA_results/"), "\n")
cat("  QC plots:        ", file.path(OUT_DIR, "QC_plots/"), "\n")
cat("  Summary table:   ", out_summary, "\n\n")

cat("Open/closed threshold: CPM >", CPM_THRESHOLD,
    "in >=", MIN_DONORS_OPEN, "donors\n")
cat("DCA threshold: padj <", PADJ_CUTOFF,
    "& |log2FC| >", LFC_CUTOFF, "\n\n")

cat("Peak status summary:\n")
for (ct in CELLTYPES) {
  for (cond in CONDITIONS) {
    key <- paste(ct, cond, sep = "_")
    col <- paste0("status_", key)
    if (col %in% colnames(summary_tbl)) {
      n_open   <- sum(summary_tbl[[col]] == "open",   na.rm = TRUE)
      n_closed <- sum(summary_tbl[[col]] == "closed", na.rm = TRUE)
      cat(sprintf("  %-15s open: %6d  closed: %6d\n", key, n_open, n_closed))
    }
  }
}

cat("\nDynamics summary (IFNG vs PBS):\n")
for (ct in CELLTYPES) {
  col <- paste0("dynamics_", ct, "_IFNG")
  if (col %in% colnames(summary_tbl)) {
    cat(sprintf("  %-8s gainAcc: %5d  loseAcc: %5d  stable: %6d  filtered: %5d\n",
                ct,
                sum(summary_tbl[[col]] == "gainAcc",  na.rm = TRUE),
                sum(summary_tbl[[col]] == "loseAcc",  na.rm = TRUE),
                sum(summary_tbl[[col]] == "stable",   na.rm = TRUE),
                sum(is.na(summary_tbl[[col]]))))
  }
}
