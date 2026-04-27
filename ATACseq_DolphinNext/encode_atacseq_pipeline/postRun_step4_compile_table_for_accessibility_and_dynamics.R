library(tidyverse)
library(DESeq2)
library(future.apply)

plan(multisession, workers = 3)
log <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")

# ─────────────────────────────────────────
# 1. Load data
# ─────────────────────────────────────────
dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks"
raw <- read.table(paste0(dir, "/master_peaks_multicov.txt"), header = TRUE, sep = "\t", check.names = FALSE)

peak_info <- raw[, 1:4]
colnames(peak_info) <- c("chr", "start", "end", "peak_name")

count_mat <- raw[, 5:ncol(raw)]
rownames(count_mat) <- peak_info$peak_name
count_mat_int <- round(as.matrix(count_mat))
storage.mode(count_mat_int) <- "integer"

log("loaded data")

# ─────────────────────────────────────────
# 2. Build metadata
# ─────────────────────────────────────────
meta <- data.frame(sample = colnames(count_mat_int)) %>%
  mutate(
    donor     = str_extract(sample, "F\\d+"),
    celltype  = str_extract(sample, "(?<=_)(FRB|KRT|MEL)(?=_)"),
    condition = case_when(
      str_detect(sample, "IFNG|IFN") ~ "IFNG",
      str_detect(sample, "PBS")      ~ "PBS",
      TRUE ~ NA_character_
    ),
    group = paste(celltype, condition, sep = "_")
  ) %>%
  column_to_rownames("sample")

stopifnot(all(rownames(meta) == colnames(count_mat_int)))
log("Built metadata")
log("Samples:", nrow(meta))
log("Groups:", paste(unique(meta$group), collapse = ", "))

# ─────────────────────────────────────────
# 3. Open status from MACS2 narrowPeak overlap
# ─────────────────────────────────────────
log("Loading narrowPeak overlap table...")

MIN_DONORS <- length(unique(meta$donor)) %/% 2  # 3

overlap_raw <- read.table(
  paste0(dir, "/peak_overlap/master_peaks_overlap.txt"),
  header = TRUE, sep = "\t", check.names = FALSE
)

# Sanity check: overlap table samples match meta
overlap_samples <- colnames(overlap_raw)[-1]  # drop peak_name column
stopifnot(all(overlap_samples %in% rownames(meta)))

# Build open status per celltype x condition
GROUPS <- unique(meta$group)

open_merged <- purrr::map_dfc(GROUPS, function(grp) {
  grp_samples <- rownames(meta)[meta$group == grp]
  donor_counts <- rowSums(overlap_raw[, grp_samples, drop = FALSE])
  tibble(!!paste0("open_", grp) := donor_counts >= MIN_DONORS)
}) %>%
  mutate(peak_name = overlap_raw$peak_name) %>%
  relocate(peak_name)

log("Open status done. Peaks open in at least one group:",
    sum(rowSums(open_merged[, -1]) > 0))

# ─────────────────────────────────────────
# 4. Per-celltype DESeq2
# ─────────────────────────────────────────
run_celltype <- function(celltype_str) {
  
  log("=============================")
  log("Processing celltype:", celltype_str)
  
  ct_samples <- rownames(meta)[meta$celltype == celltype_str]
  ct_meta    <- meta[ct_samples, ]
  ct_counts  <- count_mat_int[, ct_samples]
  
  ct_meta$condition <- factor(ct_meta$condition, levels = c("PBS", "IFNG"))
  ct_meta$donor     <- factor(ct_meta$donor)
  
  dds <- DESeqDataSetFromMatrix(
    countData = ct_counts,
    colData   = ct_meta,
    design    = ~ donor + condition
  )
  
  keep <- rowSums(counts(dds) >= 5) >= 3
  dds  <- dds[keep, ]
  log("Peaks passing filter in", celltype_str, ":", sum(keep))
  
  # VST for QC plot only
  vst_ct  <- vst(dds, blind = FALSE)
  norm_ct <- assay(vst_ct)
  
  pdf(paste0(dir, "/VST_distribution_", celltype_str, ".pdf"))
  hist(norm_ct, breaks = 100,
       main = paste("VST distribution -", celltype_str),
       xlab = "VST value",
       col  = "steelblue")
  dev.off()
  
  # DESeq2 DE
  dds <- DESeq(dds)
  
  res <- results(dds,
                 contrast = c("condition", "IFNG", "PBS"),
                 alpha    = 0.05) %>%
    as.data.frame() %>%
    rownames_to_column("peak_name") %>%
    mutate(
      celltype  = celltype_str,
      direction = case_when(
        padj < 0.05 & log2FoldChange >  1 ~ "gain_accessibility",
        padj < 0.05 & log2FoldChange < -1 ~ "lose_accessibility",
        TRUE                               ~ "no_change"
      )
    )
  
  log("Summary for", celltype_str, ":")
  log(capture.output(dplyr::count(res, direction)) %>% paste(collapse = "\n"))
  
  list(
    norm_mat   = norm_ct,
    de_results = res
  )
}

celltypes <- c("FRB", "KRT", "MEL")
results   <- future_lapply(celltypes, run_celltype) %>%
  setNames(celltypes)

# ─────────────────────────────────────────
# 5. Merge DESeq2 results across celltypes
# ─────────────────────────────────────────
de_merged <- purrr::map(results, "de_results") %>%
  bind_rows()

de_wide <- de_merged %>%
  select(peak_name, celltype, log2FoldChange, padj, direction) %>%
  pivot_wider(
    names_from  = celltype,
    values_from = c(log2FoldChange, padj, direction),
    names_glue  = "{celltype}_{.value}"
  )

# ─────────────────────────────────────────
# 6. Build final annotation table
# ─────────────────────────────────────────
annotation_table <- peak_info %>%
  left_join(open_merged, by = "peak_name") %>%
  left_join(de_wide,     by = "peak_name") %>%
  # NAs in open columns = peak filtered in all celltypes → safely FALSE
  dplyr::mutate(across(starts_with("open_"), ~replace_na(., FALSE)))

log("Final annotation table:", nrow(annotation_table), "peaks x", ncol(annotation_table), "columns")

# ─────────────────────────────────────────
# 7. Save
# ─────────────────────────────────────────
write_tsv(annotation_table, paste0(dir, "/ATACseq_peak_annotation.tsv"))
write_tsv(de_merged,        paste0(dir, "/ATACseq_DESeq2_full_results.tsv"))

purrr::iwalk(results, function(res, ct) {
  as.data.frame(res$norm_mat) %>%
    rownames_to_column("peak_name") %>%
    write_tsv(paste0(dir, "/ATACseq_VST_normalized_", ct, ".tsv"))
})

log("Done!")