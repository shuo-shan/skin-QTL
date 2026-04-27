#!/usr/bin/env Rscript
# =============================================================================
# step12_explore_colocsusie_results.R
#
# PURPOSE:
#   Aggregate coloc-susie results across genes and cytokines for one cell type,
#   merge with DE gene lists, classify genes into H1/H2/H3/H4 categories,
#   and produce summary plots.
#
# USAGE:
#   Interactive — set ct and stim at the top and run section by section.
#   Adjust thresholds and filters as needed during exploration.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(ggtext)
  library(scales)
})

# =============================================================================
# 0. PARAMETERS — change these to explore different cell types / cytokines
# =============================================================================
ct   <- "FRB"    # cell type: FRB = fibroblast, MEL = melanocyte
stim <- "IFNG"   # cytokine: IFNG, IFNB, TNF

# Coloc posterior thresholds for classification
# PP.H4 > h4_thresh  -> H4 (shared signal)
# PP.H3 > h3_thresh  -> H3 (distinct signal)
# These are the standard thresholds used in the coloc literature.
# Lower threshold (e.g. 0.5) is more permissive, 0.8 is more stringent.
h4_thresh <- 0.5
h3_thresh <- 0.5

# =============================================================================
# 1. PATHS
# =============================================================================
base_dir    <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)
coloc_dir   <- paste0(base_dir, "/coloc_susie")
deg_file    <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData"
metagene_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt"
out_dir     <- paste0(base_dir, "/coloc_susie/summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. UTILITY: META-GENE NAME CONVERSION
# (carried over from my old analysis pipeline)
# =============================================================================
metagene_dict <- data.table::fread(metagene_file) %>%
  dplyr::select(-id) %>%
  dplyr::distinct()

convert_meta_genes <- function(gene_vector) {
  this_metagenes_idx <- which(grepl("meta", gene_vector))
  temp1 <- data.frame(
    meta = gene_vector[-this_metagenes_idx],
    name = gene_vector[-this_metagenes_idx]
  )
  temp2 <- data.frame(meta = gene_vector[this_metagenes_idx]) %>%
    dplyr::left_join(metagene_dict, by = "meta")
  rbind(temp1, temp2) %>%
    magrittr::set_colnames(c("meta", "gene"))
}

# =============================================================================
# 3. LOAD AND AGGREGATE COLOC-SUSIE RESULTS
# =============================================================================
# Each gene x cytokine combination produces one _coloc_summary.tsv file.
# Aggregate all of them for this cell type across all cytokines.

cytokines <- c("IFNB", "IFNG", "TNF")

all_results_raw <- fread(paste0(coloc_dir,"/coloc_susie_summary_",stim,".txt"))

cat(sprintf("\nLoaded %d gene x cytokine rows total.\n", nrow(all_results_raw)))
cat(sprintf("Genes: %d unique | Cytokines: %s\n",
            dplyr::n_distinct(all_results_raw$gene),
            paste(unique(all_results_raw$cytokine), collapse = ", ")))

# =============================================================================
# 4. CLASSIFY EACH GENE x CYTOKINE INTO H0/H1/H2/H3/H4
# =============================================================================
# Classification logic:
#   skipped     — job errored out before coloc ran (check skip_reason)
#   no_signal   — neither condition has a significant eQTL (H0 territory)
#   H1          — PBS signal only (constitutive eQTL lost upon stimulation)
#   H2          — cytokine signal only (stimulation-induced eQTL, e.g. IRF3)
#   H3_susie    — coloc.susie PP.H3 > threshold (distinct causal variants)
#   H4_susie    — coloc.susie PP.H4 > threshold (shared causal variant)
#   H3_abf      — coloc.abf PP.H3 > threshold (SuSiE failed, abf fallback)
#   H4_abf      — coloc.abf PP.H4 > threshold (SuSiE failed, abf fallback)
#   ambiguous   — signal in both conditions but no hypothesis dominates
#
# NOTE on H1/H2 via n_cs:
#   n_cs_PBS == 0 does not strictly mean no PBS signal — SuSiE may have found
#   no credible set even with a marginal signal. But combined with coloc.abf
#   PP.H2 being dominant, it's a reliable indicator. The coloc_call is a
#   classification for exploration, not a formal statistical test.

all_results <- all_results_raw %>%
  dplyr::mutate(
    coloc_call = dplyr::case_when(
      # Job failed before coloc ran
      !is.na(skip_reason) & str_starts(skip_reason, "ERROR")
      ~ "error",
      # No signal in either condition — H0
      !is.na(skip_reason) & str_starts(skip_reason, "no_signal_either")
      ~ "H0_no_signal",
      # Signal in cytokine only — stimulation-induced eQTL (e.g. IRF3)
      # Use n_cs_PBS == 0 as primary indicator; abf PP.H2 as confirmation
      n_cs_PBS == 0 & n_cs_cyto > 0
      ~ "H2_stimulation_induced",
      # Signal in PBS only — constitutive eQTL that disappears upon stimulation
      n_cs_PBS > 0 & n_cs_cyto == 0
      ~ "H1_constitutive_lost",
      # coloc.susie results (preferred — handles multiple signals)
      !is.na(susie_PP_H3) & susie_PP_H3 > h3_thresh
      ~ "H3_distinct_susie",
      !is.na(susie_PP_H4) & susie_PP_H4 > h4_thresh
      ~ "H4_shared_susie",
      # coloc.abf fallback (SuSiE failed or no credible sets)
      !is.na(abf_PP_H3) & abf_PP_H3 > h3_thresh
      ~ "H3_distinct_abf",
      !is.na(abf_PP_H4) & abf_PP_H4 > h4_thresh
      ~ "H4_shared_abf",
      # Both conditions have signal but no clear winner
      TRUE ~ "ambiguous"
    ),
    
    # Simplified H-category for plotting (collapse susie/abf distinction)
    H_category = dplyr::case_when(
      coloc_call %in% c("H0_no_signal", "error") ~ "H0/skipped",
      coloc_call == "H1_constitutive_lost"        ~ "H1\nconstitutive\nlost",
      coloc_call == "H2_stimulation_induced"       ~ "H2\nstimulation\ninduced",
      coloc_call %in% c("H3_distinct_susie",
                        "H3_distinct_abf")         ~ "H3\ndistinct\nsignal",
      coloc_call %in% c("H4_shared_susie",
                        "H4_shared_abf")           ~ "H4\nshared\nsignal",
      TRUE                                         ~ "ambiguous"
    ),
    
    # Flag: was result from coloc.susie or coloc.abf fallback?
    method_used = dplyr::case_when(
      str_detect(coloc_call, "susie") ~ "coloc.susie",
      str_detect(coloc_call, "abf")   ~ "coloc.abf (fallback)",
      coloc_call %in% c("H1_constitutive_lost",
                        "H2_stimulation_induced")  ~ "SuSiE CS count",
      TRUE                                         ~ NA_character_
    )
  )

cat("\nClassification summary (all cytokines combined):\n")
print(table(all_results$coloc_call))
cat("\nBy cytokine:\n")
print(table(all_results$cytokine, all_results$H_category))

# =============================================================================
# 5. LOAD DE GENES AND MERGE
# =============================================================================
# DE genes come from a pre-saved RData file containing one data.frame per
# cell type x cytokine combination, named deg.{ct}.{stim} (lowercase).
# Each has columns: gene (or meta), log2FoldChange, padj, etc.

load(deg_file)

# Helper to get DE genes for a specific ct + stim
get_deg <- function(ct_lower, stim_lower, direction = c("both", "up", "down")) {
  direction <- match.arg(direction)
  var <- ls(pattern = paste0("^deg\\.", ct_lower, "\\.", stim_lower),
            envir = .GlobalEnv)
  if (length(var) == 0) {
    cat(sprintf("  No DE object found for deg.%s.%s\n", ct_lower, stim_lower))
    return(NULL)
  }
  deg_raw <- get(var)
  colnames(deg_raw)[1] <- "meta"
  
  # Convert meta-gene IDs to gene names
  meta_dict <- convert_meta_genes(deg_raw$meta)
  deg_converted <- dplyr::left_join(meta_dict, deg_raw, by = "meta")
  
  # Filter by direction
  if (direction == "up")   deg_converted <- dplyr::filter(deg_converted, log2FoldChange > 0)
  if (direction == "down")  deg_converted <- dplyr::filter(deg_converted, log2FoldChange < 0)
  
  dplyr::select(deg_converted, gene, log2FoldChange, padj) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
}

# Load DE genes for all cytokines and add is_DE / DE direction flags
ct_lower <- tolower(ct)

deg_flags <- purrr::map_dfr(tolower(cytokines), function(s) {
  deg <- get_deg(ct_lower, s)
  if (is.null(deg)) return(NULL)
  tibble::tibble(
    gene        = deg$gene,
    cytokine    = toupper(s),
    is_DE       = TRUE,
    DE_direction = dplyr::case_when(
      deg$log2FoldChange > 0 ~ "up",
      deg$log2FoldChange < 0 ~ "down",
      TRUE                   ~ "unchanged"
    ),
    log2FC      = deg$log2FoldChange,
    DE_padj     = deg$padj
  )
}) %>%
  dplyr::distinct(gene, cytokine, .keep_all = TRUE)

# Merge DE flags into coloc results
all_results <- all_results %>%
  dplyr::left_join(deg_flags, by = c("gene", "cytokine")) %>%
  dplyr::mutate(
    is_DE        = tidyr::replace_na(is_DE, FALSE),
    DE_direction = tidyr::replace_na(DE_direction, "not_DE")
  )

cat(sprintf("\nDE gene overlap:\n"))
cat(sprintf("  Total gene x cytokine rows: %d\n", nrow(all_results)))
cat(sprintf("  Rows with is_DE=TRUE: %d\n", sum(all_results$is_DE)))
cat(sprintf("  DE genes with H1: %d\n",
            sum(all_results$is_DE & str_detect(all_results$coloc_call, "H1"))))
cat(sprintf("  DE genes with H2: %d\n",
            sum(all_results$is_DE & str_detect(all_results$coloc_call, "H2"))))
cat(sprintf("  DE genes with H3: %d\n",
            sum(all_results$is_DE & str_detect(all_results$coloc_call, "H3"))))
cat(sprintf("  DE genes with H4: %d\n",
            sum(all_results$is_DE & str_detect(all_results$coloc_call, "H4"))))

# =============================================================================
# 6. SUBSET TO SPECIFIC ct + stim FOR INTERACTIVE EXPLORATION
# =============================================================================
# This is the interactive section — change ct and stim at the top to explore.

this_deg    <- get_deg(ct_lower, tolower(stim), direction = "both")
this_deg_up <- get_deg(ct_lower, tolower(stim), direction = "up")

# All coloc results for this ct x stim
this_coloc <- all_results %>%
  dplyr::filter(celltype == ct, cytokine == stim)

cat(sprintf("\n--- %s | %s ---\n", ct, stim))
cat(sprintf("  Total genes tested: %d\n", nrow(this_coloc)))
cat(sprintf("  DE genes: %d\n", sum(this_coloc$is_DE)))
cat(sprintf("  DE genes up: %d\n", nrow(this_deg_up)))
cat("\n  Coloc call distribution:\n")
print(table(this_coloc$coloc_call))

# DE genes subset — all directions
this_coloc_DE <- this_coloc %>%
  dplyr::filter(is_DE)

cat(sprintf("\n  Coloc call among DE genes:\n"))
print(table(this_coloc_DE$coloc_call))

# DE genes, up-regulated only
this_coloc_DE_up <- this_coloc %>%
  dplyr::filter(gene %in% this_deg_up$gene)

cat(sprintf("\n  Coloc call among up-regulated DE genes:\n"))
print(table(this_coloc_DE_up$coloc_call))

# Quick look at H3 DE genes — these are the most interesting
cat(sprintf("\n  H3 DE genes (distinct regulatory architecture + DE):\n"))
this_coloc_DE %>%
  dplyr::filter(str_detect(coloc_call, "H3")) %>%
  dplyr::select(gene, coloc_call, susie_PP_H3, susie_PP_H4,
                abf_PP_H3, abf_PP_H4, DE_direction, log2FC) %>%
  dplyr::arrange(dplyr::desc(dplyr::coalesce(susie_PP_H3, abf_PP_H3))) %>%
  print(n = 20)

# =============================================================================
# 7. PLOTTING
# =============================================================================
# Color palette for H categories — consistent across all plots
h_colors <- c(
  "H0/skipped"          = "grey85",
  "H1\nconstitutive\nlost"    = "#4393C3",
  "H2\nstimulation\ninduced"  = "#92C5DE",
  "H3\ndistinct\nsignal"      = "#D6604D",
  "H4\nshared\nsignal"        = "#4DAF4A",
  "ambiguous"                 = "grey60"
)

# 7a. OVERALL H-CATEGORY DISTRIBUTION — bar chart, all cytokines stacked
plot_overall_distribution <- function(results, title_suffix = "") {
  plot_df <- results %>%
    dplyr::count(cytokine, H_category) %>%
    dplyr::group_by(cytokine) %>%
    dplyr::mutate(pct = n / sum(n))
  
  ggplot(plot_df, aes(x = cytokine, y = pct, fill = H_category)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = h_colors, name = "Coloc\ncategory") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title    = sprintf("Coloc category distribution%s", title_suffix),
      subtitle = sprintf("%s | thresholds: H3/H4 PP > %.1f", ct, h3_thresh),
      x        = "Cytokine stimulation",
      y        = "Fraction of genes tested"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.title    = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
}

# 7b. H-CATEGORY DISTRIBUTION IN DE GENES ONLY
plot_DE_distribution <- function(results, direction = "both", title_suffix = "") {
  sub <- results %>%
    dplyr::filter(is_DE)
  if (direction == "up")   sub <- dplyr::filter(sub, DE_direction == "up")
  if (direction == "down")  sub <- dplyr::filter(sub, DE_direction == "down")
  
  dir_label <- switch(direction,
                      "up"   = " (up-regulated DE genes)",
                      "down" = " (down-regulated DE genes)",
                      " (all DE genes)"
  )
  
  plot_df <- sub %>%
    dplyr::count(cytokine, H_category) %>%
    dplyr::group_by(cytokine) %>%
    dplyr::mutate(pct = n / sum(n))
  
  ggplot(plot_df, aes(x = cytokine, y = pct, fill = H_category)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = h_colors, name = "Coloc\ncategory") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title    = sprintf("Coloc category distribution%s%s", dir_label, title_suffix),
      subtitle = sprintf("%s | thresholds: H3/H4 PP > %.1f", ct, h3_thresh),
      x        = "Cytokine stimulation",
      y        = "Fraction of DE genes"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.title    = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
}

# 7c. POSTERIOR PROBABILITY SCATTER — PP.H3 vs PP.H4 per gene
# Useful for seeing the full distribution, not just threshold-classified genes
plot_posterior_scatter <- function(results, use_susie = TRUE,
                                   highlight_DE = TRUE, stim_filter = NULL) {
  sub <- results
  if (!is.null(stim_filter)) sub <- dplyr::filter(sub, cytokine == stim_filter)
  
  if (use_susie) {
    sub <- sub %>%
      dplyr::mutate(
        PP_H3 = dplyr::coalesce(susie_PP_H3, abf_PP_H3),
        PP_H4 = dplyr::coalesce(susie_PP_H4, abf_PP_H4),
        method = ifelse(!is.na(susie_PP_H3), "coloc.susie", "coloc.abf")
      )
  } else {
    sub <- sub %>%
      dplyr::mutate(PP_H3 = abf_PP_H3, PP_H4 = abf_PP_H4, method = "coloc.abf")
  }
  
  sub <- sub %>% dplyr::filter(!is.na(PP_H3) & !is.na(PP_H4))
  
  p <- ggplot(sub, aes(x = PP_H4, y = PP_H3)) +
    geom_point(aes(color = is_DE, shape = method), alpha = 0.6, size = 1.8) +
    geom_hline(yintercept = h3_thresh, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = h4_thresh, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.9, y = 0.05, label = "H4\n(shared)",
             color = "#4DAF4A", fontface = "bold", size = 3.5) +
    annotate("text", x = 0.05, y = 0.9, label = "H3\n(distinct)",
             color = "#D6604D", fontface = "bold", size = 3.5) +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"),
                       labels = c("TRUE" = "DE gene", "FALSE" = "not DE"),
                       name = NULL) +
    scale_shape_manual(values = c("coloc.susie" = 16, "coloc.abf" = 17),
                       name = "Method") +
    facet_wrap(~cytokine) +
    labs(
      title    = sprintf("Posterior probabilities: H3 vs H4 | %s", ct),
      subtitle = "Dashed lines = classification threshold",
      x        = "PP.H4 (shared causal variant)",
      y        = "PP.H3 (distinct causal variants)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  p
}

# 7d. COUNT SUMMARY TABLE PLOT — clean table showing N per category
plot_count_table <- function(results, subset_label = "All genes") {
  tbl <- results %>%
    dplyr::count(cytokine, H_category) %>%
    tidyr::pivot_wider(names_from = H_category, values_from = n,
                       values_fill = 0)
  
  # Convert to long for plotting as a tile table
  tbl_long <- tbl %>%
    tidyr::pivot_longer(-cytokine, names_to = "H_category", values_to = "n") %>%
    dplyr::mutate(
      H_category = factor(H_category, levels = names(h_colors))
    )
  
  ggplot(tbl_long, aes(x = H_category, y = cytokine, fill = H_category)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = n), fontface = "bold", size = 4) +
    scale_fill_manual(values = h_colors, guide = "none") +
    labs(
      title    = sprintf("Gene counts per coloc category | %s | %s", ct, subset_label),
      x        = NULL, y        = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid  = element_blank()
    )
}

# =============================================================================
# 8. RENDER AND SAVE PDF
# =============================================================================
# One PDF with all summary plots.
# Adjust width/height per page as needed.

out_pdf <- file.path(out_dir, sprintf("coloc_susie_summary_%s.pdf", ct))

pdf(out_pdf, width = 12, height = 8)

# Page 1: Overall distribution, all genes
print(plot_overall_distribution(all_results))

# Page 2: Distribution in DE genes (all directions)
print(plot_DE_distribution(all_results, direction = "both"))

# Page 3: Distribution in up-regulated DE genes
print(plot_DE_distribution(all_results, direction = "up",
                           title_suffix = ""))

# Page 4: Distribution in down-regulated DE genes
print(plot_DE_distribution(all_results, direction = "down",
                           title_suffix = ""))

# Page 5: Posterior scatter, all genes
print(plot_posterior_scatter(all_results))

# Page 6: Posterior scatter, DE genes only
print(plot_posterior_scatter(
  dplyr::filter(all_results, is_DE),
  highlight_DE = FALSE
) + labs(subtitle = "DE genes only | dashed lines = classification threshold"))

# Page 7: Count table, all genes
print(plot_count_table(all_results, subset_label = "All genes"))

# Page 8: Count table, DE genes only
print(plot_count_table(
  dplyr::filter(all_results, is_DE),
  subset_label = "DE genes only"
))

# Page 9: Per-cytokine detail for the selected stim (interactive focus)
print(
  plot_posterior_scatter(
    dplyr::filter(all_results, cytokine == stim),
    stim_filter = stim
  ) + labs(title = sprintf("Posterior probabilities | %s | %s only", ct, stim))
)

dev.off()
cat(sprintf("\nPDF written to: %s\n", out_pdf))

# =============================================================================
# 9. QUICK ACCESS: GENE LISTS OF INTEREST
# =============================================================================
# Run these interactively to pull specific gene lists for follow-up

# H3 genes that are also DE (your primary story)
H3_DE_genes <- all_results %>%
  dplyr::filter(str_detect(coloc_call, "H3"), is_DE) %>%
  dplyr::select(gene, cytokine, coloc_call, susie_PP_H3, abf_PP_H3,
                DE_direction, log2FC, flag_susie_not_converged) %>%
  dplyr::arrange(cytokine, dplyr::desc(dplyr::coalesce(susie_PP_H3, abf_PP_H3)))

cat("\nH3 DE genes across all cytokines:\n")
print(H3_DE_genes, n = 30)

# H4 genes that are also DE (shared regulatory architecture + expression change)
H4_DE_genes <- all_results %>%
  dplyr::filter(str_detect(coloc_call, "H4"), is_DE) %>%
  dplyr::select(gene, cytokine, coloc_call, susie_PP_H4, abf_PP_H4,
                DE_direction, log2FC, flag_susie_not_converged) %>%
  dplyr::arrange(cytokine, dplyr::desc(dplyr::coalesce(susie_PP_H4, abf_PP_H4)))

cat("\nH4 DE genes across all cytokines:\n")
print(H4_DE_genes, n = 30)

# H2 genes (stimulation-induced eQTLs) — like IRF3
H2_genes <- all_results %>%
  dplyr::filter(coloc_call == "H2_stimulation_induced") %>%
  dplyr::select(gene, cytokine, n_cs_cyto, abf_PP_H2, is_DE, DE_direction) %>%
  dplyr::arrange(cytokine, dplyr::desc(abf_PP_H2))

cat("\nH2 (stimulation-induced eQTL) genes:\n")
print(H2_genes, n = 20)

# Genes with ambiguous results — worth manual inspection
ambiguous_genes <- all_results %>%
  dplyr::filter(coloc_call == "ambiguous", is_DE) %>%
  dplyr::select(gene, cytokine, susie_PP_H3, susie_PP_H4,
                abf_PP_H3, abf_PP_H4, flag_susie_not_converged) %>%
  dplyr::arrange(cytokine)

cat("\nAmbiguous DE genes (worth manual inspection):\n")
print(ambiguous_genes, n = 20)