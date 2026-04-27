##  QTL Study — Summary Figures
##  Figures:
##    1. Grouped bar chart: n_genes by condition, faceted by QTLtype
##    2. 
##    4. Heatmap: condition × celltype grid, faceted by QTLtype
##    5. reQTL composition plot: PBS-only / shared / reQTL-only
##
##  Required packages:
##    install.packages(c("tidyverse", "patchwork", "scales", "ggtext", "cowplot"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(UpSetR)
  library(ComplexUpset)
  library(ggvenn)
  library(ggtext)
  library(cowplot)
})

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models")
outdir <- paste0(dir,"/data")

# =============== Global Settings ======================
# ── Global aesthetics ────────────────────────────────────────────────────────
ct_colors <- c(fibroblast = "#CCCCFF", keratinocyte = "#DAF7A6", melanocyte = "#964B00")
ct_labels  <- c(fibroblast = "Fibroblast", keratinocyte = "Keratinocyte", melanocyte = "Melanocyte")
cond_order  <- c("PBS", "IFNB", "IFNG", "TNF")
cond_labels <- c(PBS = "PBS\n(ctrl)", IFNB = "IFNβ", IFNG = "IFNγ", TNF = "TNFα")

prop_colors <- c(
  "PBS eQTL only"           = "#4393C3",
  "Shared PBS eQTL & reQTL" = "#762A83",
  "reQTL only"              = "#D6604D"
)

base_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, color = "grey45", lineheight = 1.3),
    axis.title         = element_text(face = "bold"),
    strip.background   = element_rect(fill = "grey93", colour = NA),
    strip.text         = element_text(face = "bold", size = 11),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.text        = element_text(size = 10)
  )

# ── Load data ────────────────────────────────────────────────────────────────
# summarize gene counts
dat <- fread(paste0(dir,"/data/all_QTL_genes_FDR05.txt")) 
dat$celltype <- dplyr::recode(
  dat$celltype,
  MEL = "melanocyte",
  KRT = "keratinocyte",
  FRB = "fibroblast"
)
dat <- dat %>%
  mutate(
  celltype  = factor(celltype,  levels = c("fibroblast", "keratinocyte", "melanocyte")),
  condition = factor(condition, levels = c("PBS","IFNB","IFNG","TNF")),
  QTLtype   = factor(QTLtype,   levels = c("eQTL", "reQTL"))
)
  
celltypes  <- unique(dat$celltype)    # FRB, KRT, MEL
conditions <- unique(dat$condition)   # PBS, IFNB, IFNG, TNF
qtltypes   <- unique(dat$QTLtype)     # eQTL, reQTL

# Complete summary (zeros for missing combos, e.g. PBS reQTL, KRT IFNB reQTL)
counts <- dat %>%
  group_by(celltype, condition, QTLtype) %>%
  summarise(n_genes = n_distinct(gene), .groups = "drop") %>%
  complete(celltype, condition, QTLtype, fill = list(n_genes = 0))

## ============================================================
##  FIGURE 1. — BAR CHART
##  Conditions on x-axis, n_genes on y, colored by cell type,
##  faceted by QTLtype (eQTL / reQTL), free y-scales.
## ============================================================

p_eqtl <- counts %>%
  dplyr::filter(QTLtype=="eQTL") %>%
  droplevels() %>%  
  ggplot(., aes(x = condition, y = n_genes, fill = celltype)) +
  geom_col(
    position = position_dodge(width = 0.80),
    width = 0.74, color = "white", linewidth = 0.25
  ) +
  geom_text(
    aes(label = ifelse(n_genes > 0, comma(n_genes), "")),
    position = position_dodge(width = 0.80),
    vjust = -0.45, size = 2.7, fontface = "bold", color = "grey30"
  ) +
  scale_fill_manual(values = ct_colors, labels = ct_labels, name = NULL) +
  scale_x_discrete(labels = cond_labels) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    labels = comma
  ) +
  labs(
    title    = "QTL gene discoveries across conditions and cell types",
    subtitle = "Each bar = genes with ≥1 significant QTL (eigenMT per-gene FWER, then BH FDR < 0.05 across genes)",
    x        = "Stimulation condition",
    y        = "Number of genes with significant QTL"
  ) +
  base_theme +
  theme(legend.position = "right")

p_reqtl <- counts %>%
  dplyr::filter(QTLtype=="reQTL") %>%
  dplyr::filter(!(condition=="PBS" & QTLtype=="reQTL")) %>%
  droplevels() %>%  
  ggplot(., aes(x = condition, y = n_genes, fill = celltype)) +
  geom_col(
    position = position_dodge(width = 0.80),
    width = 0.74, color = "white", linewidth = 0.25
  ) +
  geom_text(
    aes(label = ifelse(n_genes > 0, comma(n_genes), "")),
    position = position_dodge(width = 0.80),
    vjust = -0.45, size = 2.7, fontface = "bold", color = "grey30"
  ) +
  scale_fill_manual(values = ct_colors, labels = ct_labels, name = NULL) +
  scale_x_discrete(labels = cond_labels) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    labels = comma
  ) +
  labs(
    x        = "Stimulation condition",
    y        = "Number of genes with significant QTL"
  ) +
  base_theme +
  theme(legend.position = "right")

fig1 <- p_eqtl + p_reqtl + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = paste0(outdir, "/Fig1_QTL_barplot.pdf"),
  plot = fig1,
  width = 10.5, height = 5.5,
  device = cairo_pdf
)
message("Figure 1 saved.")


## ============================================================
##  FIGURE 2 — Heatmap
##  Rows = conditions, columns = cell types,
##  separate colour scales for eQTL (blue) and reQTL (red).
## ============================================================

make_hm_panel <- function(qtl_val, pal_low, pal_high) {
  d <- counts %>%
    filter(QTLtype == qtl_val) %>%
    mutate(
      label  = ifelse(n_genes == 0, "—", comma(n_genes)),
      n_fill = ifelse(n_genes == 0, NA_real_, n_genes)
    )
  
  max_val <- max(d$n_genes)
  
  ggplot(d, aes(x = celltype, y = fct_rev(condition), fill = n_fill)) +
    geom_tile(color = "white", linewidth = 1.4) +
    geom_text(
      aes(
        label = label,
        color = ifelse(n_genes > max_val * 0.55, "light", "dark")
      ),
      size = 5, fontface = "bold"
    ) +
    scale_color_manual(
      values = c(light = "white", dark = "grey15"),
      guide  = "none"
    ) +
    scale_x_discrete(labels = ct_labels, position = "top") +
    scale_y_discrete(labels = rev(cond_labels)) +
    scale_fill_gradient(
      low      = pal_low,
      high     = pal_high,
      na.value = "grey88",
      name     = "Genes",
      labels   = comma,
      limits   = c(1, max_val)
    ) +
    labs(title = qtl_val, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 14, hjust = 0.5,
                                     margin = margin(b = 8)),
      axis.text.x     = element_text(face = "bold", size = 11),
      axis.text.y     = element_text(size = 11),
      panel.grid      = element_blank(),
      legend.position = "right"
    )
}

hm_eqtl  <- make_hm_panel("eQTL",  "#DEEBF7", "#08306B")
hm_reqtl <- make_hm_panel("reQTL", "#FEE5D9", "#99000D")

fig2 <- (hm_eqtl | hm_reqtl) +
  plot_annotation(
    title    = "QTL gene discovery heatmap",
    subtitle = "Number of genes with ≥1 significant QTL per condition × cell type  (— = none detected)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey45")
    )
  )

ggsave(paste0(outdir,"/Fig2_QTL_heatmap.pdf"), fig2, width = 10.5, height = 5.2, device = cairo_pdf)
message("Figure 2 saved.")


## ============================================================
##  FIGURE 3 — reQTL composition plot
##
##  For each celltype × cytokine:
##    Universe = PBS eQTL genes  ∪  cytokine reQTL genes
##    Three exclusive categories:
##      (A) PBS eQTL only      — constitutive eQTL, not a reQTL
##      (B) Shared             — significant under PBS eQTL AND cytokine reQTL
##      (C) reQTL only         — novel stimulus-induced QTL, absent in PBS eQTL
##
##  Two-panel design:
##    Panel A  Full stack (log10) to show scale context
##    Panel B  Zoom into reQTL universe (shared + reQTL only),
##             annotated with % novel reQTL genes
## ============================================================

cytokines <- c("IFNB", "IFNG", "TNF")

prop_df <- map_dfr(levels(dat$celltype), function(ct) {
  pbs_genes <- dat %>%
    filter(celltype == ct, condition == "PBS", QTLtype == "eQTL") %>%
    pull(gene) %>% unique()
  
  map_dfr(cytokines, function(cond) {
    rq_genes <- dat %>%
      filter(celltype == ct, condition == cond, QTLtype == "reQTL") %>%
      pull(gene) %>% unique()
    
    n_shared   <- length(intersect(pbs_genes, rq_genes))
    n_pbs_only <- length(setdiff(pbs_genes,  rq_genes))
    n_rq_only  <- length(setdiff(rq_genes,   pbs_genes))
    n_rq_total <- length(rq_genes)
    
    tibble(
      celltype                  = ct,
      condition                 = cond,
      `PBS eQTL only`           = n_pbs_only,
      `Shared PBS eQTL & reQTL` = n_shared,
      `reQTL only`              = n_rq_only,
      n_rq_total                = n_rq_total,
      pct_novel = ifelse(n_rq_total > 0,
                         round(100 * n_rq_only / n_rq_total, 1),
                         NA_real_)
    )
  })
}) %>%
  mutate(
    celltype  = factor(celltype,  levels = c("fibroblast", "keratinocyte", "melanocyte")),
    condition = factor(condition, levels = cytokines)
  )

# Long format for ggplot stacking
prop_long <- prop_df %>%
  pivot_longer(
    cols      = c("PBS eQTL only", "Shared PBS eQTL & reQTL", "reQTL only"),
    names_to  = "category",
    values_to = "n_genes"
  ) %>%
  mutate(category = factor(category,
                           levels = c("PBS eQTL only", "Shared PBS eQTL & reQTL", "reQTL only")))

# ── Panel A: full stack on log10 scale ────────────────────────
panelA <- ggplot(prop_long,
                 aes(x = condition, y = n_genes, fill = category)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  facet_wrap(~ celltype, ncol = 3,
             labeller = labeller(celltype = ct_labels)) +
  scale_fill_manual(values = prop_colors, name = NULL) +
  scale_x_discrete(
    labels = c(IFNB = "IFNβ", IFNG = "IFNγ", TNF = "TNFα")
  ) +
  scale_y_log10(
    breaks = scales::log_breaks(n = 5),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.08))
  ) +
  annotation_logticks(sides = "l", size = 0.3, color = "grey55") +
  labs(
    title = "A   Full scale (log₁₀)  —  PBS eQTL background + reQTL genes",
    x     = NULL,
    y     = "Number of genes (log₁₀)"
  ) +
  base_theme +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.45, "cm"),
    legend.spacing.x = unit(0.3, "cm")
  )

# ── Panel B: zoom into reQTL universe ─────────────────────────
zoom_long <- prop_long %>%
  filter(category != "PBS eQTL only") %>%
  mutate(
    category = factor(
      category,
      levels = c("Shared PBS eQTL & reQTL", "reQTL only")
    )
  )

# Annotation labels: % novel, placed above top of bar
label_df <- prop_df %>%
  filter(!is.na(pct_novel), n_rq_total > 0) %>%
  mutate(
    y_top = `Shared PBS eQTL & reQTL` + `reQTL only`,
    label = paste0(pct_novel, "% novel")
  )

label_n_df <- prop_df %>%
  dplyr::filter(n_rq_total > 0) %>%
  dplyr::mutate(
    y_top = `Shared PBS eQTL & reQTL` + `reQTL only`,   # bar height in Panel B
    label_n = paste0("n=", scales::comma(n_rq_total))
  )

panelB <- ggplot(zoom_long,
                 aes(x = condition, y = n_genes, fill = category)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_text(
    data        = label_df,
    aes(x = condition, y = y_top, label = label),
    inherit.aes = FALSE,
    vjust = -0.4, size = 3.2, fontface = "italic", color = "grey25"
  ) +
  facet_wrap(~ celltype, ncol = 3,
             labeller = labeller(celltype = ct_labels)) +
  scale_fill_manual(
    values = prop_colors[c("Shared PBS eQTL & reQTL", "reQTL only")],
    name   = NULL
  ) +
  scale_x_discrete(
    labels = c(IFNB = "IFNβ", IFNG = "IFNγ", TNF = "TNFα")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.22)),
    labels = comma
  ) +
  labs(
    title = "B   Zoom  —  reQTL gene composition  (% = fraction not in PBS eQTL → 'novel')",
    x     = "Cytokine stimulation",
    y     = "Number of reQTL genes"
  ) +
  base_theme +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.45, "cm"),
    legend.spacing.x = unit(0.3, "cm")
  )

fig3 <- panelA / panelB +
  plot_annotation(
    title    = "Stimulus-specific vs constitutive genetic regulation",
    subtitle = paste0(
      "reQTL genes are classified relative to the PBS eQTL gene list within each cell type.\n",
      "'Novel' reQTL genes (red) are not significant eQTLs under PBS, ",
      "indicating stimulus-induced genetic effects not detectable at baseline."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey45", lineheight = 1.3)
    )
  )

ggsave(paste0(outdir,"/Fig3_reQTL_composition.pdf"), fig3, width = 11, height = 9.5, device = cairo_pdf)
message("Figure 3 saved.")

# ── Print composition table ────────────────────────────────────
cat("\n=== reQTL composition summary ===\n")
prop_df %>%
  select(celltype, condition,
         `PBS eQTL only`, `Shared PBS eQTL & reQTL`,
         `reQTL only`, n_rq_total, pct_novel) %>%
  arrange(celltype, condition) %>%
  print(n = Inf)

## ============================================================
## FIGURE 4. UPSET PLOTS: cytokine specificity. per celltype, genes across conditions
##    Separately for eQTL and reQTL
## ============================================================

make_upset <- function(celltype_val, qtltype_val) {
  sub <- dat %>%
    dplyr::filter(celltype == celltype_val, QTLtype == qtltype_val)
  
  if (nrow(sub) == 0) {
    message(sprintf("No data for %s %s -- skipping", celltype_val, qtltype_val))
    return(invisible(NULL))
  }
  
  # Condition order you want (fixed)
  cond_order <- c("PBS", "IFNB", "IFNG", "TNF")
  
  # Build membership wide table: one row per gene, columns=conditions (TRUE/FALSE)
  membership <- sub %>%
    dplyr::select(gene, condition) %>%
    dplyr::distinct() %>%
    dplyr::mutate(value = TRUE) %>%
    tidyr::pivot_wider(
      names_from = condition,
      values_from = value,
      values_fill = FALSE
    )
  
  # Keep only conditions that actually exist as columns (in case some missing)
  sets_present <- intersect(cond_order, colnames(membership))
  if (length(sets_present) < 2) {
    message(sprintf("Only %d condition(s) present for %s %s -- skipping upset",
                    length(sets_present), celltype_val, qtltype_val))
    return(invisible(NULL))
  }
  
  title_str <- sprintf("%s %s - Genes per Condition", celltype_val, qtltype_val)
  
  # Build plot
  p <- ComplexUpset::upset(
    membership,
    sets_present,
    name = "Condition",
    sort_sets = FALSE,
    sort_intersections = "descending",
    width_ratio = 0.2,
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(
        text = list(size = 3, vjust = -0.2)
      )
    )
  )  +
    ggplot2::ggtitle(title_str) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  print(p)
  invisible(NULL)
}

# Save all upset plots
pdf(paste0(outdir,"/Fig4_QTL_UpsetPlots_CytokineOverlap.pdf"), width = 10, height = 6)
for (ct in celltypes) {
  message("eQTL upset: ", ct)
  make_upset(ct, "eQTL")
  
  message("reQTL upset: ", ct)
  make_upset(ct, "reQTL")
}
dev.off()

message("Figure 4: Upset plots saved to Fig4_QTL_UpsetPlots_CytokineOverlap.pdf")

## ============================================================
## FIGURE 5. VENN DIAGRAMS: celltype specificity. per condition x QTLtype, genes per celltype
## ============================================================
# Color palette per celltype
make_venn <- function(cond_val, qtltype_val) {
  sub <- dat %>%
    filter(condition == cond_val, QTLtype == qtltype_val)
  
  if (nrow(sub) == 0) {
    message(sprintf("No data for %s %s -- skipping", cond_val, qtltype_val))
    return(invisible(NULL))
  }
  
  cts_present <- intersect(celltypes, unique(sub$celltype))
  gene_sets <- setNames(
    lapply(cts_present, function(ct) unique(sub$gene[sub$celltype == ct])),
    cts_present
  )
  
  p <- ggvenn(gene_sets,
              fill_color = ct_colors[cts_present],
              stroke_size = 0.5,
              set_name_size = 5,
              text_size = 4) +
    labs(title = sprintf("%s %s: Gene Overlap Across Cell Types", cond_val, qtltype_val)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

# Save all venn diagrams
pdf(paste0(outdir,"/Fig5_QTL_VennDiagrams_CelltypeOverlap.pdf"), width = 7, height = 7, device = cairo_pdf)
for (qt in qtltypes) {
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    message("Venn: ", cond, " ", qt)
    p <- make_venn(cond, qt)
    if (!is.null(p)) print(p)
  }
}
dev.off()

message("Venn diagrams saved to Fig5_QTL_VennDiagrams_CelltypeOverlap.pdf")

## ============================================================
## FIGURE 6. UPSET PLOTS: celltype specificity. per condition x QTLtype, across celltypes
## ============================================================
make_celltype_upset <- function(cond_val, qtltype_val) {
  sub <- dat %>%
    dplyr::filter(condition == cond_val, QTLtype == qtltype_val) %>%
    dplyr::select(gene, celltype) %>%
    dplyr::distinct()
  
  if (nrow(sub) == 0) return(NULL)
  
  ct_order <- c("melanocyte", "keratinocyte", "fibroblast")
  
  membership <- sub %>%
    dplyr::mutate(value = TRUE) %>%
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = value,
      values_fill = FALSE
    )
  
  sets_present <- intersect(ct_order, colnames(membership))
  if (length(sets_present) < 2) return(NULL)
  
  title_str <- sprintf("%s %s: Gene overlap across cell types", cond_val, qtltype_val)
  
  p <- ComplexUpset::upset(
    membership,
    sets_present,
    name = "Cell type",
    sort_sets = FALSE,
    sort_intersections = "descending",
    width_ratio = 0.25,
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size()
      # (set-size bars are already shown; labeling them depends on your ComplexUpset version)
    )
  ) +
    ggplot2::ggtitle(title_str) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p
}

pdf(paste0(outdir,"/Fig6_QTL_UpSetPlots_CelltypeOverlap.pdf"), width = 10, height = 6)
for (qt in qtltypes) {
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    p <- make_celltype_upset(cond, qt)
    if (!is.null(p)) print(p)
  }
}
dev.off()

