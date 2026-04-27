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
  library(ggplot2)
  library(qpdf)
  library(tidyverse)
  library(cowplot)
  library(ggtext)
})

# ---- make heatmaps and bar charts ----
#====================================================
# 1. Input data
#====================================================
df <- tribble(
  ~celltype, ~cytokine, ~DE_direction, ~n_DEGs, ~n_DEG_cyt_eQTL, ~n_baseline_shared, ~n_stim_specific, ~ratio,
  "KRT", "IFNB", "up",  646,  31,   5,  26, 0.04798762,
  "KRT", "IFNG", "up", 1469, 144,  32, 112, 0.09802587,
  "KRT", "TNF",  "up",  706,  57,   8,  49, 0.08073654,
  "MEL", "IFNB", "up", 1357,  36,  16,  20, 0.02652911,
  "MEL", "IFNG", "up", 3641, 251,  39, 212, 0.06893711,
  "MEL", "TNF",  "up", 1350, 108,  15,  93, 0.08000000,
  "FRB", "IFNB", "up", 1632,  68,  34,  34, 0.04166667,
  "FRB", "IFNG", "up",  862,  78,  35,  43, 0.09048724,
  "FRB", "TNF",  "up", 1648, 172,  36, 136, 0.10436893
) %>%
  mutate(
    celltype = factor(celltype, levels = c("KRT", "MEL", "FRB")),
    cytokine = factor(cytokine, levels = c("IFNB", "IFNG", "TNF"))
  )
df <- df %>%
  mutate(
    cytokine = recode(
      cytokine,
      IFNB = "IFN\u03B2",
      IFNG = "IFN\u03B3",
      TNF  = "TNF\u03B1"
    )
  )

# Optional: longer display names for slide plots
celltype_labels <- c(
  KRT = "Keratinocyte",
  MEL = "Melanocyte",
  FRB = "Fibroblast"
)

# Common theme
theme_slide <- function(base_size = 18) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold"),
      axis.title = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
}

#====================================================
# 2. Plot 1: Heatmap of #DEGs
#====================================================
p_deg <- ggplot(df, aes(x = cytokine, y = celltype, fill = n_DEGs)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = comma(n_DEGs)), size = 6, fontface = "bold") +
  scale_y_discrete(labels = celltype_labels) +
  scale_fill_gradient(
    low = "#deebf7",
    high = "#08519c",
    labels = comma
  ) +
  labs(
    title = "Cytokine-upregulated genes",
    fill = "# DEGs"
  ) +
  theme_slide(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

#====================================================
# 3. Plot 2: Heatmap of % DEGs with cytokine eQTL
#====================================================
df_pct <- df %>%
  mutate(
    pct_label = percent(ratio, accuracy = 0.1),
    fraction_label = paste0(n_DEG_cyt_eQTL, "/", comma(n_DEGs))
  )

p_pct <- ggplot(df_pct, aes(x = cytokine, y = celltype, fill = ratio)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(
    aes(label = paste0(pct_label, "\n(", fraction_label, ")")),
    size = 5.2,
    fontface = "bold",
    lineheight = 0.95
  ) +
  scale_y_discrete(labels = celltype_labels) +
  scale_fill_gradient(
    low = "#efedf5",
    high = "#6a51a3",
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Genetic contribution to cytokine-responsive genes",
    fill = "% of upregulated DEGs\nwith cytokine eQTL"
  ) +
  theme_slide(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

#====================================================
# 4. Plot 3: Stacked bar plot
#====================================================

plot_df <- df %>%
  mutate(
    celltype = factor(celltype, levels = c("KRT", "MEL", "FRB")),
    cytokine = factor(cytokine, levels = c("IFN\u03B2", "IFN\u03B3", "TNF\u03B1"))
  ) %>%
  pivot_longer(
    cols = c(n_baseline_shared, n_stim_specific),
    names_to = "gene_class",
    values_to = "n_genes"
  ) %>%
  mutate(
    gene_class = factor(
      gene_class,
      levels = c("n_baseline_shared", "n_stim_specific"),
      labels = c("baseline-shared", "stimulation-specific")
    )
  )

facet_label_df <- df %>%
  mutate(
    celltype = factor(celltype, levels = c("KRT", "MEL", "FRB")),
    cytokine = factor(cytokine, levels = c("IFN\u03B2", "IFN\u03B3", "TNF\u03B1")),
    
    ratio_label = scales::percent(ratio, accuracy = 0.1),
    
    facet_label = paste0(
      cytokine,
      "<br><span style='font-size:8pt;'>",
      n_DEG_cyt_eQTL, " of ", format(n_DEGs, big.mark=","), " DEGs (", ratio_label, ")",
      "</span>"
    )
    
  ) %>%
  select(celltype, cytokine, facet_label)


plot_df <- plot_df %>%
  left_join(facet_label_df, by = c("celltype", "cytokine"))

# named vectors for labellers within each celltype panel
facet_labels_list <- split(facet_label_df, facet_label_df$celltype)
facet_labeller_list <- lapply(facet_labels_list, function(x) {
  labs <- x$facet_label
  names(labs) <- as.character(x$cytokine)
  labs
})

ymax <- max(plot_df$n_genes)
y_upper <- ymax * 1.13

make_celltype_panel <- function(ct) {
  
  sub <- plot_df %>%
    filter(celltype == ct)
  
  this_labeller <- as_labeller(facet_labeller_list[[as.character(ct)]])
  
  ggplot(sub, aes(x = gene_class, y = n_genes, fill = gene_class)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.35) +
    geom_text(aes(label = n_genes), vjust = -0.35, size = 3) +
    facet_wrap(
      ~ cytokine,
      nrow = 1,
      labeller = this_labeller
    ) +
    scale_fill_manual(values = c(
      "baseline-shared" = "#088199",
      "stimulation-specific" = "#AA4A44"
    )) +
    scale_y_continuous(
      limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = ggtext::element_markdown(face = "bold", size = 11, lineheight = 1.1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

p_krt <- make_celltype_panel("KRT") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
p_mel <- make_celltype_panel("MEL") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
p_frb <- make_celltype_panel("FRB") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

p_legend <- make_celltype_panel("KRT") +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

legend <- cowplot::get_legend(p_legend)

label_krt <- ggdraw() +
  draw_label("Keratinocyte", fontface = "bold", size = 10, angle = 90)

label_mel <- ggdraw() +
  draw_label("Melanocyte", fontface = "bold", size = 10, angle = 90)

label_frb <- ggdraw() +
  draw_label("Fibroblast", fontface = "bold", size = 10, angle = 90)

row_krt <- plot_grid(
  label_krt, p_krt,
  ncol = 2,
  rel_widths = c(0.05, 1),
  align = "h"
)

row_mel <- plot_grid(
  label_mel, p_mel,
  ncol = 2,
  rel_widths = c(0.05, 1),
  align = "h"
)

row_frb <- plot_grid(
  label_frb, p_frb,
  ncol = 2,
  rel_widths = c(0.05, 1),
  align = "h"
)

panels <- plot_grid(
  row_krt,
  row_mel,
  row_frb,
  ncol = 1,
  rel_heights = c(1, 1, 1),
  align = "v"
)

title <- ggdraw() +
  draw_label(
    "Cytokine eQTL genes among DEGs",
    fontface = "bold",
    size = 10
  )

final_plot <- plot_grid(
  title,
  legend,
  panels,
  ncol = 1,
  rel_heights = c(0.08, 0.07, 1)
)

print(final_plot)

#====================================================
# Save all plots
#====================================================
ggsave("~/Downloads/nl/human/skin/eQTLs/DREG/2025models/data/heatmap_nDEGs_induced.pdf", p_deg, width = 7.3, height = 5.2, device = cairo_pdf)
ggsave("~/Downloads/nl/human/skin/eQTLs/DREG/2025models/data/heatmap_cytokine_eQTL_genes_among_DEGs_induced.pdf", p_pct, width = 9, height = 5.2, device = cairo_pdf)
ggsave("~/Downloads/nl/human/skin/eQTLs/DREG/2025models/data/barplot_cytokine_eQTL_genes_classes.pdf", final_plot, width = 6.3, height = 6.5, device = cairo_pdf)

