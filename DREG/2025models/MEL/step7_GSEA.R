#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(magrittr)
  library(tibble)
  library(fgsea)
  library(msigdbr)
  library(forcats)
  library(ggrepel)
})


# ---------------------- Set-up --------------------- ####
ct <- "MEL"
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
keyword <- "MEL_TNF_reQTL"
set.seed(1)

META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
CPM_FILE="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"

# get QTL modeling statistics
stats  <- fread(paste0(dir, "/eigenMT/results/", keyword,".eigenMT.txt"), header = TRUE)
stats$qval <- p.adjust(stats$p_gene_eigenMT, method="BH")
stats <- stats[, c("gene","lead_snp","pmin","p_gene_eigenMT","qval")]

# calculate ranks: named numeric vector, decreasing ####
ranks <- setNames(-log10(pmax(stats$p_gene_eigenMT, 1e-300)), stats$gene)
ranks <- sort(ranks, decreasing = TRUE)

# GSEA on KEGG pathways ####
msig_kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>%
  select(gs_name, gene_symbol)

pathways <- split(msig_kegg$gene_symbol, msig_kegg$gs_name)

res_kegg <- fgseaMultilevel(pathways = pathways, stats = ranks,
  minSize = 15, maxSize = 500, nPermSimple = 10000) %>% arrange(pval) %>%
  as.data.frame() %>%
  mutate(collection = "KEGG")

# GSEA on Hallmark pathways ####
msig_h <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(gs_name, gene_symbol)

pathways <- split(msig_h$gene_symbol, msig_h$gs_name)

res_hallmark <- fgseaMultilevel(pathways = pathways, stats = ranks,
  minSize = 15, maxSize = 500, nPermSimple = 10000) %>% arrange(pval) %>%
  as.data.frame() %>%
  mutate(collection = "Hallmark")

# GSEA on GO BP pathways ####
msig_bp <- msigdbr(
  species = "Homo sapiens",
  collection = "C5",
  subcollection = "BP"
)
pathways <- split(msig_bp$gene_symbol, msig_bp$gs_name)

res_GOBP <- fgseaMultilevel(pathways = pathways, stats = ranks,
  minSize = 15, maxSize = 500, nPermSimple = 10000) %>% arrange(pval) %>%
  as.data.frame() %>%
  mutate(collection = "GOBP")

res_all <- bind_rows(res_kegg, res_hallmark, res_GOBP)
nrow(res_all[which(res_all$padj<0.1),])
View(head(res_all %>% arrange(padj), n=30))

# Summarize “best evidence” per collection ####
best_tbl <- res_all %>%
  group_by(collection) %>%
  summarise(
    n_sets = n(),
    min_padj = min(padj, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    top_pathway = pathway[which.min(padj)],
    top_NES = NES[which.min(padj)],
    .groups = "drop"
  ) %>%
  mutate(neglog10_minpadj = -log10(pmax(min_padj, 1e-300)))

View(best_tbl)

# Plot: qqplot ####
make_qq <- function(df, label){
  df <- df %>% filter(!is.na(pval)) %>% arrange(pval)
  df$exp <- -log10(ppoints(nrow(df)))
  df$obs <- -log10(pmax(df$pval, 1e-300))
  df$collection <- label
  df
}

qq_df <- bind_rows(
  make_qq(res_kegg, "KEGG"),
  make_qq(res_hallmark, "Hallmark"),
  make_qq(res_GOBP, "GO BP")
)

pqq <- ggplot(qq_df, aes(x = exp, y = obs)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~collection, scales = "free") +
  theme_classic(base_size = 12) +
  labs(title="GSEA pathway p-values are near-null at baseline (MEL PBS eQTLs)",
       x="Expected -log10(p)", y="Observed -log10(p)")








# View top 20 of each GSEA collection
View(head(res_GOBP[, c("pathway","pval","NES","size","leadingEdge")], n = 20))
# Dot plot (NES vs –log10 FDR) ####
plot_df <- res_all %>%
  mutate(
    log10FDR = -log10(pmax(padj, 1e-300)),
    sig = padj < 0.1
  ) %>%
  arrange(desc(log10FDR)) %>%
  slice_head(n = 20) %>%                # top 20 pathways
  mutate(
    pathway = fct_reorder(pathway, NES) # reorder for plotting
  )

p <- ggplot(plot_df, aes(x = NES, y = log10FDR)) +
  geom_point(aes(size = size, color = sig), alpha = 0.85) +
  geom_text_repel(
    aes(label = pathway),
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.2,
    segment.color = "grey60",
    segment.size = 0.2
  ) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D55E00")) +
  theme_classic(base_size = 13) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = expression(-log[10]("FDR")),
    size = "Gene set size",
    color = "FDR < 0.1",
    title = "IFNB reQTLs"
  )

ggsave(paste0(dir,"/plots/IFNB_reQTL_GSEA_dotplot.png"), p, width = 8, height = 5, dpi = 300)

# Focused GSEA pathway plot for 1 pathway at a time ####
top_row <- res_all %>%
  arrange(padj, pval) %>%
  slice(1)
top_pathway <- top_row$pathway
top_pathway

ranks <- setNames(-log10(stats$pmin), stats$gene)
ranks <- sort(ranks, decreasing = TRUE)

p_enrich <- plotEnrichment(
  pathways[[top_pathway]],
  ranks
) +
  labs(
    title = paste0("TNF reQTL GSEA: ", top_pathway),
    subtitle = paste0(
      "NES=", round(top_row$NES, 3),
      " | FDR=", signif(top_row$padj, 3),
      " | p=", signif(top_row$pval, 3)
    ),
    x = "Ranked genes",
    y = "Running enrichment score"
  ) +
  theme_classic(base_size = 13)

print(p_enrich)

# In cluster/RStudio-server, save to file (most reliable)
png(paste0(dir,"/plots/TNF_reQTL_topPathway_enrichment.png"), width = 3200, height = 2000, res = 300)
print(p_enrich)
dev.off()
