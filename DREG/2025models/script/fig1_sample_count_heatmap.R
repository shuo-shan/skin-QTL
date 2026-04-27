library(ggplot2)
library(dplyr)

dir="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models"


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
ct       <- "KRT"
this_gene <- "IRF3"
cytokine  <- "TNF"

# Paths ----
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)

chunk_id_lookup <- data.table::fread(paste0(dir, "/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[chunk_id_lookup$gene == this_gene, ]$chunk)
chunk_id <- sprintf("%03d", chunk_id)

geno_file <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
meta_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
PEER_FILE <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_", ct, "_PBS-IFNG-IFNB-TNF.tsv")

meta_all <- readr::read_tsv(meta_file, show_col_types = FALSE)
peers    <- data.table::fread(PEER_FILE)
colnames(peers)[1] <- "sample"

nPEER    <- 10
all_peer_cols <- setdiff(colnames(peers), "sample")
peer_use <- all_peer_cols[5:(5 + nPEER - 1)]

get_n_donors <- function(condition_name) {
  meta_cond <- meta_all %>%
    dplyr::filter(celltype == ct, condition == condition_name) %>%
    dplyr::select(sample, donor)
  
  peer_cond <- as.data.frame(peers) %>%
    dplyr::filter(sample %in% meta_cond$sample) %>%
    dplyr::select(sample, all_of(peer_use)) %>%
    dplyr::left_join(meta_cond, by = "sample")
  
  nrow(peer_cond)
}

N_PBS  <- get_n_donors("PBS")
N_cyto <- get_n_donors(cytokine)

cat(sprintf("N for SuSiE — PBS: %d | %s: %d\n", N_PBS, cytokine, N_cyto))


# plot ----
df <- data.frame(
  CellType = c("Melanocyte","Melanocyte","Melanocyte","Melanocyte",
               "Keratinocyte","Keratinocyte","Keratinocyte","Keratinocyte",
               "Fibroblast","Fibroblast","Fibroblast","Fibroblast"),
  Cytokine = c("PBS","IFNG","IFNB","TNF",
               "PBS","IFNG","IFNB","TNF",
               "PBS","IFNG","IFNB","TNF"),
  value = c(49,48,31,30,
            71,73,42,42,
            87,86,55,55)
)

# value = c(51,51,33,32,
#           72,75,43,43,
#           90,90,58,58)
# reorder axes to match your requirement
df$CellType <- factor(df$CellType, levels = c("Fibroblast", "Keratinocyte", "Melanocyte"))
df$Cytokine <- factor(df$Cytokine, levels = c("PBS", "IFNB", "IFNG", "TNF"))

pdf(file="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/fig1_sample_count_heatmap.pdf",
    width=4, height=3)
cytokine_labels <- c(
  "IFNG"     = expression(IFN*gamma),
  "IFNB"     = expression(IFN*beta),
  "TNFalpha" = expression(TNF*alpha)
)
p <- ggplot(df, aes(x = CellType, y = Cytokine, fill = value)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(
    label = value,
    color = value > (min(df$value) + (max(df$value) - min(df$value)) * 0.5)
  ), size = 5) +
  scale_fill_gradientn(
    colours = c("#EEF2F7", "#B8C4DA", "#6F8FB8", "#2E5F95"),
    limits = c(min(df$value), max(df$value))
  ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  scale_y_discrete(limits = rev(levels(df$Cytokine)), labels = cytokine_labels) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )
print(p)
dev.off()