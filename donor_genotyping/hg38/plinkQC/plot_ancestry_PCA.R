library(plinkQC)
indir="/Users/crystal/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata"
name <- "data" # Because your files are test.bed, test.bim, test.fam
path2plink <- "/Users/crystal/Downloads/plink_mac_20250615/plink"
refname <- 'HapMapIII_CGRCh38'
prefixMergedDataset <- paste(name, ".", refname, sep="")

exclude_ancestry <-
  evaluate_check_ancestry(indir=indir, name=name,
                          prefixMergedDataset=prefixMergedDataset,
                          refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                               sep=""), 
                          refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                              sep=""),
                          verbose=TRUE,
                          interactive=TRUE)


library(plinkQC)
library(ggplot2)
library(data.table)
library(tidyverse)

indir <- "/Users/crystal/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata"
name  <- "data"
refname <- "HapMapIII_CGRCh38"
prefixMergedDataset <- paste(name, ".", refname, sep = "")

# Read the merged PCA eigenvec file
eigenvec <- fread(file.path(indir, paste0(prefixMergedDataset, ".eigenvec")),
                  header = FALSE)
colnames(eigenvec)[1:4] <- c("FID", "IID", "PC1", "PC2")

# Read HapMap population labels
hapmap_ids <- fread(file.path(indir, "HapMap_ID2Pop.txt"),
                    header = TRUE)  # columns: IID, Pop (or similar)

# Read HapMap population colors (optional — used to match colors)
hapmap_colors <- fread(file.path(indir, "HapMap_PopColors.txt"),
                       header = TRUE)  # columns: Pop, Color

# Merge population labels onto eigenvec
# Study samples won't match HapMap IDs, so they get NA → label as "data"
pca_df <- merge(eigenvec, hapmap_ids, by = "IID", all.x = TRUE)
pca_df$Population <- ifelse(is.na(pca_df$Pop), "data", pca_df$Pop)

# Hapmap + My data clusters ----
hapmap_pops <- unique(hapmap_colors$Pop)
pca_df$Population <- factor(pca_df$Population,
                            levels = c(hapmap_pops, "data"))

# Sort by factor level so ggplot draws in that order
pca_df <- pca_df[order(pca_df$Population), ]

# Color vector: HapMap colors + blue for study samples
pop_colors <- setNames(hapmap_colors$Color, hapmap_colors$Pop)
pop_colors["data"] <- "black"
  
pop_sizes  <- setNames(rep(1.5, length(hapmap_pops)), hapmap_pops)
pop_sizes["data"]  <- 2.5   # your samples slightly larger

pop_alphas <- setNames(rep(0.6, length(hapmap_pops)), hapmap_pops)
pop_alphas["data"] <- 1.0   # your samples fully opaque


# Hapmap clusters ----
p_pca <- pca_df %>%
  dplyr::filter(Population!="data") %>%
  ggplot(., aes(x = PC1, y = PC2,
                color = Population,
                size  = Population,
                alpha = Population)) +
  geom_point(shape = 21,             # filled circle WITH border
             aes(fill = Population),
             color = "white",        # border color for ALL dots
             stroke = 0.3,
             size = 4) +         # border thickness — keep subtle
  scale_fill_manual(values = pop_colors) +   # fill now drives color
  scale_color_manual(values = pop_colors) +  # needed to avoid legend issues
  scale_size_manual(values = pop_sizes) +
  scale_alpha_manual(values = pop_alphas) +
  labs(title = "PCA on reference genotypes",
       x = "PC1", y = "PC2") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    legend.title     = element_text(size = 10),
    plot.title       = element_text(size = 12)
  ) +
  lims(x=c(-0.05, 0.03), y=c(-0.045, 0.045) ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1)),
         size  = "none",   # suppress duplicate legends
         alpha = "none")

ggsave("~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_plot_HapMap.pdf", p_pca, width = 5, height = 4)

# Hapmap plus my donors clusters ----
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                   color = Population,
                   size  = Population,
                   alpha = Population)) +
  geom_point(shape = 21,             # filled circle WITH border
             aes(fill = Population),
             color = "white",        # border color for ALL dots
             stroke = 0.3,
             size = 4) +         # border thickness — keep subtle
  scale_fill_manual(values = pop_colors) +   # fill now drives color
  scale_color_manual(values = pop_colors) +  # needed to avoid legend issues
  scale_size_manual(values = pop_sizes) +
  scale_alpha_manual(values = pop_alphas) +
  labs(title = "PCA on combined reference and study genotypes",
       x = "PC1", y = "PC2") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    legend.title     = element_text(size = 10),
    plot.title       = element_text(size = 12)
  ) +
  lims(x=c(-0.05, 0.03), y=c(-0.045, 0.045) ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1)),
         size  = "none",   # suppress duplicate legends
         alpha = "none")
ggsave("~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/ancestry_PCA_plot_HapMapPlusMySamples.pdf", p_pca, width = 5, height = 4)
