library(tidyverse)
library(DESeq2)
library(magrittr)
library(stringr)
library(gplots)
library(ggrepel)
library(ggfortify)
library(gridExtra)
library(ggpubr)
library(ComplexHeatmap)
library(cowplot)
library(sva)
library(clusterProfiler)
library(plotly)
library(Mfuzz)
library(Biobase)

# load basic metadata table 
f = data.table::fread("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/all_samples_06102025.txt")
f$sample_long = paste0(f$path,"/",f$sample)
f.dup = f[duplicated(f$sample_short) | duplicated(f$sample_short, fromLast = TRUE), ]
f$has_rep = ifelse(f$sample_long %in% f.dup$sample_long, "yes", "no")

# append line count and read count info
f2 = data.table::fread("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/linecount.txt") %>% set_colnames(c("sample_long","linecount"))
f2$readcount = f2$linecount / 4
f = left_join(f, f2, by="sample_long") %>% dplyr::select(-sample_long)

# explore Azenta phase2 3ctk sample read counts
cs2 = data.table::fread("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/cs2_barcode.txt", header=F) %>% set_colnames(c("id","barcode"))
f.Azenta = f %>% dplyr::filter(grepl("Azenta",path))
f.Azenta$barcode = gsub("_.*","",f.Azenta$sample)
f.Azenta = left_join(f.Azenta, cs2, by="barcode")
f.Azenta$id = ordered(f.Azenta$id, levels=c("1S","2S","3S","4S","5S","9S",
                                            "10S","11S","12S","13S","17S","18S",
                                            "23S","25S","26S","28S","29S",
                                            "30S","31S","38S","39S",
                                            "40S","44S","46S"))
# readcount vs cs2 barcode
ggplot(f.Azenta) +
  geom_point(aes(x=id, y=readcount)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 5e6, color = "red", linetype = "dashed", linewidth = 1) +  # red line at 5M
  geom_hline(yintercept = 10e6, color = "green", linetype = "dashed", linewidth = 1) + # green line at 10M
  ggtitle("read count of phase2 skin-eQTL-3ctk celseq2 libs grouped by CS2 barcode\nred=5M (12/600)\ngreen=10M (36/600)")

# readcount vs pool
f.Azenta %>%
  mutate(pool = gsub(".*/","",f.Azenta$path)) %>%
  mutate(pooln = as.numeric(gsub("_.*","", gsub(".*-","", pool)))) %>%
  ggplot( . ) +
  geom_point(aes(x=pooln, y=readcount)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(1, 25, by = 1)) +
  geom_hline(yintercept = 5e6, color = "red", linetype = "dashed", linewidth = 1) +  # red line at 5M
  geom_hline(yintercept = 10e6, color = "green", linetype = "dashed", linewidth = 1) + # green line at 10M
  ggtitle("read count of phase2 skin-eQTL-3ctk celseq2 libs grouped by lib prep pool\nred=5M (12/600)\ngreen=10M (36/600)")

# write metadata table
data.table::fwrite(f, file="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/metadata_06102025.txt", quote=F, sep="\t")

# read metadata table

