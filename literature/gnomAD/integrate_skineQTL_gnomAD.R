#!/usr/bin/env Rscript
# written by Crystal Shan 01/2024
# explore the loss-of-function (LoF) intolerance in genes linked to an eQTL or reQTL

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggplot2)

dir <- "~/Downloads/nl/human/skin/eQTLs"
res.fname <- paste0(dir,"/website/data/modeling_results/modeling_results_featureSelectedModel_phenotypeRankNormCPM.txt")
res.tab <- data.table::fread(res.fname)
res.tab$PBSeQTL_pval <- as.numeric(res.tab$PBSeQTL_pval)
res.tab$IFNeQTL_pval <- as.numeric(res.tab$IFNeQTL_pval)
res.tab$reQTL_pval <- as.numeric(res.tab$reQTL_pval)
gnomad.tab <- data.table::fread(paste0(dir,"/literature/gnomAD/gnomad.v4.0.constraint_metrics.tsv"))

# remove result entries whose SNPs are not in HWE

# filter gnomad table by QTL-linked genes
genes.PBSeQTL <- res.tab %>% dplyr::filter(PBSeQTL_pval<0.000001 & !is.na(PBSeQTL_pval)) %>% pull(GENE) %>% unique()
genes.IFNeQTL <- res.tab %>% dplyr::filter(IFNeQTL_pval<0.000001 & !is.na(IFNeQTL_pval)) %>% pull(GENE) %>% unique()
genes.reQTL <- res.tab %>% dplyr::filter(reQTL_pval<0.000001 & !is.na(reQTL_pval)) %>% pull(GENE) %>% unique()
genes.QTL <- unique(c(genes.PBSeQTL, genes.IFNeQTL, genes.reQTL))

# get the pLI and LOEUF scores of QTL genes. 
# For genes with more than one trancript, pick the transcript with the highest pLI score.
gnomad.QTLgenes <- gnomad.tab %>% 
  dplyr::filter(gene %in% genes.reQTL) %>%
  group_by(gene) %>%
  arrange(desc(lof.pLI)) %>%
  slice(1)
gnomad.QTLgenes.highpLI <- gnomad.QTLgenes %>% dplyr::filter(lof.pLI >= 0.9)
hist(gnomad.QTLgenes$lof.pLI, breaks=100,
     main="LoF pLI of all QTL genes")
abline(v=0.95, col="red", lwd=2)

# distribution of pLI in vitiligo GWAS susceptibility genes
vit.genes <- data.table::fread(paste0(dir,"/literature/vitiligo_susceptibility_genes.txt"))
gnomad.vitgenes <- gnomad.tab %>% 
  dplyr::filter(gene %in% vit.genes$gene) %>%
  group_by(gene) %>%
  arrange(desc(lof.pLI)) %>%
  slice(1)
hist(gnomad.vitgenes$lof.pLI, breaks=100,
     main="LoF pLI of all vitiligo susceptibility genes by GWAS")
abline(v=0.95, col="red", lwd=2)
