#!/usr/bin/env Rscript
# written by Crystal Shan 01/2024
# create a histogram for the HWE p-value distribution of European donors biallelic SNPs that are genotyped.
# genotype file processing code from original genotype merged vcf can be found: skin/eQTLs/donor_genotyping/hg38/HWE
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggplot2)

dir <- "~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/HWE"
fname <- paste0(dir,"/european_biallelic_snps_HWEpval.txt")
tab <- data.table::fread(fname)

ggplot(tab, aes(x=P_HWE)) +
  geom_histogram(bins=1000)

tab$Padj_HWE <- p.adjust(tab$P_HWE, method="BH")
ggplot(tab, aes(x=Padj_HWE)) +
  geom_histogram(bins=1000)

tab_inHWE <- tab %>% dplyr::filter(Padj_HWE > 0.05)

nrow(tab_inHWE) / nrow(tab)
