#!/usr/bin/env Rscript
# written by Crystal Shan 01/2024
# manhattan plot : https://r-graph-gallery.com/101_Manhattan_plot.html

library(tidyverse)
library(ggplot2)
library(qqman)
library(ggrepel)

# load modeling result
dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/manhattan_plot/eQTL_vs_GWAS"

# load eQTL modeling result and pvalue
eqtl <- data.table::fread(paste0(dir,"/../best_associated_PBSeQTL_pairs_and_pval_and_position.txt"))
colnames(eqtl) <- c("SNP","CHR","BP","GENE","PVAL")
eqtl$CHR <- gsub("chr", "", eqtl$CHR) %>% as.numeric()
eqtl$BP <- as.numeric(eqtl$BP)

# load GWAS result
gwas <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/vitiligo/JinY2016/all_GWAS123cmh.txt")
gwas$SNP <- gsub("RS","rs",gwas$SNP)
gwas$P <- as.numeric(gwas$P)
gwas$padj.bonf <- p.adjust(gwas$P, method="bonferroni")
gwas$padj.fdr <- p.adjust(gwas$P, method="fdr")
gwas$padj.bh <- p.adjust(gwas$P, method="BH")

# only select eQTL results of common SNPs tested in GWAS
df <- inner_join(eqtl,gwas[,c("SNP","P")],by="SNP")
colnames(df) <- c("SNP","CHR","BP","GENE","MEL_PBSeQTL_PVAL","GENERAL_VITILIGO_PVAL")
df$MEL_PBSeQTL_PVAL <- as.numeric(df$MEL_PBSeQTL_PVAL)
df$GENERAL_VITILIGO_PVAL <- as.numeric(df$GENERAL_VITILIGO_PVAL)

# check for top GWAS SNPs and its eQTL pvalue
View(df %>% dplyr::filter(GENERAL_VITILIGO_PVAL < 0.00000005) %>% arrange(GENERAL_VITILIGO_PVAL))

# plot correlation plots of pvalue
p <- ggplot(df, aes(x=-log10(MEL_PBSeQTL_PVAL), y=-log10(GENERAL_VITILIGO_PVAL))) +
  geom_point(size=0.1, alpha=0.1) +
  ggtitle("all common tested SNPs")

png(paste0(dir,"/correlation_all_common_SNPs.png"))
p
dev.off()

# plot correlation plots of significant eQTLs
this.df <- df %>% dplyr::filter(MEL_PBSeQTL_PVAL < 0.00000001)
this.df$CHR <- as.factor(this.df$CHR)
p <- ggplot(this.df, aes(x=-log10(MEL_PBSeQTL_PVAL), y=-log10(GENERAL_VITILIGO_PVAL), color=CHR)) +
  geom_point(size=0.1) +
  ggtitle("MEL_PBSeQTL_PVAL < 0.00000001")
png(paste0(dir,"/correlation_sig_MELPBSeQTL_1E-8_SNPs.png"))
p
dev.off()

# plot correlation plots of significant eQTLs
this.df <- df %>% dplyr::filter(GENERAL_VITILIGO_PVAL < 0.00000005)
this.df$CHR <- as.factor(this.df$CHR)
p <- ggplot(this.df, aes(x=-log10(MEL_PBSeQTL_PVAL), y=-log10(GENERAL_VITILIGO_PVAL), color=CHR)) +
  geom_point(size=0.1) +
  ggtitle("GENERAL_VITILIGO_PVAL < 0.00000005")
png(paste0(dir,"/correlation_sig_VITILIGO_GWAS_SNPs_5E-8.png"))
p
dev.off()


# correlation between sig GWAS SNPs' GWAS Pval and eQTL Pval
this.df <- df %>% dplyr::filter(GENERAL_VITILIGO_PVAL < 0.00000005)
cor(this.df$MEL_PBSeQTL_PVAL, this.df$GENERAL_VITILIGO_PVAL, method="spearman")
