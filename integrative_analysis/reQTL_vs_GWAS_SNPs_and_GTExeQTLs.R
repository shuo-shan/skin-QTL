library(shiny)
library(DT)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
# gsea and go
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggupset)
#### load data ####
load("~/Downloads/nl/human/skin/eQTLs/website/data/skineQTL_website_plotting_data.RData")

bigtable.frb <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_FRB_reQTL.txt",
                      quote=F, sep="\t")
bigtable.krt <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_KRT_reQTL.txt",
                                  quote=F, sep="\t")
bigtable.mel <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_MEL_reQTL.txt",
                                  quote=F, sep="\t")


#### calculate the number of reQTL genes ####
gene.krt <- bigtable.krt %>% 
  dplyr::filter(p.betaComp10KPermut < 0.001) %>%
  pull(gene) %>%
  unique()

gene.mel <- bigtable.mel %>% 
  dplyr::filter(p.betaComp10KPermut < 0.001) %>%
  pull(gene) %>%
  unique()

gene.frb <- bigtable.frb %>% 
  dplyr::filter(p.betaComp10KPermut < 0.001) %>%
  pull(gene) %>%
  unique()

unique.frb <- setdiff(gene.frb, union(gene.krt, gene.mel))
unique.mel <- setdiff(gene.mel, union(gene.krt, gene.frb))
unique.krt <- setdiff(gene.krt, union(gene.mel, gene.frb))
genes <- union(gene.krt, union(gene.mel, gene.frb))

#### calculate the number of reQTLs ####
reqtl.krt <- bigtable.krt %>% 
  dplyr::filter(p.betaComp10KPermut < 0.001) %>%
  pull(QTL) %>%
  unique()

reqtl.mel <- bigtable.mel %>% 
  dplyr::filter(p.betaComp10KPermut < 0.001) %>%
  pull(QTL) %>%
  unique()

reqtl.frb <- bigtable.frb %>% 
  dplyr::filter(p.betaComp10KPermut < 0.001) %>%
  pull(QTL) %>%
  unique()

reqtl <- union(reqtl.krt, union(reqtl.mel, reqtl.frb))

#### overlap with GWAS SNP mapped genes ####
gwas.table <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/GWAS-catalog-all-associations-v1.0.tsv")
gwas.table.filt <- gwas.table %>%
  dplyr::filter(SNPS %in% reqtl)

df <- gwas.table.filt %>% dplyr::select(`DISEASE/TRAIT`,SNPS, `REPORTED GENE(S)`, MAPPED_GENE) %>%
  tidyr::separate_rows(MAPPED_GENE,sep=" - ") %>%
  tidyr::separate_rows(`REPORTED GENE(S)`, sep=", ") %>%
  left_join( . , bigtable.krt[,c("QTL","gene")], by=c("SNPS"="QTL")) %>%
  na.omit()

df2 <- df[which(df$gene==df$`REPORTED GENE(S)` | df$gene==df$MAPPED_GENE),]

# filter for skin and autoimmune disease related traits
df3 <- df2 %>% dplyr::filter(`DISEASE/TRAIT` %in% unique(df2$`DISEASE/TRAIT`)[c(53,7,47,1,2,3,5,6,7,8,11,15,51,52,18)])


#### compare eQTL result to GTEx skin eQTL result ####
gtex <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/literature/GTEx/GTEx_Analysis_v8_eQTL/Skin_Not_Sun_Exposed_Suprapubic.v8.egenes.txt")
gtex.filt <- gtex[,c("rs_id_dbSNP151_GRCh38p7","ref","alt","slope","pval_nominal","gene_name")]

eqtl.krt <- bigtable.krt %>% dplyr::select(REF,ALT,QTL,gene,PBSeQTL_pval,PBSeQTL_beta)
eqtl.mel <- bigtable.mel %>% dplyr::select(REF,ALT,QTL,gene,PBSeQTL_pval,PBSeQTL_beta)
eqtl.frb <- bigtable.frb %>% dplyr::select(REF,ALT,QTL,gene,PBSeQTL_pval,PBSeQTL_beta)
eQTL_vs_gtex.krt <- left_join(eqtl.krt, gtex.filt, by=c("QTL"="rs_id_dbSNP151_GRCh38p7")) %>% na.omit()
eQTL_vs_gtex.mel <- left_join(eqtl.mel, gtex.filt, by=c("QTL"="rs_id_dbSNP151_GRCh38p7")) %>% na.omit()
eQTL_vs_gtex.frb <- left_join(eqtl.frb, gtex.filt, by=c("QTL"="rs_id_dbSNP151_GRCh38p7")) %>% na.omit()
eQTL_vs_gtex <- rbind(eQTL_vs_gtex.krt,eQTL_vs_gtex.mel,eQTL_vs_gtex.frb) %>% distinct() %>%
  mutate(samegene=(gene==gene_name)) %>%
  dplyr::filter(samegene==TRUE)

eQTL_vs_gtex[which(sign(eQTL_vs_gtex$PBSeQTL_beta)!=sign(eQTL_vs_gtex$slope)),] %>% View()

# out of the 30 snp:gene pair that also have data in GTEx data, 29 pairs have the same direction. 
# 1 pair rs7825980:LYNX1 have the opposite direction, but our data shows very convincing results in all three cell types.