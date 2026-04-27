library(tidyverse)
library(magrittr)
library(ggpubr)

dir="~/Downloads/nl/human/skin/eQTLs/edQTL/output"
setwd(dir)

# load edQTL
edQTL.all = data.table::fread(input = paste0(dir,"/foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal")) %>%
  dplyr::select(,c(V1,V6,V9,V11))
colnames(edQTL.all) = c("edSite","edQTL","slope","padj")

# load GWAS table
gwas.all = data.table::fread(input = "~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/gwas_catalog_v1.0_all_associations_snps.txt")
gwas.snps = unique(gwas.all$SNP_ID_CURRENT) %>% gsub(" ","",.)

# join edQTL and GWAS table
tab.joined = left_join(edQTL.all, gwas.all, by = c("edQTL" = "SNP_ID_CURRENT"))
edQTL.in.gwas = edQTL.all %>% dplyr::filter(edQTL %in% gwas.snps)
