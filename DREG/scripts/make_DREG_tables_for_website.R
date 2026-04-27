#!/usr/bin/env Rscript
# written by Crystal Shan 09/2023

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
library(LDlinkR)

# parsing input arguments
args = commandArgs(trailingOnly=TRUE)
outDir=args[1]
snp=args[2]
gene=args[3]

snp="rs4360063"
gene="ERAP2"
outDir=paste0("~/Downloads/nl/human/skin/eQTLs/website/data/plots","/",gene,"/",snp)
system(paste("mkdir", paste0("~/Downloads/nl/human/skin/eQTLs/website/data/plots","/",gene)))
system(paste("mkdir", paste0("~/Downloads/nl/human/skin/eQTLs/website/data/plots","/",gene,"/",snp)))

###### functions and loading data tables #####
rankNorm <- function (y) {
  # input y: numeric vector of CPM across all genes per donor.
  k <- 0.375 # an offset to ensure the z-score is finite. from Blom transform.
  n <- length(y)
  # Ranks.
  r <- rank(y, ties.method = "average") # if same value, same rank
  # Apply transformation.
  r.prob <- (r - k) / (n - 2 * k + 1)
  y.rankNorm <- stats::qnorm(r.prob) # because qnorm(1) is infinite and qnorm(-1) is -Inf
  return(y.rankNorm)
}

###### load-data-for-plotting ###### 
# load enhancer data
enhancerF <- "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method4/promoters_and_enhancers_surrounding_genes_all3cts.bed"
enhancerATACF <- "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ATACseq_allcts_allregions_1kbp.bed"
enhancerK27acF = "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ChIPseq_allcts_allregions_1kbp.bed"

system(paste0("awk -v g=",gene," '{if ($6==g) print $0}' ",enhancerF," > ",outDir,"/enhancerDict.txt"))
enhancerDict <- read.table(paste0(outDir,"/enhancerDict.txt"), header=F, sep="\t") %>%
  set_colnames(c("chr","start","end","name","ID","gene","celltype","regionType"))
system(paste0("rm ",outDir,"/enhancerDict.txt"))

count_ATAC <- data.table::fread(enhancerATACF) %>%
  dplyr::select(-c(chr,start,end)) %>% column_to_rownames("name")
col.sum=apply(count_ATAC,2,sum); 
CPM_atac=sweep(count_ATAC, 2, col.sum, `/`)*1000000
CPM_atac <- CPM_atac[,c("F25_KRT_PBS","F49_KRT_PBS","F55_KRT_PBS",
                        "F25_MEL_PBS","F49_MEL_PBS","F55_MEL_PBS",
                        "F25_FRB_PBS","F49_FRB_PBS","F55_FRB_PBS",
                        "F25_KRT_IFN","F49_KRT_IFN","F55_KRT_IFN",
                        "F25_MEL_IFN","F49_MEL_IFN","F55_MEL_IFN",
                        "F25_FRB_IFN","F49_FRB_IFN","F55_FRB_IFN")]

enhancerATAC <-  CPM_atac %>% rownames_to_column("name") %>% dplyr::filter(name %in% enhancerDict$ID)
rm(count_ATAC)

count_K27ac <- data.table::fread(enhancerK27acF) %>%
  dplyr::select(-c(chr,start,end)) %>% column_to_rownames("name")
col.sum=apply(count_K27ac,2,sum); 
CPM_K27ac=sweep(count_K27ac, 2, col.sum, `/`)*1000000
CPM_K27ac <- CPM_K27ac[,c("F25_KRT_PBS","F49_KRT_PBS","F55_KRT_PBS",
                          "F25_MEL_PBS","F49_MEL_PBS","F55_MEL_PBS",
                          "F25_FRB_PBS","F49_FRB_PBS","F55_FRB_PBS",
                          "F25_KRT_IFN","F49_KRT_IFN","F55_KRT_IFN",
                          "F25_MEL_IFN","F49_MEL_IFN","F55_MEL_IFN",
                          "F25_FRB_IFN","F49_FRB_IFN","F55_FRB_IFN")]

enhancerK27ac <-  CPM_K27ac %>% rownames_to_column("name") %>% dplyr::filter(name %in% enhancerDict$ID)
rm(count_K27ac)

