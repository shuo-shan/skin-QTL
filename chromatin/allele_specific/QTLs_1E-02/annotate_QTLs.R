#!/usr/bin/env Rscript
# written by Crystal Shan 04/2024

library(tidyverse)
library(magrittr)

dir="~/Downloads/nl/human/skin/eQTLs/chromatin/allele_specific/QTLs_1E-02"
qtl = data.table::fread(paste0(dir,"/QTL.bed"),header=F,sep="\t") 
colnames(qtl) = c("chr","start","end","ID","REF","ALT")

flist = list.files( dir , pattern="QTL_1E-2*")
mark = flist %>%
  gsub("QTL_1E-2_overlapping_AlleleSpecific_", "", .) %>%
  gsub("HM-ChIP-seq_", "", .) %>%
  gsub("TF-ChIP-seq_", "", .) %>%
  gsub("hetSNVs_", "", .) %>%
  gsub(".bed", "", .) %>%
  gsub("_pooledtissue", "", . ) %>%
  gsub("_skin", "", .) %>%
  gsub("ATAC-seq", "ATACseq", . ) %>%
  gsub("RNA-seq", "RNAseq", .)

result_df = data.frame(ID = character(), AlleleSpecificMark = character())
for (i in 1:length(flist)) {
  f = data.table::fread(paste0(dir,"/",flist[i]), header=TRUE, sep="\t")
  this_mark = mark[i]
  this_df = data.frame(ID = unique(f$ID), AlleleSpecificMark = this_mark)
  result_df = rbind(result_df, this_df)
}
result_df2 <- result_df %>% 
  group_by(ID) %>% 
  summarize(AlleleSpecificMark=paste(AlleleSpecificMark,collapse=','))
count <- result_df %>% group_by(ID) %>% mutate(AlleleSpecificMarkCount = n()) %>% dplyr::select(-AlleleSpecificMark) %>% distinct()
result_df3 <- left_join(count, result_df2)

result_final <- left_join(qtl, result_df3, by="ID")
result_final$AlleleSpecificMarkCount[is.na(result_final$AlleleSpecificMarkCount)] = 0

data.table::fwrite(result_final, file=paste0(dir,"/QTL_1E-02_annotated_AlleleSpecificMark_summary.bed"), 
                   quote=FALSE, sep="\t")
