# written by shuo.shan@umassmed.edu 03/2024
library(tidyverse)
library(magrittr)

# load reQTL and TF CHIP intersect result
dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/"
f = data.table::fread(paste0(dir,"/MELreQTL1E-05_overlapping_ChIPpeaks.bed"))
colnames(f) <- c("chr","start","end","SNP","REF","ALT","TF","ENCODE_ID")
data = f[,c("SNP","TF")] %>% distinct()

data2 = data %>% 
write(data, file="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/classifier/MELreQTL1E-05")
# Calculate Jaccard similarity scores for all TF pairs
snp_list = unique(data$SNP)
tf_list = unique(data$TF)

idx_pairs = as.data.frame( t(combn(tf_list, 2)) )
idx_pairs$jaccard = rep(0,nrow(idx_pairs))
idx_pairs$intersect_overlapping_SNP = rep(0,nrow(idx_pairs))
idx_pairs$union_overlapping_SNP = rep(0,nrow(idx_pairs))
for (i in 1:nrow(idx_pairs)) {
  TF1 = idx_pairs[i,]$V1
  TF2 = idx_pairs[i,]$V2
  SNP_TF1 = data[which(data$TF==TF1),]$SNP
  SNP_TF2 = data[which(data$TF==TF2),]$SNP
  
  # intersect
  this_intersect = length( intersect(SNP_TF1, SNP_TF2) )
  
  # union
  this_union = length( union(SNP_TF1, SNP_TF2) )
  
  # Jaccard SImilary Index
  jaccard_index <- this_intersect / this_union
  idx_pairs$jaccard[i] <- jaccard_index
  idx_pairs$intersect_overlapping_SNP[i] <- this_intersect
  idx_pairs$union_overlapping_SNP[i] <- this_union
}

  
overlap_freq <- table(idx_pairs$V1) %>% as.data.frame()
colnames(overlap_freq) <- c("TF","Freq")  





dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/"
f = data.table::fread(paste0(dir,"/MELreQTL1E-05_overlapping_TFpeaks_with_modelingresults.txt"))
  