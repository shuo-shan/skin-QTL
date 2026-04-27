# written by shuo.shan@umassmed.edu 03/2024
library(tidyverse)
library(magrittr)

# load QTL and TF CHIPseq intersect result
dir="~/Downloads/nl/human/skin/eQTLs/chromatin/TF_peak_overlapping/QTL_1E-02"
f = data.table::fread(paste0(dir,"/QTL_overlapping_TF_peaks.bed"))
colnames(f) <- c("chr","start","end","SNP","REF","ALT","TF","ENCODE_ID")


# collapse the TF column for each unique SNP
data = f[,c("SNP","TF")] %>% distinct()

data.collapsed.n = data %>% 
  group_by(SNP) %>% 
  mutate(count = n()) %>% 
  distinct( . , SNP, .keep_all=TRUE) %>%
  select(-TF)

data.collapsed = data %>% 
  arrange(TF) %>%
  group_by(SNP) %>% 
  summarize(TF=paste(TF,collapse=','))

data.joined <- left_join(data.collapsed.n , data.collapsed, by="SNP")

# compile bed file and write out
bedf <- data.table::fread(paste0(dir,"/QTL.bed"))
colnames(bedf) <- c("chr","start","end","SNP","REF","ALT")
output <- left_join(bedf, data.joined, by="SNP")
colnames(output) <- c("chr","start","end","ID","REF","ALT","num_of_peak_overlapping_TF","peak_overlapping_TF")
output$num_of_peak_overlapping_TF[which(is.na(output$num_of_peak_overlapping_TF))] <- 0
data.table::fwrite(output, file=paste0(dir,"/QTL_overlapping_TF_peaks_collapsed.bed"), 
                   quote=F, sep="\t", row.names=F, col.names=T)
