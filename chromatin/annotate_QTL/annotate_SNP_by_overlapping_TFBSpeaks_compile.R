# written by shuo.shan@umassmed.edu 03/2024
library(tidyverse)
library(magrittr)

# take in arguments
args <- commandArgs(trailingOnly = TRUE)
workingdir <- args[1]
inFile <- args[2]
prefix <- args[3]

# load QTL and TF CHIPseq intersect result
f <- data.table::fread(inFile)

if (nrow(f) > 0) {
  colnames(f) <- c("chr","start","end","SNP","REF","ALT","TF","ENCODE_ID")
  
  # collapse the TF column for each unique SNP
  data <- f[,c("SNP","TF")] %>% distinct()
  
  data.collapsed.n <- data %>% 
    group_by(SNP) %>% 
    mutate(n_TFBS = n()) %>% 
    distinct(., SNP, .keep_all = TRUE) %>%
    select(-TF)
  
  data.collapsed <- data %>% 
    arrange(TF) %>%
    group_by(SNP) %>% 
    summarize(TFBS_list = paste(TF, collapse = ','))
  
  data.joined <- left_join(data.collapsed.n, data.collapsed, by = "SNP")
  
} else {
  data.joined <- tibble(SNP = NA_character_, n_TFBS = 0, TFBS_list = "")
}

# write out
outF=paste0(workingdir, "/QTL_overlapping_TF_peaks_collapsed_", prefix, ".bed")
data.table::fwrite(data.joined, 
                   file = outF,
                   quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

cat(paste0("wrote output to ", outF,"\n"))