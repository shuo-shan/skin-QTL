# decide gene chunks basesd on SNP:gene pairs
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("
Usage:
  Rscript step0_decide_gene_chunks.R <ct>

Example:
  Rscript step0_decide_gene_chunks.R FRB", "\n")
  quit(save = "no", status = 1)
}

ct <- args[[1]]
#ct <- "FRB"
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)

# read file 
f <- fread(paste0(dir,"/data/SNPs_near_TSS.bed"))

# summarize number of SNPs per gene
summary <- table(f$gene_name) %>% 
  as.data.frame() %>%
  set_colnames(c("gene","nSNP"))
rm(f)


# decide gene chunks
threshold <- 100000
i <- 1
chunks <- list()

df <- summary

while (nrow(df) > 0) {
  cs <- cumsum(df$nSNP)
  
  j_exceed <- which(cs > threshold)[1]
  
  if (is.na(j_exceed)) {
    idx <- seq_len(nrow(df))
  } else {
    idx <- seq_len(j_exceed - 1)
    if (length(idx) == 0) idx <- 1
  }
  
  gene_chunk_i <- df$gene[idx]
  
  chunk_name <- sprintf("gene_chunk_%03d", i)
  chunks[[chunk_name]] <- gene_chunk_i
  
  df <- df[-idx, , drop = FALSE]
  i <- i + 1
}

# Write to file
for (i in 1:length(chunks)) {
  this_name <- sprintf("gene_chunk_%03d",i)
  this_gene_chunk <- chunks[[this_name]]
  write(this_gene_chunk, file=paste0(dir,"/chunks/",this_name))
}





