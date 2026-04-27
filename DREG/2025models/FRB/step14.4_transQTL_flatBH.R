# performs Benjamini Hochberg multiple testing correction across all genomic loci using bonferroni locally corrected p.values
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
  Rscript step14_transQTL_flatBH.R <infile>

Example:
  Rscript step14_transQTL_flatBH.R FRB_PBS_eQTL.result.txt", "\n")
  quit(save = "no", status = 1)
}

#infile <- "FRB_PBS_eQTL.result.txt"
infile         <- args[[1]]
basename <- unlist(strsplit(infile, ".result.txt"))
ct <- unlist(strsplit(basename, "_"))[1]
condition  <- unlist(strsplit(basename, "_"))[2]
QTLtype    <- unlist(strsplit(basename, "_"))[3]

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/transQTL/resultsBHcorrected")

# read file 
f <- fread(paste0(dir,"/",infile))

# performs BH multiple testing correction
f$q_gene = p.adjust(f$p, method="BH")
f = f %>% arrange(q_gene,p)
f_gene_fdr05 <- f[which(f$q_gene<0.05),] %>% arrange(q_gene,p)

if (nrow(f_gene_fdr05)>0){
  # write results that passed fdr < 0.05
  outname <- paste0(basename,"_gene_fdr05_table.txt")  
  data.table::fwrite(f_gene_fdr05, 
                     file = paste0(dir,"/",outname), 
                     quote=F, sep="\t")
  
  write(f_gene_fdr05$gene, file=paste0(dir,"/",basename,"_gene_fdr05_genelist.txt"))
  
  message(paste0("wrote output to ", paste0(dir,"/",outname)))
  message(paste0("wrote output to ", paste0(dir,"/",basename,"_gene_fdr05_genelist.txt")))
}

# overwrite concatenated result table by adding the BH-corrected q-value back
data.table::fwrite(
  f,
  file = paste0(dir, "/", infile),
  quote = FALSE,
  sep = "\t"
)


message(paste0("wrote output to ", paste0(dir, "/", infile)))
