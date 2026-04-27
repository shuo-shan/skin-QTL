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
  Rscript step14_transQTL_BHcorrection.R <infile>

Example:
  Rscript step14_transQTL_BHcorrection.R FRB_PBS_eQTL.bonferroni.txt", "\n")
  quit(save = "no", status = 1)
}

#infile <- "FRB_PBS_eQTL.bonferroni.txt"
infile         <- args[[1]]
basename <- unlist(strsplit(infile, ".bonferroni.txt"))
ct <- unlist(strsplit(basename, "_"))[1]
condition  <- unlist(strsplit(basename, "_"))[2]
QTLtype    <- unlist(strsplit(basename, "_"))[3]

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/transQTL/bonferroni/results")

# read file 
f <- fread(paste0(dir,"/",infile))

# performs BH multiple testing correction
f$q_gene = p.adjust(f$p_gene_bonferroni, method="BH")
f = f %>% arrange(q_gene,pmin)
f_gene_fdr05 <- f[which(f$q_gene<0.05),] %>% arrange(q_gene,pmin)

# write results that passed fdr < 0.05
outname <- paste0(basename,"_gene_fdr05_table.txt")  
data.table::fwrite(f_gene_fdr05, 
                   file = paste0(dir,"/",outname), 
                   quote=F, sep="\t")

write(f_gene_fdr05$gene, file=paste0(dir,"/",basename,"_gene_fdr05_genelist.txt"))

# overwrite bonferroni table by adding the BH-corrected q-value back
data.table::fwrite(
  f,
  file = paste0(dir, "/", infile),
  quote = FALSE,
  sep = "\t"
)

message(paste0("wrote output to ", paste0(dir,"/",outname)))
message(paste0("wrote output to ", paste0(dir,"/",basename,"_gene_fdr05_genelist.txt")))
message(paste0("wrote output to ", paste0(dir, "/", infile)))