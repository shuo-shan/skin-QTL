#!/usr/bin/env Rscript
library(dplyr)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
file = args[1]

#file="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/FIMO/GATA2/GATA2_GATA2HighCorrelationGenesEnhancers/fimo_GATA2.bed"
dir=gsub("/fimo_.*","",file)
fname=gsub(".*/","",file)
df = data.table::fread(file)
df$padj <- p.adjust(df$V8, method="BH")
data.table::fwrite(df, file=paste0(dir,"/padj_",fname),quote=F,sep="\t",row.name=F,col.names=F)
