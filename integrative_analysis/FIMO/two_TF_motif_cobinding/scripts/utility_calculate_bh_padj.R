#!/usr/bin/env Rscript
library(dplyr)
library(magrittr)
# calculate BH p.adjust and fill the existing column q-value in the same input file.

args = commandArgs(trailingOnly=TRUE)
file = args[1]

#file="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/FIMO/GATA2/GATA2_GATA2HighCorrelationGenesEnhancers/fimo_GATA2.bed"
#file="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_cobinding/results/GATA2/GATA2_GATA2HighCorrelationGenesPromoters/fimo_all_output.txt"

dir=gsub("/fimo_.*","",file)
fname=gsub(".*/","",file)
df = data.table::fread(file)
df$`q-value` <- p.adjust(as.numeric(df$`p-value`), method="BH")
df = df %>% arrange(`q-value`)
df$`p-value` = format(df$`p-value`, scientific = F)
df$`q-value` = format(df$`q-value`, scientific = F)
data.table::fwrite(df, file=paste0(dir,"/",fname),quote=F,sep="\t",row.name=F,col.names=T)
