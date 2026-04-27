#!/bin/Rscript
library(dplyr)
library(tidyverse)

Dir="/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL"
df <- read.table(file=paste0(Dir,"/filtered.dosage.txt"),sep="\t", header=TRUE) 
dosage <- df[,-c(1:2)] %>% dplyr::filter(duplicated(ID) == FALSE) %>% column_to_rownames(.,"ID")

#calculate principal components
results <- prcomp(dosage, scale = TRUE)
#reverse the signs
results$rotation <- -1*results$rotation
#display principal components
top.5.pc <- t(results$rotation[,c(1:5)]) %>% as.data.frame()
#write to file
write.table(top.5.pc,file=paste0(Dir,"/genotype_dosage_topPCs.txt"),quote=F,sep="\t",row.names=T,col.names=T)
