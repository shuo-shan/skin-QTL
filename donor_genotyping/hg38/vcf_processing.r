#!/bin/Rscript
# author: Shuo Shan
# July 2021
# turn phased genotype data into matrix suitable for clustering

library(dplyr)
library(tidyverse)
library(BiocParallel)
library(magrittr)
library(stringr)
library(vcfR)
library(cluster) # for gower similarity and pam
library(ggplot2) # for visualization
dir="/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38"
setwd("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/")
Sys.time()
vcf <- read.vcfR(paste0(dir,"/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.16donors.filtered2.pruned.vcf.gz"))
vcf_table <- as_tibble(cbind(vcf@fix,vcf@gt))
#write.table(vcf_table, file=paste0(dir,"/vcf/bd841628-fcc2-487a-8460-f5428237f0c9.13donors.phased_vcf_table.txt"), row.names=T)
write("loaded vcf", stdout())
Sys.time()

##################
#vcf_table2 <- mutate_if(vcf_table, is.character, str_replace_all, pattern="1\\|1", replacement="HomoAlt")
#vcf_table2 <- mutate_if(vcf_table2, is.character, str_replace_all, pattern="0\\|0", replacement="HomoRef")
#vcf_table2 <- mutate_if(vcf_table2, is.character, str_replace_all, pattern="1\\|0", replacement="Het")
#vcf_table2 <- mutate_if(vcf_table2, is.character, str_replace_all, pattern="0\\|1", replacement="Het")
#write.table(vcf_table2, file=paste0(dir,"/vcf/bd841628-fcc2-487a-8460-f5428237f0c9.13donors.phased_vcf_table2.txt"), row.names=T)
#write("compiled vcf_table2", stdout())
#Sys.time()
#
#################
vcf_table3 <- mutate_if(vcf_table, is.character, str_replace_all, pattern="1\\|1", replacement="1")
vcf_table3 <- mutate_if(vcf_table3, is.character, str_replace_all, pattern="0\\|0", replacement="0")
vcf_table3 <- mutate_if(vcf_table3, is.character, str_replace_all, pattern="1\\|0", replacement="0.5")
vcf_table3 <- mutate_if(vcf_table3, is.character, str_replace_all, pattern="0\\|1", replacement="0.5")
write.table(vcf_table3, file=paste0(dir,"/vcf/bd841628-fcc2-487a-8460-f5428237f0c9.16donors.pruned.phased_vcf_table3.txt"), row.names=T, quote=F)
write("compiled vcf_table3", stdout())
Sys.time()

##############
write("start processing", stdout())
Sys.time()
#vcf_table3 <- read.table(file=paste0(dir,"/13donors.vcf_table3.txt"), sep="\t")
#write("loading done", stdout())
#Sys.time()
# transpose vcf2_table_select
vcf_table3_transposed <- as_tibble(cbind(names=colnames(vcf_table3),t(vcf_table3)))
# turn values into numeric
vcf_table3_transposed[,2:ncol(vcf_table3_transposed)] <- 
  apply(vcf_table3_transposed[,2:ncol(vcf_table3_transposed)], 2, function(x) as.numeric(as.character(x)))
write("transposing done", stdout())
Sys.time()
manhattan_dist <- daisy(vcf_table3_transposed[, -1], metric = "manhattan")
write("manhattan dist calculation done", stdout())
Sys.time()
aggl.clust.c.manhattan <- hclust(manhattan_dist, method = "complete")
write("clustering done", stdout())
Sys.time()
aggl.clust.c.manhattan$labels=as.character(vcf_table3_transposed$names)
manhattan_cluster_cut3<- cutree(aggl.clust.c.manhattan, k=3)
donors_manhattan_cluster_cut3 <- tibble(vcf_table3_transposed$names,manhattan_cluster_cut3) %>% arrange(manhattan_cluster_cut3)
print(donors_manhattan_cluster_cut3, n=Inf)
write.table(donors_manhattan_cluster_cut3, file=paste0(dir,"/16donors.donors_manhattan_cluster_cut3.txt"), quote=F, sep="\t")

manhattan_cluster_cut3<- cutree(aggl.clust.c.manhattan, k=4)
donors_manhattan_cluster_cut3 <- tibble(vcf_table3_transposed$names,manhattan_cluster_cut3) %>% arrange(manhattan_cluster_cut3)
print(donors_manhattan_cluster_cut3, n=Inf)
write.table(donors_manhattan_cluster_cut3, file=paste0(dir,"/16donors.donors_manhattan_cluster_cut4.txt"), quote=F, sep="\t")

pdf(paste0(dir,"/16donors.manhattan.clust.complete.pdf")) 
plot(aggl.clust.c.manhattan, main = "Agglomerative, complete linkages")
dev.off()
write("plotting done", stdout())
Sys.time()
