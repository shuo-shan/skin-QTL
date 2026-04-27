library(dplyr)
library(magrittr)
library(tidyverse)

################# SET_UP
celltype="FRB"
Dir="~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method4"
setwd(Dir)
# read in the look up table of each regulatory region. (bed file)
bedF <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/regulatory_regions.bed")
colnames(bedF) <- c("chr","start","end","region")
# read in the DESeq2 result for ATACseq
dir.deseq <- "~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2"
atacF <- data.table::fread(paste0(dir.deseq,"/DESeq2_results_ATACseq_peaks_1kb_",celltype,".txt"))
colnames(atacF) <- c("region","atac_baseMean","atac_log2FoldChange","lfcSE","stat","pvalue","atac_padj","atac_dynamic")
# read in the DESeq2 result for H3K27ac ChIPseq
dir.deseq <- "~/Downloads/nl/human/skin/eQTLs/ChIP-seq/DESeq2"
chipF <- data.table::fread(paste0(dir.deseq,"/DESeq2_results_H3K27acChIPseq_peaks_1kb_",celltype,".txt"))
colnames(chipF) <- c("region","chip_baseMean","chip_log2FoldChange","lfcSE","stat","pvalue","chip_padj","chip_dynamic")

################# MAIN
# read in the gene-cRE link (dictionary_closest_gene_enhancer_links_FRB or dictionary_gene_promoter_links_FRB)
fname <- paste0("dictionary_closest_gene_enhancer_links_",celltype)
gene_cRE_link <- data.table::fread(paste0(Dir,"/",fname,".txt"))
colnames(gene_cRE_link) <- c("gene","region","dist","celltype")

# join the many tables
resultF1 <- right_join(bedF, gene_cRE_link, by="region")
resultF1$name <- paste(resultF1$region, resultF1$gene, resultF1$dist, resultF1$celltype, sep="_") %>% gsub("merged_all_skin-eQTL_ATACseq_files_peak_","peak",.)
resultF1 <- resultF1[,c("chr","start","end","name","region","gene","dist","celltype")]

resultF2 <- left_join(resultF1, atacF[,c("region","atac_baseMean","atac_log2FoldChange","atac_padj","atac_dynamic")], by="region")
resultF2$atac_dynamic[resultF2$atac_dynamic == ""] <- "none"

resultF3 <- left_join(resultF2, chipF[,c("region","chip_baseMean","chip_log2FoldChange","chip_padj","chip_dynamic")], by="region")
resultF3$chip_dynamic[resultF3$chip_dynamic == ""] <- "none"

# write to file
colnames(resultF3)[1] <- paste0("#",colnames(resultF3)[1]) # this step comments out the header row so bedtools could read it
data.table::fwrite(resultF3, paste0(Dir,"/",fname,".bed"), quote=F, sep="\t",col.names=T, row.names=F)
rm(gene_cRE_link, resultF1, resultF2, resultF3, atacF, chipF)
