library(dplyr)
library(magrittr)
# for a given TF of interest, find a list of highly correlated genes and a list of poorly correlated genes.
# all genes come from the RNAseq heatmap (DE genes).
# correlation method is: pearson correlation coefficient using z-score of log2(CPM+1)
# threshold is: abs(PCC) > 0.6 and abs(PCC) < 0.2

# first load log2FC1.5 heatmap .rdata
load("/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/myEnvironment_heatmap_RNAseq_DEgenes_14kmm-hc_04182023_padj0.05_log2FC1.5_avgCPM10.RData")

# pick DE genes and TFs
args = commandArgs(trailingOnly=TRUE)
TF=args[1] # RUNX1
celltype=args[2] # MEL, KRT, FRB, 3cts
Dir=paste0("/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_cobinding/results/",TF)

#### use zlogCPM to calculate pearson correlation coefficient
zlogCPM <- t(scale(t( log2( CPM_rna + 1 ))))
CPM_TF <- zlogCPM[TF,]

gene_list=c(rownames(zlogCPM.heatmap.new),"FOXA1")
CPM_picked_genes <- zlogCPM[gene_list,]

# generate correlation matrix
correlation_list <- c()
for (i in 1:length(gene_list)) {
  correlation_list[i] <- cor(as.numeric(CPM_picked_genes[i,]), as.numeric(CPM_TF) , method="pearson")
}
names(correlation_list) <- gene_list

# find highly correlated genes
write(names(which(abs(correlation_list)>0.6)), paste0(Dir,"/",TF,"_highly_correlated_genes.txt"))

# find lowly correlated genes
write(names(which(abs(correlation_list)<0.2)), 
      paste0(Dir,"/",TF,"_poorly_correlated_genes.txt"))


# #### archive don't use
# CPM_TF <- CPM_rna[TF,]
# 
# gene_list=c(rownames(zlogCPM.heatmap.new),"FOXA1")
# CPM_picked_genes <- CPM_rna[gene_list,]
# 
# # generate correlation matrix
# correlation_list <- c()
# for (i in 1:length(gene_list)) {
#   correlation_list[i] <- cor(as.numeric(CPM_picked_genes[i,]), as.numeric(CPM_TF) , method="spearman")
# }
# names(correlation_list) <- gene_list
# 
# # find highly correlated genes
# write(names(which(abs(correlation_list)>0.5)), paste0(Dir,"/",TF,"_highly_correlated_genes.txt"))
# 
# # find lowly correlated genes
# write(names(which(abs(correlation_list)<0.2)), 
#       paste0(Dir,"/",TF,"_poorly_correlated_genes.txt"))
