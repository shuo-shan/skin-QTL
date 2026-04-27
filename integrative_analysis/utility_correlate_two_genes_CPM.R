library(dplyr)
library(magrittr)
library(tidyverse)

# first load log2FC1.5 heatmap .rdata
Dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis"
load(paste0(Dir,"/heatmap_RNA/myEnvironment_7kmm-hc_padj0.05_log2FC1.5_avgCPM10.RData"))

# pick DE genes and TFs
TF="RUNX1"
CPM_TF <- CPM_rna[TF,]

gene_list=c(cluster1_genes, cluster2_genes, cluster3_genes, 
            cluster4_genes, cluster5_genes, cluster6_genes, cluster7_genes,"FOXA1")
CPM_picked_genes <- CPM_rna[gene_list,]

# generate correlation matrix
correlation_list <- c()
for (i in 1:length(gene_list)) {
  correlation_list[i] <- cor(as.numeric(CPM_picked_genes[i,]), as.numeric(CPM_TF) , method="pearson")
}
names(correlation_list) <- gene_list

# find highly correlated genes
names(which(correlation_list>0.5))
write(names(which(correlation_list>0.5)), paste0(Dir,"/FIMO/",TF,"_highly_correlated_genes.txt"))

# find lowly correlated genes
names(which(correlation_list>0 & correlation_list<0.2))
write(names(which(correlation_list>0 & correlation_list<0.2)), 
      paste0(Dir,"/FIMO/",TF,"_poorly_correlated_genes.txt"))

#names(which(correlation_list < -0))
#write(names(which(correlation_list < -0)), paste0(Dir,"/FIMO/",TF,"_negatively_correlated_genes.txt"))



induced.genes=rownames(CPM_rna_DEG)
non.induced.genes=rownames(CPM_rna_new)[which(!rownames(CPM_rna_new) %in% rownames(CPM_rna_DEG))]
set.seed(42)
non.induced.genes.sampled=sample(non.induced.genes,length(induced.genes))

write(induced.genes, paste0(Dir,"/FIMO/IFNg_induced_genes.txt"))
write(non.induced.genes.sampled, paste0(Dir,"/FIMO/non_IFNg_induced_genes_number_matched.txt"))
