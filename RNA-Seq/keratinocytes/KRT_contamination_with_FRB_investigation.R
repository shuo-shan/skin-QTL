#!/usr/bin/env Rscript
# written by shuo.shan@umassmed.edu 04/2024

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)


Dir="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes"
###### load-data-for-plotting ###### 
# load KRT data
CPM_table_PBS.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt"
PBS_knownVar_file.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/metadata_PBS.txt"
IFN_knownVar_file.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/metadata_IFN.txt"
CPM.pbs.krt <- read.table(CPM_table_PBS.krt,header=TRUE,sep="\t") %>% column_to_rownames("X")
CPM.ifn.krt <- read.table(CPM_table_IFN.krt,header=TRUE,sep="\t") %>% column_to_rownames("X")
PBS_knownVar.krt <- read.table(PBS_knownVar_file.krt, sep="\t",header=TRUE) 
IFN_knownVar.krt <- read.table(IFN_knownVar_file.krt, sep="\t",header=TRUE) 
# load FRB data
CPM_table_PBS.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt"
PBS_knownVar_file.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/metadata_PBS.txt"
IFN_knownVar_file.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/metadata_IFN.txt"
CPM.pbs.frb <- read.table(CPM_table_PBS.frb,header=TRUE,sep="\t") %>% column_to_rownames("X")
CPM.ifn.frb <- read.table(CPM_table_IFN.frb,header=TRUE,sep="\t") %>% column_to_rownames("X")
PBS_knownVar.frb <- read.table(PBS_knownVar_file.frb, header=TRUE,sep="\t") 
IFN_knownVar.frb <- read.table(IFN_knownVar_file.frb, header=TRUE,sep="\t") 

###### plot FRB-specific genes ######
gene_list <- c("PDGFRA","COL1A1","COL3A1","THY1")
# Initialize an empty list to store plots
plots_list <- list()
# Iterate over each gene to create individual plots
for (g in gene_list) {
  df.paired.krt <- data.frame(PBS=as.data.frame(t(CPM.pbs.krt[g,])),
                              IFN=as.data.frame(t(CPM.ifn.krt[g,])),
                              celltype="keratinocytes") %>% 
    set_colnames(c("PBS","IFN","celltype"))
  df.paired.frb <- data.frame(PBS=as.data.frame(t(CPM.pbs.frb[g,])),
                              IFN=as.data.frame(t(CPM.ifn.frb[g,])),
                              celltype="fibroblasts") %>% 
    set_colnames(c("PBS","IFN","celltype"))
  
  df.paired <- rbind(df.paired.krt, df.paired.frb)
  
  # Generate plot for current gene
  p <- ggpaired(df.paired, cond1="PBS", cond2="IFN", 
                color="celltype", palette=c("grey","red"),
                xlab="Condition", ylab="CPM",
                title = g) +
    theme(legend.position = "none") +
    facet_wrap(~celltype)
  
  # Add the plot to the list
  plots_list[[g]] <- p
}

# Save all plots in one PDF page
pdf(paste0(Dir,"/RNAseq_plot_FRB_specific_genes_combined.pdf"))
grid.arrange(grobs = plots_list, ncol = 2) # Adjust ncol to fit your preference
dev.off()