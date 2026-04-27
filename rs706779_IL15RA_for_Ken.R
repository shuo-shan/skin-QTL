#!/usr/bin/env Rscript
# written by Crystal Shan 09/2023
# This script makes a customized paired-plot for rs706779 IL15RA gene pair

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)

outDir="~/Downloads"
# load genotype
genotype_table="~/Downloads/nl/human/skin/eQTLs/DREG/master_filtered_genotype.bed"
genotype = data.table::fread(genotype_table,header=TRUE,sep="\t")
###### functions and loading data tables #####
rankNorm <- function (y) {
  # input y: numeric vector of CPM across all genes per donor.
  k <- 0.375 # an offset to ensure the z-score is finite. from Blom transform.
  n <- length(y)
  # Ranks.
  r <- rank(y, ties.method = "average") # if same value, same rank
  # Apply transformation.
  r.prob <- (r - k) / (n - 2 * k + 1)
  y.rankNorm <- stats::qnorm(r.prob) # because qnorm(1) is infinite and qnorm(-1) is -Inf
  return(y.rankNorm)
}
###### load-data-for-plotting ###### 
# load MEL data
CPM_table_PBS.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/CPM_expressedGenes_IFN.txt"
PBS_knownVar_file.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/metadata_PBS.txt"
IFN_knownVar_file.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/metadata_IFN.txt"
CPM.pbs.mel <- read.table(CPM_table_PBS.mel,sep="\t",header=TRUE) %>% column_to_rownames("X")
CPM.ifn.mel <- read.table(CPM_table_IFN.mel,sep="\t",header=TRUE) %>% column_to_rownames("X")
PBS_knownVar.mel <- read.table(PBS_knownVar_file.mel, header=TRUE,sep="\t") 
IFN_knownVar.mel <- read.table(IFN_knownVar_file.mel, header=TRUE,sep="\t") 
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

# column names for result_table
result_table_colnames <- c()
result_table_colnames[c(1,2)] <- c("SNP","GENE")
result_table_colnames[c(3,4,5,6)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("log2FC_minimalModel_reQTL_genotype_",x)))
result_table_colnames[c(7,8,9,10)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("log2FC_featureSelected_reQTL_genotype_",x)))
result_table_colnames[c(11,12,13,14)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormlog2FC_minimalModel_reQTL_genotype_",x)))
result_table_colnames[c(15,16,17,18)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormlog2FC_featureSelected_reQTL_genotype_",x)))

result_table_colnames[c(19,20,21,22)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_minimalModel_PBSeQTL_genotype_",x)))
result_table_colnames[c(23,24,25,26)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_featureSelected_PBSeQTL_genotype_",x)))
result_table_colnames[c(27,28,29,30)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_minimalModel_PBSeQTL_genotype_",x)))
result_table_colnames[c(31,32,33,34)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_featureSelected_PBSeQTL_genotype_",x)))

result_table_colnames[c(35,36,37,38)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_minimalModel_IFNeQTL_genotype_",x)))
result_table_colnames[c(39,40,41,42)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_featureSelected_IFNeQTL_genotype_",x)))
result_table_colnames[c(43,44,45,46)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_minimalModel_IFNeQTL_genotype_",x)))
result_table_colnames[c(47,48,49,50)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_featureSelected_IFNeQTL_genotype_",x)))

###### plot paired-plot between rs706779 and IL15RA. ######
# combine the donors of TT and TC into one panel, 
# and CC in a separate panel.
make_reQTL_plot_CPM <- function(this.snp,this.gene,CPM.pbs,CPM.ifn, celltype) {
  
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0","TT_or_TC", . ) %>%
    gsub("0/1","TT_or_TC", . ) %>%
    gsub("1/1","CC", . ) %>%
    gsub("./.",NA, . )
  this.genotype$genotype <- ordered(this.genotype$genotype, 
                                    levels <- c("TT_or_TC",
                                                "CC"))
  this.genotype = this.genotype %>% na.omit()
  rm(this.snp.genotype)
  
  # pick the gene
  this.CPM.PBS <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
  this.CPM.IFN <- CPM.ifn[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","IFN"))
  this.CPM <- inner_join(this.CPM.PBS, this.CPM.IFN, by="donor") %>% 
    inner_join( . , this.genotype, by="donor")
  
  # make paired plots
  position <- "identity"; width = 0.5; point.size = 1; line.size = 0.5
  line.color = "grey"; linetype = "solid"; palette="bright"
  
  df <- this.CPM %>% pivot_longer(. , c("PBS","IFN"),names_to="condition",values_to = "CPM")
  df$condition <- ordered(df$condition, levels=c("PBS","IFN"))
  highlightdata1 <- df %>% filter(donor=="F25")
  highlightdata2 <- df %>% filter(donor=="F49")
  highlightdata3 <- df %>% filter(donor=="F55")
  genotype_counts <- df %>%
    group_by(genotype) %>%
    summarise(donor_count = n()) %>%
    mutate(label = paste0(genotype, "\n(", donor_count/2, ")"))
  facet_labels <- setNames(genotype_counts$label, genotype_counts$genotype)
  df$genotype <- factor(df$genotype, levels = names(facet_labels))
  
  p <- ggplot(df, create_aes(list(x = "condition", y = "CPM"))) +
    geom_exec(geom_line, data = df, group = "donor",
              color = line.color, size = line.size, linetype = linetype,
              position = position) +
    geom_exec(geom_point, data = df, color = "condition", size = point.size,
              position = position) +
    scale_color_manual(values=c("#999999", "#E69F00")) +
    geom_exec(geom_line, data = highlightdata1, group = "donor",
              color = "blue", size = 0.5, linetype = linetype,
              position = position) +
    geom_exec(geom_line, data = highlightdata2, group = "donor",
              color = "brown", size = 0.5, linetype = linetype,
              position = position) +
    geom_exec(geom_line, data = highlightdata3, group = "donor",
              color = "black", size = 0.5, linetype = linetype,
              position = position) +
    ggtitle(celltype) +
    theme(legend.position="none", 
          plot.background = element_rect(colour = "grey", fill=NA, linewidth=0.5), 
          axis.line = element_line(colour = "grey"),
          panel.background = element_rect(fill=NA),
          panel.spacing.x = unit(1, "points"),
          plot.title = element_text(margin=margin(0,0,0,0))) +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=8),
          plot.title = element_text(size =8, face = "bold"),
          strip.text = element_text(size =8)) +
    facet_wrap(~ genotype, labeller = as_labeller(facet_labels))
  return(p)
}
make_reQTL_plot_CPM_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.mel, CPM.ifn.mel,"MEL")
  plot.krt <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.krt, CPM.ifn.krt,"KRT")
  plot.frb <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.frb, CPM.ifn.frb,"FRB")
  
  combined_plot <- plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3)
  title_plot <- ggdraw() + draw_label(paste0(this.snp,":",this.gene," CPM"), size=10,fontface = 'bold', x = 0.01, hjust = 0)
  
  final_plot <- plot_grid(title_plot, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
  return(final_plot)
}

make_reQTL_plot_CPM_PBSonly <- function(this.snp,this.gene,CPM.pbs,CPM.ifn, celltype) {
  
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0","TT_or_TC", . ) %>%
    gsub("0/1","TT_or_TC", . ) %>%
    gsub("1/1","CC", . ) %>%
    gsub("./.",NA, . )
  this.genotype$genotype <- ordered(this.genotype$genotype, 
                                    levels <- c("TT_or_TC",
                                                "CC"))
  this.genotype = this.genotype %>% na.omit()
  rm(this.snp.genotype)
  
  # pick the gene
  this.CPM.PBS <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
  this.CPM <- inner_join( this.CPM.PBS , this.genotype, by="donor")
  
  # make paired plots
  position <- "identity"; width = 0.5; point.size = 3; line.size = 0.5
  line.color = "grey"; linetype = "solid"; palette="bright"
  
  df <- this.CPM %>% pivot_longer(. , "PBS",names_to="condition",values_to = "CPM")
  highlightdata1 <- df %>% filter(donor=="F25")
  highlightdata2 <- df %>% filter(donor=="F49")
  highlightdata3 <- df %>% filter(donor=="F55")
  genotype_counts <- df %>%
    group_by(genotype) %>%
    summarise(donor_count = n()) %>%
    mutate(label = paste0(genotype, "\n(", donor_count, ")"))
  facet_labels <- setNames(genotype_counts$label, genotype_counts$genotype)
  df$genotype <- factor(df$genotype, levels = names(facet_labels))
  
  p <- ggplot(df, create_aes(list(x = "condition", y = "CPM"))) +
    geom_exec(geom_line, data = df, group = "donor",
              color = line.color, size = line.size, linetype = linetype,
              position = position) +
    geom_exec(geom_jitter, data = df, color = "condition", size = point.size,
              width=0.1, alpha=0.8) +
    scale_color_manual(values="#999999") +
    geom_exec(geom_point, data = highlightdata1, group = "donor",
              color = "blue", size = 3,
              position = position) +
    geom_exec(geom_point, data = highlightdata2, group = "donor",
              color = "brown", size = 3,
              position = position) +
    geom_exec(geom_point, data = highlightdata3, group = "donor",
              color = "black", size = 3,
              position = position) +
    ggtitle(celltype) +
    theme(legend.position="none", 
          plot.background = element_rect(colour = "grey", fill=NA, linewidth=0.5), 
          axis.line = element_line(colour = "grey"),
          panel.background = element_rect(fill=NA),
          panel.spacing.x = unit(1, "points"),
          plot.title = element_text(margin=margin(0,0,0,0))) +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=8),
          plot.title = element_text(size =8, face = "bold"),
          strip.text = element_text(size =8)) +
    facet_wrap(~ genotype, labeller = as_labeller(facet_labels))
  return(p)
}
make_reQTL_plot_CPM_PBSonly_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_reQTL_plot_CPM_PBSonly(this.snp,this.gene,CPM.pbs.mel, CPM.ifn.mel,"MEL")
  plot.krt <- make_reQTL_plot_CPM_PBSonly(this.snp,this.gene,CPM.pbs.krt, CPM.ifn.krt,"KRT")
  plot.frb <- make_reQTL_plot_CPM_PBSonly(this.snp,this.gene,CPM.pbs.frb, CPM.ifn.frb,"FRB")
  
  combined_plot <- plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3)
  title_plot <- ggdraw() + draw_label(paste0(this.snp,":",this.gene," CPM"), size=10,fontface = 'bold', x = 0.01, hjust = 0)
  
  final_plot <- plot_grid(title_plot, combined_plot, ncol = 1, rel_heights = c(0.1, 1))
  
  return(final_plot)
}
###### make plots ######
pdf(paste0(outDir,"/pairedPlot_3donorColored_CPM_PBSonly_","rs706779","_","IL15RA",".pdf"))
make_reQTL_plot_CPM_PBSonly_3cts("rs706779","IL15RA")
dev.off()
