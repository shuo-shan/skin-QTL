#!/usr/bin/env Rscript
# written by Crystal Shan 03/2022. modified on 01/2023.

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
library(LDlinkR)

outDir="~/Downloads/nl/human/skin/eQTLs/DREG/MEL_minimal"
resultF="~/Downloads/nl/human/skin/eQTLs/DREG/MEL_minimal/MEL_PBSeQTL_pval1E-05.txt"
n_gPCs=0
n_ePCs=0
genotype_table="~/Downloads/nl/human/skin/eQTLs/DREG/MEL_minimal/snps_near_expressed_genes.bed"

#setwd(outDir)
# functions and loading data tables (skip if myEnvironment.RData is loaded)
# useful functions
make_PBSeQTL_plot_CPM <- function(this.snp,this.gene,CPM.pbs,PBS_knownVar) {
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(snp.ref,snp.ref), . ) %>%
    gsub("0/1",paste0(snp.ref,snp.alt), . ) %>%
    gsub("1/1",paste0(snp.alt,snp.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype = this.genotype %>% na.omit()
  rm(this.snp.genotype)
  
  # pick the gene
  this.CPM.PBS <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS")) %>% 
    inner_join( . , this.genotype, by="donor")
  this.CPM.PBS$genotype <- ordered(this.CPM.PBS$genotype, 
                                   levels <- c(paste0(snp.ref,snp.ref),
                                               paste0(snp.ref,snp.alt),
                                               paste0(snp.alt,snp.alt)))
  this.CPM.PBS <- left_join(this.CPM.PBS, PBS_knownVar, by="donor")
  # make scatter plots
  p1<-ggplot(this.CPM.PBS, aes(x=genotype, y=PBS, shape=seq_model, color=ancestry)) +
    geom_jitter(size=3,width=0.1) +
    ggtitle(paste0(this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM")) +
    scale_color_brewer(palette="Dark2")
  p1
  return(p1)
}
make_IFNeQTL_plot_CPM <- function(this.snp,this.gene,CPM.ifn,IFN_knownVar) {
  
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(snp.ref,snp.ref), . ) %>%
    gsub("0/1",paste0(snp.ref,snp.alt), . ) %>%
    gsub("1/1",paste0(snp.alt,snp.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype = this.genotype %>% na.omit()
  rm(this.snp.genotype)
  
  # pick the gene
  this.CPM.IFN <- CPM.ifn[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","IFN")) %>% 
    inner_join( . , this.genotype, by="donor")
  this.CPM.IFN$genotype <- ordered(this.CPM.IFN$genotype, 
                                   levels <- c(paste0(snp.ref,snp.ref),
                                               paste0(snp.ref,snp.alt),
                                               paste0(snp.alt,snp.alt)))
  this.CPM.IFN <- left_join(this.CPM.IFN, IFN_knownVar, by="donor")
  # make scatter plots
  p1<-ggplot(this.CPM.IFN, aes(x=genotype, y=IFN, shape=seq_model, color=ancestry)) +
    geom_jitter(size=3,width=0.1) +
    ggtitle(paste0(this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM")) +
    scale_color_brewer(palette="Dark2")
  p1
  return(p1)
}
make_PBSeQTL_plot_CPM_nocolor <- function(this.snp,this.gene,CPM.pbs,PBS_knownVar) {
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(snp.ref,snp.ref), . ) %>%
    gsub("0/1",paste0(snp.ref,snp.alt), . ) %>%
    gsub("1/1",paste0(snp.alt,snp.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype = this.genotype %>% na.omit()
  rm(this.snp.genotype)
  
  # pick the gene
  this.CPM.PBS <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS")) %>% 
    inner_join( . , this.genotype, by="donor")
  this.CPM.PBS$genotype <- ordered(this.CPM.PBS$genotype, 
                                   levels <- c(paste0(snp.ref,snp.ref),
                                               paste0(snp.ref,snp.alt),
                                               paste0(snp.alt,snp.alt)))
  highlightdata0 <- this.CPM.PBS %>% filter(!donor %in% c("F25","F49","F55"))
  highlightdata1 <- this.CPM.PBS %>% filter(donor=="F25")
  highlightdata2 <- this.CPM.PBS %>% filter(donor=="F49")
  highlightdata3 <- this.CPM.PBS %>% filter(donor=="F55")
  
  # make scatter plots
  p1<-ggplot(this.CPM.PBS, aes(x=genotype, y=PBS)) +
    geom_exec(geom_point, data = highlightdata0, group = "donor",
              color = "#5A5A5A", size = 3) +
    geom_exec(geom_point, data = highlightdata1, group = "donor",
              color = "blue", size = 4) +
    geom_exec(geom_point, data = highlightdata2, group = "donor",
              color = "brown", size = 4) +
    geom_exec(geom_point, data = highlightdata3, group = "donor",
              color = "black", size = 4) +
    ggtitle(paste0(this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM"))
  p1
  return(p1)
}
make_IFNeQTL_plot_CPM_nocolor <- function(this.snp,this.gene,CPM.ifn,IFN_knownVar) {
  
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(snp.ref,snp.ref), . ) %>%
    gsub("0/1",paste0(snp.ref,snp.alt), . ) %>%
    gsub("1/1",paste0(snp.alt,snp.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype = this.genotype %>% na.omit()
  rm(this.snp.genotype)
  
  # pick the gene
  # pick the gene
  this.CPM.IFN <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS")) %>% 
    inner_join( . , this.genotype, by="donor")
  this.CPM.IFN$genotype <- ordered(this.CPM.IFN$genotype, 
                                   levels <- c(paste0(snp.ref,snp.ref),
                                               paste0(snp.ref,snp.alt),
                                               paste0(snp.alt,snp.alt)))
  highlightdata0 <- this.CPM.IFN %>% filter(!donor %in% c("F25","F49","F55"))
  highlightdata1 <- this.CPM.IFN %>% filter(donor=="F25")
  highlightdata2 <- this.CPM.IFN %>% filter(donor=="F49")
  highlightdata3 <- this.CPM.IFN %>% filter(donor=="F55")
  
  # make scatter plots
  p1<-ggplot(this.CPM.IFN, aes(x=genotype, y=PBS)) +
    geom_exec(geom_point, data = highlightdata0, group = "donor",
              color = "#5A5A5A", size = 3) +
    geom_exec(geom_point, data = highlightdata1, group = "donor",
              color = "blue", size = 4) +
    geom_exec(geom_point, data = highlightdata2, group = "donor",
              color = "brown", size = 4) +
    geom_exec(geom_point, data = highlightdata3, group = "donor",
              color = "black", size = 4) +
    ggtitle(paste0(this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM"))
  p1
  return(p1)
}
make_reQTL_plot_CPM <- function(this.snp,this.gene,CPM.pbs,CPM.ifn) {
  
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(snp.ref,snp.ref), . ) %>%
    gsub("0/1",paste0(snp.ref,snp.alt), . ) %>%
    gsub("1/1",paste0(snp.alt,snp.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype$genotype <- ordered(this.genotype$genotype, 
                                    levels <- c(paste0(snp.ref,snp.ref),
                                                paste0(snp.ref,snp.alt),
                                                paste0(snp.alt,snp.alt)))
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
  position <- "identity"; width = 0.5; point.size = 2; line.size = 0.5
  line.color = "grey"; linetype = "solid"; palette="bright"
  
  df <- this.CPM %>% pivot_longer(. , c("PBS","IFN"),names_to="condition",values_to = "CPM")
  df$condition <- ordered(df$condition, levels=c("PBS","IFN"))
  highlightdata1 <- df %>% filter(donor=="F25")
  highlightdata2 <- df %>% filter(donor=="F49")
  highlightdata3 <- df %>% filter(donor=="F55")
  
  p <- ggplot(df, create_aes(list(x = "condition", y = "CPM"))) +
    geom_exec(geom_line, data = df, group = "donor",
              color = line.color, size = line.size, linetype = linetype,
              position = position) +
    geom_exec(geom_point, data = df, color = "condition", size = point.size,
              position = position) +
    scale_color_manual(values=c("#999999", "#E69F00")) +
    geom_exec(geom_line, data = highlightdata1, group = "donor",
              color = "blue", size = 0.7, linetype = linetype,
              position = position) +
    geom_exec(geom_line, data = highlightdata2, group = "donor",
              color = "brown", size = 0.7, linetype = linetype,
              position = position) +
    geom_exec(geom_line, data = highlightdata3, group = "donor",
              color = "black", size = 0.7, linetype = linetype,
              position = position) +
    ggtitle(paste0(this.snp,":",this.gene," CPM")) +
    theme(legend.position="none") +
    facet_wrap(~ genotype)
  
  return(p)
}
make_paired_plot_CPM_3cts <- function(this.gene) {
  # p1: MEL
  this.CPM.PBS <- CPM.pbs.mel[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
  this.CPM.IFN <- CPM.ifn.mel[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","IFN"))
  this.CPM.mel <- inner_join(this.CPM.PBS, this.CPM.IFN, by="donor") %>% mutate(celltype=rep("MEL",nrow(this.CPM.PBS)))
  # p1: KRT
  this.CPM.PBS <- CPM.pbs.krt[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
  this.CPM.IFN <- CPM.ifn.krt[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","IFN"))
  this.CPM.krt <- inner_join(this.CPM.PBS, this.CPM.IFN, by="donor") %>% mutate(celltype=rep("KRT",nrow(this.CPM.PBS)))
  # p1: FRB
  this.CPM.PBS <- CPM.pbs.frb[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
  this.CPM.IFN <- CPM.ifn.frb[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","IFN"))
  this.CPM.frb <- inner_join(this.CPM.PBS, this.CPM.IFN, by="donor") %>% mutate(celltype=rep("FRB",nrow(this.CPM.PBS)))
  this.CPM <- rbind(this.CPM.mel, this.CPM.krt, this.CPM.frb)
  this.CPM$celltype <- ordered(this.CPM$celltype, levels=c("MEL","KRT","FRB"))
  
  # make paired plots
  p1<-ggpaired(this.CPM,cond1="PBS",cond2="IFN",facet.by="celltype",
               color="condition",palette="aaas",line.color = "gray", line.size = 0.4,
               title=this.gene,width=0,
               xlab="condition",ylab=paste0(this.gene," CPM"))
  p1$layers<-p1$layers[-1] # remove geom_boxplot layer 
  p1
  return(p1)  
}
make_reQTL_plot_CPM_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.mel, CPM.ifn.mel)
  plot.krt <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.krt, CPM.ifn.krt)
  plot.frb <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.frb, CPM.ifn.frb)
  p = plot_grid(plot.mel, plot.krt, plot.frb, labels=c("MEL", "KRT", "FRB"), ncol = 3, nrow = 1)
  p
  return(p)
}
make_PBSeQTL_plot_CPM_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_PBSeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel) + theme(legend.position="none")
  plot.krt <- make_PBSeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt) + theme(legend.position="none")
  plot.frb <- make_PBSeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb) + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, labels=c("MEL", "KRT", "FRB"), ncol = 3, nrow = 1)
  p
  return(p)
}
make_IFNeQTL_plot_CPM_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_IFNeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.ifn.mel,IFN_knownVar.mel) + theme(legend.position="none")
  plot.krt <- make_IFNeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.ifn.krt,IFN_knownVar.krt) + theme(legend.position="none")
  plot.frb <- make_IFNeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.ifn.frb,IFN_knownVar.frb) + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, labels=c("MEL", "KRT", "FRB"), ncol = 3, nrow = 1)
  p
  return(p)
}
make_log2FC_plot <- function(this.snp,this.gene,CPM.pbs,CPM.ifn) {
  
  # pick the snp
  this.snp.genotype <- genotype %>% filter(ID==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(snp.ref,snp.ref), . ) %>%
    gsub("0/1",paste0(snp.ref,snp.alt), . ) %>%
    gsub("1/1",paste0(snp.alt,snp.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype$genotype <- ordered(this.genotype$genotype, 
                                    levels <- c(paste0(snp.ref,snp.ref),
                                                paste0(snp.ref,snp.alt),
                                                paste0(snp.alt,snp.alt)))
  this.genotype = this.genotype %>% na.omit()
  rm(this.snp.genotype)
  
  # pick the gene
  this.CPM.PBS <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
  this.CPM.IFN <- CPM.ifn[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","IFN"))
  this.CPM <- inner_join(this.CPM.PBS, this.CPM.IFN, by="donor") %>% 
    inner_join( . , this.genotype, by="donor")
  this.CPM$log2FC <- log2((this.CPM$IFN + 10) / (this.CPM$PBS + 10))
  
  # pivot to long format and then plot
  this.CPM.long <- this.CPM %>% pivot_longer(.,cols=c(PBS,IFN,log2FC))
  this.CPM.long$name <- ordered(this.CPM.long$name, levels=c("PBS","IFN","log2FC"))
  
  highlightdata1 <- this.CPM.long %>% filter(donor=="F25")
  highlightdata2 <- this.CPM.long %>% filter(donor=="F49")
  highlightdata3 <- this.CPM.long %>% filter(donor=="F55")
  
  p1<-ggplot(this.CPM.long,aes(x=genotype,y=value,color=name)) +
    geom_jitter(size=3,width=0.1,alpha=0.5) +
    scale_color_manual(values=c("#999999","#999999","#999999")) +
    geom_exec(geom_point, data = highlightdata1, group = "donor",
              color = "blue", size=3,position = "identity") +
    geom_exec(geom_point, data = highlightdata2, group = "donor",
              color = "brown",size=3,position = "identity") +
    geom_exec(geom_point, data = highlightdata3, group = "donor",
              color = "black", size=3,position = "identity") +
    facet_wrap(~name,scales="free") +
    theme(legend.position="none") +
    ylab("CPM or log2(CPM+10/CPM+10)")
  return(p1)
}
###### load-data-for-plotting
# load genotype
genotype = data.table::fread(genotype_table,header=TRUE,sep="\t")
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

##### load modeling result
# load  modeling result and name the header
f.mel = data.table::fread(resultF, fill=TRUE) %>% unique() %>% set_rownames(paste0(.$V1, "_", .$V2))
colnames(f.mel)[c(1,2)] <- c("SNP","GENE")
colnames(f.mel)[c(3,4,5,6)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("minimalModel_reQTL_genotype_",x)))
colnames(f.mel)[c(7,8,9,10)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("featureSelected_reQTL_genotype_",x)))
colnames(f.mel)[c(11,12,13,14)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("minimalModel_PBSeQTL_genotype_",x)))
colnames(f.mel)[c(15,16,17,18)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("featureSelected_PBSeQTL_genotype_",x)))
colnames(f.mel)[c(19,20,21,22)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("minimalModel_IFNeQTL_genotype_",x)))
colnames(f.mel)[c(23,24,25,26)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("featureSelected_IFNeQTL_genotype_",x)))
###### filter by effect size
#f.mel.filter1 <- f.mel[which(abs(f.mel$reQTL_genotype_beta) > 0.35),]

################## plot reQTL interesting results ##################
snplst <- c("rs7142855")
this.df <- f.mel[which(f.mel$SNP %in% snplst),]

f.mel.reQTL <- f.mel.filter1 %>% dplyr::filter(reQTL_selectedFeatures_pval < 0.00001) %>% dplyr::arrange(reQTL_selectedFeatures_pval)
#write.table(f.mel.reQTL,paste0(outDir,"/reQTL_significant_results_pvalE-05_KRT.txt"),quote=F,row.names = F)

this.df <- f.mel
plot_list = list()
options(scipen = 0,digits=3) # turn on scientific notation with 3 digits
for (i in 1:nrow(this.df)) {
  this.snp=this.df[i,"SNP"] %>% as.character
  this.gene=this.df[i,"GENE"] %>% as.character
  plot_PBSeqtl_CPM <- make_PBSeQTL_plot_CPM(this.snp, this.gene, CPM.pbs.mel, PBS_knownVar.mel)
  plot_IFNeqtl_CPM <- make_IFNeQTL_plot_CPM(this.snp, this.gene, CPM.ifn.mel, IFN_knownVar.mel)
  plot_reQTL_3cts <- make_reQTL_plot_CPM_3cts(this.snp, this.gene)
  plot_reQTL_log2FC <- make_log2FC_plot(this.snp, this.gene, CPM.pbs.mel, CPM.ifn.mel)
  
  tt <- ttheme_minimal(base_size = 5, core=list(fg_params=list(cex = 2.0)))
  tab=data.frame(metric=c("p.val","beta","stdErr"),
                log2FC=as.character(this.df[i,c("reQTL_genotype_pval","reQTL_genotype_beta","reQTL_genotype_se")]), 
                 PBS_CPM=as.character(this.df[i,c("PBSeQTL_genotype_pval","PBSeQTL_genotype_beta","PBSeQTL_genotype_se")]), 
                 IFN_CPM=as.character(this.df[i,c("IFNeQTL_genotype_pval","IFNeQTL_genotype_beta","IFNeQTL_genotype_se")]))
  plot_table = plot_grid(tableGrob(tab, cols = colnames(tab), rows=NULL, theme = tt), ncol=1, nrow=1)
  
  plotrow1 <- plot_grid(plot_PBSeqtl_CPM, plot_IFNeqtl_CPM, ncol=2, nrow=1)
  plotrow3 <- plot_grid(plot_reQTL_3cts, plot_table, ncol=2,nrow=1)
  p = plot_grid(plotrow1, 
                plot_reQTL_log2FC,
                plotrow3, ncol = 1, nrow = 3)
  plot_list[[i]] = p
}
pdf(paste0(outDir,"/plot/","temp_aa_to_al_pval1E-5.pdf"),width=15,height=8)
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()
options(scipen = 100) # turn back to normal format
cat("done plotting reQTL plots! \n")


#### troubleshoot
rankNorm <- function (y) {
  # input y: CPM across all genes per donor.
  k <- 0.375 # an offset to ensure the z-score is finite. from Blom transform.
  n <- length(y)
  # Ranks.
  r <- rank(y, ties.method = "average")
  # Apply transformation.
  y.rankNorm <- stats::qnorm((r - k) / (n - 2 * k + 1))
  return(y.rankNorm)
}
this.snp="rs7142855"
this.gene="NEDD8"
modelDir="~/Downloads/nl/human/skin/eQTLs/DREG/MEL_minimal"
PBSFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/PBS_",g,".txt")
IFNFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/IFN_",g,".txt")
genotypeFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/genotype.txt")
# pick the snp
this.snp.genotype <- genotype %>% filter(ID==this.snp)
snp.ref=this.snp.genotype$REF
snp.alt=this.snp.genotype$ALT

this.genotype=this.snp.genotype %>% 
  dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
  rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))

this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
this.genotype$genotype = this.genotype$genotype %>%
  gsub("0/0",0, . ) %>%
  gsub("0/1",1, . ) %>%
  gsub("1/1",2, . ) %>%
  gsub("./.",NA, .)
this.genotype$genotype <- ordered(this.genotype$genotype, 
                                  levels <- c(0,1,2))
this.genotype = this.genotype %>% na.omit()
rm(this.snp.genotype)

# pick the gene
this.CPM.PBS <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
  rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
this.CPM.IFN <- CPM.ifn[this.gene,] %>% t() %>% as.data.frame() %>%
  rownames_to_column("donor") %>% set_colnames(c("donor","IFN"))
this.CPM <- inner_join(this.CPM.PBS, this.CPM.IFN, by="donor") %>% 
  inner_join( . , this.genotype, by="donor")
this.CPM$log2FC <- log2((this.CPM$IFN + 10) / (this.CPM$PBS + 10))
this.CPM$rankNormLog2FC <- rankNorm(this.CPM$log2FC)

# 1. build a minimal model with log2FC. does it match the modeling result table?
model.fs <- lm(log2FC ~ genotype, this.CPM)
p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
p.adj <- NA

x <- as.data.frame(summary(model.fs)$coefficients)
x$feature <- rownames(x)

# 2. build a minimal model using rankNorm_log2FC. how does the modeling result change?

# 3. build a model with log2FC and a10b10 covariates. 
combination="a10b10"
covariates_log2FC=data.table::fread(paste0(outDir,"/covariates_phenotype_Log2FCwithDummy10_",combination,".txt"), header=TRUE)

# 4. build a model with rankNorm_log2FC and a10b10 covariates.
combination="a10b10"
covariates_rankNormlog2FC=data.table::fread(paste0(outDir,"/covariates_phenotype_rankNormLog2FCwithDummy10_",combination,".txt"), header=TRUE)
log2fc.dummy = log2((CPM.ifn + 10) / (CPM.pbs + 10))
rankNorm.log2fc <- apply(log2fc.dummy,2,rankNorm) %>% as.data.frame()

##### archive code

# ### 2. PBS eQTL
# f.mel.PBSeQTL <- f.mel[which(f.mel$PBSeQTL_selectedFeatures_pval < 0.01),] %>% dplyr::arrange(PBSeQTL_selectedFeatures_pval)
# write.table(f.mel.PBSeQTL,paste0(outDir,"/PBSeQTL_significant_results_0.01_KRT.txt"),quote=F,row.names = F)
# this.df <- f.mel.PBSeQTL
# plot_list = list()
# options(scipen = 0,digits=3) # turn on scientific notation with 3 digits
# for (i in 1:nrow(this.df)) {
#   this.snp=this.df[i,"SNP"] %>% as.character
#   this.gene=this.df[i,"GENE"] %>% as.character
#   plot_PBSeqtl_CPM <- make_PBSeQTL_plot_CPM(this.snp, this.gene, CPM.pbs.krt, PBS_knownVar.krt)
#   plot_IFNeqtl_CPM <- make_IFNeQTL_plot_CPM(this.snp, this.gene, CPM.ifn.krt, IFN_knownVar.krt)
#   plot_PBSeQTL_3cts = make_PBSeQTL_plot_CPM_3cts(this.snp, this.gene)
# 
#   tt <- ttheme_minimal(base_size = 5, core=list(fg_params=list(cex = 2.0)))
#   tab=data.frame(metric=c("p.val","beta","stdErr"),
#                  log2FC=as.character(this.df[i,c("reQTL_selectedFeatures_pval","reQTL_selectedFeatures_beta","reQTL_selectedFeatures_se")]),
#                  PBS_CPM=as.character(this.df[i,c("PBSeQTL_selectedFeatures_pval","PBSeQTL_selectedFeatures_beta","PBSeQTL_selectedFeatures_se")]),
#                  IFN_CPM=as.character(this.df[i,c("IFNeQTL_selectedFeatures_pval","IFNeQTL_selectedFeatures_beta","IFNeQTL_selectedFeatures_se")]))
#   plot_table = plot_grid(tableGrob(tab, cols = colnames(tab), rows=NULL, theme = tt), ncol=1, nrow=1)
# 
#   p = plot_grid(plot_PBSeqtl_CPM, plot_IFNeqtl_CPM,
#                 plot_PBSeQTL_3cts, plot_table, ncol = 2, nrow = 2)
#   plot_list[[i]] = p
# }
# pdf(paste0(outDir,"/plot/","KRT_PBSeQTL_plots_significant_results_pval0.01.pdf"),width=15,height=8)
# for (i in 1:length(plot_list)) {
#   print(plot_list[[i]])
# }
# dev.off()
# options(scipen = 100) # turn back to normal format
# cat("done plotting PBS eQTL plots! \n")
# 
# ### 3. IFN eQTL
# f.mel.IFNeQTL <- f.mel[which(f.mel$IFNeQTL_selectedFeatures_pval < 0.01),] %>% dplyr::arrange(IFNeQTL_selectedFeatures_pval)
# write.table(f.mel.IFNeQTL,paste0(outDir,"/IFNeQTL_significant_results_0.01_KRT.txt"),quote=F,row.names = F)
# this.df <- f.mel.IFNeQTL
# plot_list = list()
# options(scipen = 0,digits=3) # turn on scientific notation with 3 digits
# for (i in 1:nrow(this.df)) {
#   this.snp=this.df[i,"SNP"] %>% as.character
#   this.gene=this.df[i,"GENE"] %>% as.character
#   plot_PBSeqtl_CPM <- make_PBSeQTL_plot_CPM(this.snp, this.gene, CPM.pbs.krt, PBS_knownVar.krt)
#   plot_IFNeqtl_CPM <- make_IFNeQTL_plot_CPM(this.snp, this.gene, CPM.ifn.krt, IFN_knownVar.krt)
#   plot_IFNeQTL_3cts = make_IFNeQTL_plot_CPM_3cts(this.snp, this.gene)
# 
#   tt <- ttheme_minimal(base_size = 5, core=list(fg_params=list(cex = 2.0)))
#   tab=data.frame(metric=c("p.val","beta","stdErr"),
#                  log2FC=as.character(this.df[i,c("reQTL_selectedFeatures_pval","reQTL_selectedFeatures_beta","reQTL_selectedFeatures_se")]),
#                  PBS_CPM=as.character(this.df[i,c("PBSeQTL_selectedFeatures_pval","PBSeQTL_selectedFeatures_beta","PBSeQTL_selectedFeatures_se")]),
#                  IFN_CPM=as.character(this.df[i,c("IFNeQTL_selectedFeatures_pval","IFNeQTL_selectedFeatures_beta","IFNeQTL_selectedFeatures_se")]))
#   plot_table = plot_grid(tableGrob(tab, cols = colnames(tab), rows=NULL, theme = tt), ncol=1, nrow=1)
# 
#   p = plot_grid(plot_PBSeqtl_CPM, plot_IFNeqtl_CPM,
#                 plot_IFNeQTL_3cts, plot_table, ncol = 2, nrow = 2)
#   plot_list[[i]] = p
# }
# pdf(paste0(outDir,"/plot/","KRT_IFNeQTL_plots_significant_results_pval0.01.pdf"),width=15,height=8)
# for (i in 1:length(plot_list)) {
#   print(plot_list[[i]])
# }
# dev.off()
# options(scipen = 100) # turn back to normal format
# cat("done plotting IFN eQTL plots! \n")
################## plot PBS eQTL interesting results ##################
### 1. PBS eQTL 
snplst <- c("rs9362232")
this.df <- f.mel[which(f.mel$SNP %in% snplst),]

#this.df <- f.mel %>% arrange(GENE,minimalModel_PBSeQTL_genotype_pval)
plot_list = list()
options(scipen = 0,digits=3) # turn on scientific notation with 3 digits
for (i in 1:nrow(this.df)) {
  this.snp=this.df[i,"SNP"] %>% as.character
  this.gene=this.df[i,"GENE"] %>% as.character
  plot_PBSeqtl_CPM <- make_PBSeQTL_plot_CPM(this.snp, this.gene, CPM.pbs.mel, PBS_knownVar.mel)
  plot_IFNeqtl_CPM <- make_IFNeQTL_plot_CPM(this.snp, this.gene, CPM.ifn.mel, IFN_knownVar.mel)
  plot_PBSeQTL_3cts = make_PBSeQTL_plot_CPM_3cts(this.snp, this.gene)
  plot_reQTL_3cts <- make_reQTL_plot_CPM_3cts(this.snp, this.gene)
  
  tt <- ttheme_minimal(base_size = 5, core=list(fg_params=list(cex = 2.0)))
  tab=data.frame(metric=c("p.val","beta","stdErr"),
                 log2FC=as.character(this.df[i,c("minimalModel_reQTL_genotype_pval","minimalModel_reQTL_genotype_beta","minimalModel_reQTL_genotype_se")]), 
                 PBS_CPM=as.character(this.df[i,c("minimalModel_PBSeQTL_genotype_pval","minimalModel_PBSeQTL_genotype_beta","minimalModel_PBSeQTL_genotype_se")]), 
                 IFN_CPM=as.character(this.df[i,c("minimalModel_IFNeQTL_genotype_pval","minimalModel_IFNeQTL_genotype_beta","minimalModel_IFNeQTL_genotype_se")]))
  plot_table = plot_grid(tableGrob(tab, cols = colnames(tab), rows=NULL, theme = tt), ncol=1, nrow=1)
  
  plotrow1 <- plot_grid(plot_PBSeqtl_CPM, plot_IFNeqtl_CPM, ncol=2, nrow=1)
  plotrow3 <- plot_grid(plot_reQTL_3cts, plot_table, ncol=2,nrow=1)
  p = plot_grid(plotrow1, 
                plot_PBSeQTL_3cts,
                plotrow3, ncol = 1, nrow = 3)
  plot_list[[i]] = p
}
pdf(paste0(outDir,"/plot/","rs9362232.pdf"),width=15,height=8)
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()
options(scipen = 100) # turn back to normal format
cat("done plotting PBSeQTL plots! \n")

################## COMPILE TABLES FOR GENE OF INTEREST OR IN SHINY FOR SNP OF INTEREST ################################
enhancerF <- "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method4/promoters_and_enhancers_surrounding_genes_all3cts.bed"
enhancerATACF <- "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ATACseq_allcts_allregions_1kbp.bed"
enhancerK27acF = "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ChIPseq_allcts_allregions_1kbp.bed"

enhancerDict <- read.table(enhancerF,header=F,sep="\t") %>% 
  set_colnames(c("chr","start","end","name","ID","gene","celltype","regionType")) %>%
  dplyr::filter(gene %in% unique(f.mel$GENE))

count_ATAC <- data.table::fread(enhancerATACF) %>%
  dplyr::select(-c(chr,start,end)) %>% column_to_rownames("name")
col.sum=apply(count_ATAC,2,sum); 
CPM_atac=sweep(count_ATAC, 2, col.sum, `/`)*1000000
CPM_atac <- CPM_atac[,c("F25_KRT_PBS","F49_KRT_PBS","F55_KRT_PBS",
                        "F25_MEL_PBS","F49_MEL_PBS","F55_MEL_PBS",
                        "F25_FRB_PBS","F49_FRB_PBS","F55_FRB_PBS",
                        "F25_KRT_IFN","F49_KRT_IFN","F55_KRT_IFN",
                        "F25_MEL_IFN","F49_MEL_IFN","F55_MEL_IFN",
                        "F25_FRB_IFN","F49_FRB_IFN","F55_FRB_IFN")]

enhancerATAC <-  CPM_atac %>% rownames_to_column("name") %>% dplyr::filter(name %in% enhancerDict$ID)
rm(count_ATAC)

count_K27ac <- data.table::fread(enhancerK27acF) %>%
  dplyr::select(-c(chr,start,end)) %>% column_to_rownames("name")
col.sum=apply(count_K27ac,2,sum); 
CPM_K27ac=sweep(count_K27ac, 2, col.sum, `/`)*1000000
CPM_K27ac <- CPM_K27ac[,c("F25_KRT_PBS","F49_KRT_PBS","F55_KRT_PBS",
                          "F25_MEL_PBS","F49_MEL_PBS","F55_MEL_PBS",
                          "F25_FRB_PBS","F49_FRB_PBS","F55_FRB_PBS",
                          "F25_KRT_IFN","F49_KRT_IFN","F55_KRT_IFN",
                          "F25_MEL_IFN","F49_MEL_IFN","F55_MEL_IFN",
                          "F25_FRB_IFN","F49_FRB_IFN","F55_FRB_IFN")]

enhancerK27ac <-  CPM_K27ac %>% rownames_to_column("name") %>% dplyr::filter(name %in% enhancerDict$ID)
rm(count_K27ac)

######### compile SNP-gene-enhancer table
compile_SnpGeneEnhancerTable <- function(this.snp, this.gene){
  SnpGeneEnhancerTable <- dplyr::right_join(f.mel[,c("SNP","GENE")], 
                                            enhancerDict[,c("ID","gene","regionType","start","end")], 
                                            by=c("GENE"="gene"), relationship = "many-to-many") %>% unique()
  SnpGeneEnhancerTable$enhancer_midpoint <- SnpGeneEnhancerTable$start + 1 + 150
  
  this.snp.position <- as.numeric(genotype[which(genotype$ID==this.snp),"END"])
  this.snp.REF <- as.character(genotype[which(genotype$ID==this.snp),"REF"])
  this.snp.ALT <- as.character(genotype[which(genotype$ID==this.snp),"ALT"])
  this.snp.genotype.numeric <- genotype[which(genotype$ID==this.snp),c("F25:GT","F49:GT","F55:GT")] %>%
    as.character() %>%
    gsub("0/0",0, . ) %>%
    gsub("0/1",1, . ) %>%
    gsub("1/1",2, . ) %>%
    gsub("./.",NA, . ) %>%
    as.numeric()
  this.snp.genotype.char <- genotype[which(genotype$ID==this.snp),c("F25:GT","F49:GT","F55:GT")] %>%
    as.character() %>%
    gsub("0/0",paste0(this.snp.REF,this.snp.REF), . ) %>%
    gsub("0/1",paste0(this.snp.REF,this.snp.ALT), . ) %>%
    gsub("1/1",paste0(this.snp.ALT,this.snp.ALT), . ) %>%
    gsub("./.",NA, . ) %>%
    as.character() %>%
    paste( . , collapse="_")
  
  this_SnpGeneEnhancerTable <- SnpGeneEnhancerTable %>% dplyr::filter(GENE==this.gene & SNP==this.snp) 
  this_SnpGeneEnhancerTable$SNP_position <- this.snp.position
  this_SnpGeneEnhancerTable$dist <- this_SnpGeneEnhancerTable$enhancer_midpoint - this_SnpGeneEnhancerTable$SNP_position
  this_SnpGeneEnhancerTable <- this_SnpGeneEnhancerTable %>% arrange(dist)
  
  # make a table of this snp's nearby enhancers of all ATACseq and H3K27ac ChIPseq signals
  result.table = data.frame()
  for (i in 1:nrow(this_SnpGeneEnhancerTable)) {
    this.region <- as.character(this_SnpGeneEnhancerTable[i,"ID"])
    
    this.region.celltype <- enhancerDict[which(enhancerDict$ID==this.region),"celltype"] %>%
      paste(. , collapse=",")
    
    this.region.dist <- as.numeric(this_SnpGeneEnhancerTable[i,"dist"])
    
    this.ATAC.PBS <- enhancerATAC %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_PBS")) %>%
      as.numeric()
    
    this.ATAC.IFN <- enhancerATAC %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_IFN")) %>%
      as.numeric()
    
    this.K27ac.PBS <- enhancerK27ac %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_PBS")) %>%
      as.numeric()
    
    this.K27ac.IFN <- enhancerK27ac %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_IFN")) %>%
      as.numeric()
    
    z.ATAC = as.numeric(scale(c(this.ATAC.PBS, this.ATAC.IFN)))
    z.diff.ATAC = z.ATAC[4:6] - z.ATAC[1:3]
    
    z.K27ac = as.numeric(scale(c(this.K27ac.PBS, this.K27ac.IFN)))
    z.diff.K27ac = z.K27ac[4:6] - z.K27ac[1:3]
    
    cor.ATAC <- cor(as.numeric(scale(this.snp.genotype.numeric)), z.diff.ATAC, method="pearson")
    cor.K27ac <- cor(as.numeric(scale(this.snp.genotype.numeric)), z.diff.K27ac, method="pearson")
    
    this.result <- c(this.snp, this.region, this.region.celltype, this.region.dist,
                     round(cor.ATAC,3), round(cor.K27ac,3),
                     paste(round(this.ATAC.PBS,2),collapse="_"), 
                     paste(round(this.ATAC.IFN,2),collapse="_"),
                     paste(round(this.K27ac.PBS,2),collapse="_"),
                     paste(round(this.K27ac.IFN,2),collapse="_"))
    
    result.table = rbind(result.table, this.result)
  }
  colnames(result.table) <- c("SNP","regionID","celltype","dist_to_SNP",
                              "PCC_Genotype_ATAC_zDiff","PCC_Genotype_K27ac_zDiff",
                              "ATAC_PBS_3donors_CPM", "ATAC_IFN_3donors_CPM",
                              "K27ac_PBS_3donors_CPM", "K27ac_IFN_3donors_CPM") 
  
  closest.idx <- which(result.table$dist_to_SNP == min(abs(as.numeric(result.table$dist_to_SNP))))
  correlatedATAC.idx <- which(abs(as.numeric(result.table$PCC_Genotype_ATAC_zDiff)) > 0.8)
  correlatedK27ac.idx <- which(abs(as.numeric(result.table$PCC_Genotype_K27ac_zDiff)) > 0.8)
  most.correlatedATAC.idx <- which(abs(as.numeric(result.table$PCC_Genotype_ATAC_zDiff)) == 
                                     max(abs(as.numeric(result.table$PCC_Genotype_ATAC_zDiff))))
  most.correlatedK27ac.idx <- which(abs(as.numeric(result.table$PCC_Genotype_K27ac_zDiff)) == 
                                      max(abs(as.numeric(result.table$PCC_Genotype_K27ac_zDiff))))
  result.table$tag1 <- ""
  result.table$tag1[closest.idx] <- "\u2021"
  result.table$tag2 <- ""
  result.table$tag2[setdiff(correlatedATAC.idx,most.correlatedATAC.idx)] <- "\u25B3" # empty square
  result.table$tag2[most.correlatedATAC.idx] <- "\u25B2" # filled square
  result.table$tag3 <- ""
  result.table$tag3[setdiff(correlatedK27ac.idx, most.correlatedK27ac.idx)] <- "\u25EF" # empty circle
  result.table$tag3[most.correlatedK27ac.idx] <- "\u2B24" # filled circle
  result.table$tag <- paste0(result.table$tag1, result.table$tag2, result.table$tag3)
  result.table <- result.table %>% dplyr::select(-c("tag1","tag2","tag3"))
  result.table$dist_to_SNP <- prettyNum(result.table$dist_to_SNP, big.mark=",")
  result.table <- result.table[,c(1,11,2:10)]
  
  result.table$regionID <- gsub("merged_all_skin-eQTL_ATACseq_files_","",result.table$regionID)
  return(result.table)
}
compile_SnpGeneEnhancerTable_PBS <- function(this.snp, this.gene){
  SnpGeneEnhancerTable <- dplyr::right_join(f.mel[,c("SNP","GENE")], 
                                            enhancerDict[,c("ID","gene","regionType","start","end")], 
                                            by=c("GENE"="gene"), relationship = "many-to-many") %>% unique()
  SnpGeneEnhancerTable$enhancer_midpoint <- SnpGeneEnhancerTable$start + 1 + 150
  
  this.snp.position <- as.numeric(genotype[which(genotype$ID==this.snp),"END"])
  this.snp.REF <- as.character(genotype[which(genotype$ID==this.snp),"REF"])
  this.snp.ALT <- as.character(genotype[which(genotype$ID==this.snp),"ALT"])
  this.snp.genotype.numeric <- genotype[which(genotype$ID==this.snp),c("F25:GT","F49:GT","F55:GT")] %>%
    as.character() %>%
    gsub("0/0",0, . ) %>%
    gsub("0/1",1, . ) %>%
    gsub("1/1",2, . ) %>%
    gsub("./.",NA, . ) %>%
    as.numeric()
  this.snp.genotype.char <- genotype[which(genotype$ID==this.snp),c("F25:GT","F49:GT","F55:GT")] %>%
    as.character() %>%
    gsub("0/0",paste0(this.snp.REF,this.snp.REF), . ) %>%
    gsub("0/1",paste0(this.snp.REF,this.snp.ALT), . ) %>%
    gsub("1/1",paste0(this.snp.ALT,this.snp.ALT), . ) %>%
    gsub("./.",NA, . ) %>%
    as.character() %>%
    paste( . , collapse="_")
  
  this_SnpGeneEnhancerTable <- SnpGeneEnhancerTable %>% dplyr::filter(GENE==this.gene & SNP==this.snp) 
  this_SnpGeneEnhancerTable$SNP_position <- this.snp.position
  this_SnpGeneEnhancerTable$dist <- this_SnpGeneEnhancerTable$enhancer_midpoint - this_SnpGeneEnhancerTable$SNP_position
  this_SnpGeneEnhancerTable <- this_SnpGeneEnhancerTable %>% arrange(dist)
  
  # make a table of this snp's nearby enhancers of all ATACseq and H3K27ac ChIPseq signals
  result.table = data.frame()
  for (i in 1:nrow(this_SnpGeneEnhancerTable)) {
    this.region <- as.character(this_SnpGeneEnhancerTable[i,"ID"])
    
    this.region.celltype <- enhancerDict[which(enhancerDict$ID==this.region),"celltype"] %>%
      paste(. , collapse=",")
    
    this.region.dist <- as.numeric(this_SnpGeneEnhancerTable[i,"dist"])
    
    this.ATAC.PBS <- enhancerATAC %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_PBS")) %>%
      as.numeric()
    
    this.ATAC.IFN <- enhancerATAC %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_IFN")) %>%
      as.numeric()
    
    this.K27ac.PBS <- enhancerK27ac %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_PBS")) %>%
      as.numeric()
    
    this.K27ac.IFN <- enhancerK27ac %>% dplyr::filter(name==this.region) %>% 
      dplyr::select(ends_with("KRT_IFN")) %>%
      as.numeric()
    
    z.ATAC = as.numeric(scale(c(this.ATAC.PBS)))
    z.diff.ATAC = z.ATAC[1:3]
    
    z.K27ac = as.numeric(scale(c(this.K27ac.PBS,)))
    z.diff.K27ac = z.K27ac[1:3]
    
    cor.ATAC <- cor(as.numeric(scale(this.snp.genotype.numeric)), z.diff.ATAC, method="pearson")
    cor.K27ac <- cor(as.numeric(scale(this.snp.genotype.numeric)), z.diff.K27ac, method="pearson")
    
    this.result <- c(this.snp, this.region, this.region.celltype, this.region.dist,
                     round(cor.ATAC,3), round(cor.K27ac,3),
                     paste(round(this.ATAC.PBS,2),collapse="_"), 
                     paste(round(this.ATAC.IFN,2),collapse="_"),
                     paste(round(this.K27ac.PBS,2),collapse="_"),
                     paste(round(this.K27ac.IFN,2),collapse="_"))
    
    result.table = rbind(result.table, this.result)
  }
  colnames(result.table) <- c("SNP","regionID","celltype","dist_to_SNP",
                              "PCC_Genotype_ATAC_zPBS","PCC_Genotype_K27ac_zPBS",
                              "ATAC_PBS_3donors_CPM", "ATAC_IFN_3donors_CPM",
                              "K27ac_PBS_3donors_CPM", "K27ac_IFN_3donors_CPM") 
  
  closest.idx <- which(result.table$dist_to_SNP == min(abs(as.numeric(result.table$dist_to_SNP))))
  correlatedATAC.idx <- which(abs(as.numeric(result.table$PCC_Genotype_ATAC_zDiff)) > 0.8)
  correlatedK27ac.idx <- which(abs(as.numeric(result.table$PCC_Genotype_K27ac_zDiff)) > 0.8)
  most.correlatedATAC.idx <- which(abs(as.numeric(result.table$PCC_Genotype_ATAC_zDiff)) == 
                                     max(abs(as.numeric(result.table$PCC_Genotype_ATAC_zDiff))))
  most.correlatedK27ac.idx <- which(abs(as.numeric(result.table$PCC_Genotype_K27ac_zDiff)) == 
                                      max(abs(as.numeric(result.table$PCC_Genotype_K27ac_zDiff))))
  result.table$tag1 <- ""
  result.table$tag1[closest.idx] <- "\u2021"
  result.table$tag2 <- ""
  result.table$tag2[setdiff(correlatedATAC.idx,most.correlatedATAC.idx)] <- "\u25B3" # empty square
  result.table$tag2[most.correlatedATAC.idx] <- "\u25B2" # filled square
  result.table$tag3 <- ""
  result.table$tag3[setdiff(correlatedK27ac.idx, most.correlatedK27ac.idx)] <- "\u25EF" # empty circle
  result.table$tag3[most.correlatedK27ac.idx] <- "\u2B24" # filled circle
  result.table$tag <- paste0(result.table$tag1, result.table$tag2, result.table$tag3)
  result.table <- result.table %>% dplyr::select(-c("tag1","tag2","tag3"))
  result.table$dist_to_SNP <- prettyNum(result.table$dist_to_SNP, big.mark=",")
  result.table <- result.table[,c(1,11,2:10)]
  
  result.table$regionID <- gsub("merged_all_skin-eQTL_ATACseq_files_","",result.table$regionID)
  return(result.table)
}
######### compile SNP-trait table
compile_LDtraitTable <- function(this.snp) {
  #library(LDlinkR)
  tryCatch({
    my_token = "908c3efbf915"
    this.res <- LDtrait(snps = this.snp,
                        pop = "CEU",
                        r2d = "r2",
                        r2d_threshold = 0.6,
                        win_size = 5e+05,
                        token = my_token,
                        genome_build = "grch38")
    this.res <- this.res %>% arrange(desc(R2)) 
    return(this.res)
  }, 
  error = function(err) {
    cat(paste0(this.snp,": no associated trait was found\n"))
    this.res <- as.data.frame(t(c(this.snp,rep("NA",11))))
    return(this.res)
  })
}
######### compile SNP-HWE table
compile_HWETable <- function(this.snp) {
  
  this.snp.REF <- genotype[which(genotype$ID==this.snp),"REF"]
  this.snp.ALT <- genotype[which(genotype$ID==this.snp),"ALT"]
  this.snp.genotype.vector <- genotype[which(genotype$ID==this.snp),] %>% 
    dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>%
    as.character() %>%
    gsub("0/0",0, . ) %>%
    gsub("0/1",1, . ) %>%
    gsub("1/1",2, . ) %>%
    gsub("./.",NA, . )
  this.snp.genotype.summary <- c(MM=length(which(this.snp.genotype.vector==0)),
                                 MN=length(which(this.snp.genotype.vector==1)),
                                 NN=length(which(this.snp.genotype.vector==2))) 
  
  HW.test.pval <- HardyWeinberg::HWChisq(this.snp.genotype.summary,verbose=TRUE)$pval
  HW.test.verdict <- ifelse(HW.test.pval<0.05, "not in equilibrium", "in HW equilibrium")
  HWE.summary <- data.frame(REF=this.snp.REF, 
                               ALT=this.snp.ALT, 
                               counts=paste(this.snp.genotype.summary,collapse="_"),
                               HW_Chi2_test_pval=HW.test.pval, 
                               conclusion=HW.test.verdict)
  rownames(HWE.summary) <- this.snp
  return(HWE.summary)
}
####### GENE-centric approach: compile table for all SNPs of a gene ######### 
this.gene <- "SYNE2"
this.gene.snps <- f.mel %>% dplyr::filter(GENE==this.gene) %>% pull(SNP)

this.gene.SnpGeneEnhancerTable <- data.frame()
for (i in 1:length(this.gene.snps)) {
  this.result <- compile_SnpGeneEnhancerTable(this.gene.snps[i], this.gene)
  this.gene.SnpGeneEnhancerTable <- rbind(this.gene.SnpGeneEnhancerTable, this.result)
}

this.gene.LDtraitTable <- data.frame()
for (i in 1:length(this.gene.snps)) {
  this.result <- compile_LDtraitTable(this.gene.snps[i])
  this.gene.LDtraitTable <- rbind(this.gene.LDtraitTable, this.result)
}

this.HWETable <- data.frame()
for (i in 1:length(this.gene.snps)) {
  print(i)
  this.result <- compile_HWETable(this.gene.snps[i])
  this.HWETable <- rbind(this.HWETable, this.result)
}

View(this.gene.LDtraitTable)
View(this.gene.SnpGeneEnhancerTable)
View(this.HWETable)

####### global approach: compile LDtrait table for top eQTL for each gene ####### 
all_genes <- unique(f.mel$GENE)
LDtraitTable <- data.frame()
for (i in 1:length(all_genes)) {
  this.gene <- all_genes[i]
  this.top.snp <- f.mel %>% dplyr::filter(GENE==this.gene) %>% top_n(., -1, minimalModel_PBSeQTL_genotype_pval) %>% pull(SNP)
  this.result <- compile_LDtraitTable(this.top.snp[1])
  this.result$Gene <- this.gene
  colnames(this.result) <- colnames(LDtraitTable)
  LDtraitTable <- rbind(LDtraitTable, this.result)
}
colnames(LDtraitTable) <- c("Query","GWAS_Trait","PMID","RS_Number","Position_GRCh38","Alleles","R2","D'","Risk_Allele","Effect_Size_95_CI","Beta_or_OR","P_value","Gene")
View(LDtraitTable)
LDtraitTable.filtered <- LDtraitTable %>% dplyr::filter(PMID!="NA")
sort(table(LDtraitTable.filtered$GWAS_Trait))
View(LDtraitTable.filtered)
####### put every table in shiny ####### 
library(shiny)
library(DT) # for rendering data.tables
library(HardyWeinberg)
library(LDlinkR)
library(dplyr)
library(htmltools)
this.df <- f.mel %>% arrange(GENE,reQTL_genotype_pval)
for (i in 1:nrow(this.df[1:5,])) {
  this.snp=this.df[i,"SNP"] %>% as.character
  this.gene=this.df[i,"GENE"] %>% as.character
  # data.frames
  data_frame1 <- compile_SnpGeneEnhancerTable(this.snp, this.gene)
  data_frame2 <- compile_LDtraitTable(this.snp)
  data_frame3 <- compile_HWETable(this.snp)
  
  result.table.title <- paste0(this.gene," : ",this.snp," [", this.snp.genotype.char,"] ")
  result.table.subtitle <- paste0("\u2021"," : enhancer midpoint closest to SNP\n",
                                  "\u25B3"," : ATACseq signal correlates with SNP genotype (PCC>0.8)\n",
                                  "\u25B2"," : most correlated ATACseq signal\n",
                                  "\u25EF"," : K27ac signal correlates with SNP genotype (PCC>0.8)\n",
                                  "\u2B24"," : most correlated K27ac signal")
  ui <- fluidPage(
    titlePanel(result.table.title),
    mainPanel(
      tabsetPanel(
        tabPanel("Enhancers", DTOutput("table1"), p(result.table.subtitle)),
        tabPanel("LDtrait", DTOutput("table2"), p("LDtrait results, population = CEU, LD R2 cutoff = 0.6")),
        tabPanel("HWE", DTOutput("table3"), p("Hardy Weinberg Equilibrium Chi-square test results"))
      )
    )
  )
  
  server <- function(input, output, session) {
    output$table1 <- renderDT({
      datatable(data_frame1, options = list(scrollY = "400px", scrollX="1200px"))
    })
    output$table2 <- renderDT({
      datatable(data_frame2, options = list(scrollY = "400px", scrollX="1200px"))
    })
    output$table3 <- renderDT({
      datatable(data_frame3, options = list(scrollY = "400px"))
    })
  }
  
  app <- shinyApp(ui = ui, server = server)
  # Save the Shiny app as an HTML file
  htmltools::save_html(app, paste0(outDir,"/plot/",this.snp,"_",this.gene,"_information",".html"))
  stopApp()a
}

################## prioritize SNP - gene candidates ##############
###### filter by overlapping expanded enhancers
# filter by abs(PCC) > 0.9
result.table.filtered <- result.table %>% na.omit() %>%
  dplyr::filter(PCC_genotype_ATACFC >= 0.9 | PCC_genotype_K27acFC >= 0.9 | 
                  PCC_genotype_ATACFC <= -0.9 | PCC_genotype_K27acFC <= -0.9) %>%
  left_join(., genotype[,c("CHROM","START","ID","REF","ALT")], by=c("snp"="ID")) %>%
  dplyr::filter(PCC_genotype_K27acFC!="NaN") %>%
  dplyr::filter(regionID %in% enhancerDict[which(enhancerDict$celltype=="KRT"),"ID"])

###### filter by overlapping skin disease GWAS SNPs
snpTagF <- "~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/skinDis_melanoma_autoImmuneDis_snps_and_LD0.8plus_tags.bed"
snpTag <- data.table::fread(snpTagF) %>% set_colnames(c("chr","start","end","snp","R","anchorsnp"))
snp_list <- unique(snpTag$snp)

result.table.filtered2 <- result.table.filtered %>% dplyr::filter(snp %in% snp_list)
snpTag.filtered <- snpTag[which(snpTag$snp %in% unique(f.mel$SNP)),] %>% dplyr::arrange(desc(R))




###### filter by overlapping expanded enhancers ####
enhancerF <- "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method4/promoters_and_enhancers_surrounding_genes_all3cts.bed"
enhancerATACF <- "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ATACseq_allcts_allregions_1kbp.bed"
enhancerK27acF = "~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ChIPseq_allcts_allregions_1kbp.bed"

enhancerDict <- read.table(enhancerF,header=F,sep="\t") %>% 
  set_colnames(c("chr","start","end","name","ID","gene","celltype","regionType")) %>%
  dplyr::filter(gene %in% unique(f.mel$GENE))

count_ATAC <- data.table::fread(enhancerATACF) %>%
  dplyr::select(-c(chr,start,end)) %>% column_to_rownames("name")
col.sum=apply(count_ATAC,2,sum); 
CPM_atac=sweep(count_ATAC, 2, col.sum, `/`)*1000000
CPM_atac <- CPM_atac[,c("F25_KRT_PBS","F49_KRT_PBS","F55_KRT_PBS",
                        "F25_MEL_PBS","F49_MEL_PBS","F55_MEL_PBS",
                        "F25_FRB_PBS","F49_FRB_PBS","F55_FRB_PBS",
                        "F25_KRT_IFN","F49_KRT_IFN","F55_KRT_IFN",
                        "F25_MEL_IFN","F49_MEL_IFN","F55_MEL_IFN",
                        "F25_FRB_IFN","F49_FRB_IFN","F55_FRB_IFN")]

enhancerATAC <-  CPM_atac %>% rownames_to_column("name") %>% dplyr::filter(name %in% enhancerDict$ID)
rm(count_ATAC)

count_K27ac <- data.table::fread(enhancerK27acF) %>%
  dplyr::select(-c(chr,start,end)) %>% column_to_rownames("name")
col.sum=apply(count_K27ac,2,sum); 
CPM_K27ac=sweep(count_K27ac, 2, col.sum, `/`)*1000000
CPM_K27ac <- CPM_K27ac[,c("F25_KRT_PBS","F49_KRT_PBS","F55_KRT_PBS",
                          "F25_MEL_PBS","F49_MEL_PBS","F55_MEL_PBS",
                          "F25_FRB_PBS","F49_FRB_PBS","F55_FRB_PBS",
                          "F25_KRT_IFN","F49_KRT_IFN","F55_KRT_IFN",
                          "F25_MEL_IFN","F49_MEL_IFN","F55_MEL_IFN",
                          "F25_FRB_IFN","F49_FRB_IFN","F55_FRB_IFN")]

enhancerK27ac <-  CPM_K27ac %>% rownames_to_column("name") %>% dplyr::filter(name %in% enhancerDict$ID)
rm(count_K27ac)

# compile SNP-gene-enhancer table
SnpGeneEnhancerTable <- right_join(f.mel[,c("SNP","GENE")], 
                                   enhancerDict[,c("ID","gene","regionType")], 
                                   by=c("GENE"="gene")) %>% unique()
result.table = data.frame()
for (i in 1:nrow(SnpGeneEnhancerTable)) {
  this.snp <- as.character(SnpGeneEnhancerTable[i,"SNP"])
  this.gene <- as.character(SnpGeneEnhancerTable[i,"GENE"])
  this.region <- as.character(SnpGeneEnhancerTable[i,"ID"])
  
  #this.snp <- "rs2910686"
  #this.gene <- "ERAP2"
  #this.region <- "merged_all_skin-eQTL_ATACseq_files_peak_313795"
  
  this.genotype <- genotype[which(genotype$ID==this.snp),c("F25:GT","F49:GT","F55:GT")] %>%
    as.character() %>%
    gsub("0/0",0, . ) %>%
    gsub("0/1",1, . ) %>%
    gsub("1/1",2, . ) %>%
    gsub("./.",NA, . ) %>%
    as.numeric()
  
  this.ATAC.PBS <- enhancerATAC %>% dplyr::filter(name==this.region) %>% 
    dplyr::select(ends_with("KRT_PBS")) %>%
    as.numeric()
  
  this.ATAC.IFN <- enhancerATAC %>% dplyr::filter(name==this.region) %>% 
    dplyr::select(ends_with("KRT_IFN")) %>%
    as.numeric()
  
  this.K27ac.PBS <- enhancerK27ac %>% dplyr::filter(name==this.region) %>% 
    dplyr::select(ends_with("KRT_PBS")) %>%
    as.numeric()
  
  this.K27ac.IFN <- enhancerK27ac %>% dplyr::filter(name==this.region) %>% 
    dplyr::select(ends_with("KRT_IFN")) %>%
    as.numeric()
  
  this.df <- data.frame(genotype=this.genotype, ATAC_PBS=this.ATAC.PBS, ATAC_IFN=this.ATAC.IFN, K27ac_PBS=this.K27ac.PBS, K27ac_IFN=this.K27ac.IFN)
  this.df$ATAC_log2FC <- log2( this.df$ATAC_IFN / this.df$ATAC_PBS )
  this.df$K27ac_log2FC <- log2( this.df$K27ac_IFN / this.df$K27ac_PBS )
  
  cor1 <- cor(this.df$genotype, this.df$ATAC_log2FC, method="pearson")
  cor2 <- cor(this.df$genotype, this.df$K27ac_log2FC, method="pearson")
  
  this.result <- c(this.snp, this.gene, this.genotype, this.region, cor1, cor2)
  result.table = rbind(result.table, this.result)
}
colnames(result.table) <- c("snp","gene","NA1","NA2","NA3","regionID","PCC_genotype_ATACFC","PCC_genotype_K27acFC") 
result.table <- result.table %>%  dplyr::select(-c("NA1","NA2","NA3"))
rm(CPM_atac, CPM_K27ac)

# filter by abs(PCC) > 0.9
result.table.filtered <- result.table %>% na.omit() %>%
  dplyr::filter(PCC_genotype_ATACFC >= 0.9 | PCC_genotype_K27acFC >= 0.9 | 
                  PCC_genotype_ATACFC <= -0.9 | PCC_genotype_K27acFC <= -0.9) %>%
  left_join(., genotype[,c("CHROM","START","ID","REF","ALT")], by=c("snp"="ID"))

###### filter by overlapping skin disease GWAS SNPs ######
snpTagF <- "~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/skinDis_melanoma_autoImmuneDis_snps_and_LD0.8plus_tags.bed"
snpTag <- data.table::fread(snpTagF) %>% set_colnames(c("chr","start","end","snp","R","anchorsnp"))
snp_list <- unique(snpTag$snp)

result.table.filtered2 <- result.table.filtered %>% dplyr::filter(snp %in% snp_list)
snpTag.filtered <- snpTag[which(snpTag$snp %in% unique(f.mel$SNP)),] %>% dplyr::arrange(desc(R))


