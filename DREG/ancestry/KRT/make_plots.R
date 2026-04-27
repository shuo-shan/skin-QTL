#!/usr/bin/env Rscript
# written by Crystal Shan 09/2023

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
#library(LDlinkR)
#library(HardyWeinberg)

# don't forget to replace "~/Downloads/nl/" with "~/Downloads/nl" while toggling between cluster job submission and local interactive Rstudio
# parsing input arguments
args = commandArgs(trailingOnly=TRUE)
outDir=args[1]
gene=args[2]
snp=args[3]
# 
# gene="DHX58"
# snp="rs739636"
# outDir=paste0("~/Downloads/nl/human/skin/eQTLs/website/data/plots","/",gene,"/",snp)

# load genotype
genotype_table="~/Downloads/nl/human/skin/eQTLs/DREG/ancestry/KRT/selected_snps_genotype.bed"
genotype = data.table::fread(genotype_table,sep="\t")


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
# functions for ancestry colored scatter plots in CPM
make_PBSeQTL_plot_CPM_ancestryColored <- function(this.snp,this.gene,CPM.pbs,PBS_knownVar,celltype) {
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
    ggtitle(paste0("PBS","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM")) +
    scale_color_brewer(palette="Dark2")
  p1
  return(p1)
}
make_IFNeQTL_plot_CPM_ancestryColored <- function(this.snp,this.gene,CPM.ifn,IFN_knownVar, celltype) {
  
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
    ggtitle(paste0("IFN","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM")) +
    scale_color_brewer(palette="Dark2")
  p1
  return(p1)
}
make_PBSeQTL_plot_CPM_ancestryColored_3cts <- function(this.snp, this.gene){
  plot.mel <- make_PBSeQTL_plot_CPM_ancestryColored(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel,"MEL")
  plot.krt <- make_PBSeQTL_plot_CPM_ancestryColored(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt,"KRT")
  plot.frb <- make_PBSeQTL_plot_CPM_ancestryColored(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb,"FRB")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
  
}
make_IFNeQTL_plot_CPM_ancestryColored_3cts <- function(this.snp, this.gene){
  plot.mel <- make_IFNeQTL_plot_CPM_ancestryColored(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel,"MEL")
  plot.krt <- make_IFNeQTL_plot_CPM_ancestryColored(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt,"KRT")
  plot.frb <- make_IFNeQTL_plot_CPM_ancestryColored(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb,"FRB")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}

# functions for ancestry colored scatter plots in rank normalized CPM
make_PBSeQTL_plot_rankNormCPM_ancestryColored <- function(this.snp,this.gene,CPM.pbs,PBS_knownVar,celltype) {
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
  this.CPM.PBS <- CPM.pbs[this.gene,] %>% rankNorm() %>% as.data.frame() %>%
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
    ggtitle(paste0("PBS","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," rank normalized value")) +
    scale_color_brewer(palette="Dark2")
  p1
  return(p1)
}
make_IFNeQTL_plot_rankNormCPM_ancestryColored <- function(this.snp,this.gene,CPM.ifn,IFN_knownVar, celltype) {
  
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
  this.CPM.IFN <- CPM.ifn[this.gene,] %>% rankNorm() %>% as.data.frame() %>%
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
    ggtitle(paste0("IFN","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," rank normalized value")) +
    scale_color_brewer(palette="Dark2")
  p1
  return(p1)
}
make_PBSeQTL_plot_rankNormCPM_ancestryColored_3cts <- function(this.snp, this.gene){
  plot.mel <- make_PBSeQTL_plot_rankNormCPM_ancestryColored(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel,"MEL")
  plot.krt <- make_PBSeQTL_plot_rankNormCPM_ancestryColored(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt,"KRT")
  plot.frb <- make_PBSeQTL_plot_rankNormCPM_ancestryColored(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb,"FRB")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}
make_IFNeQTL_plot_rankNormCPM_ancestryColored_3cts <- function(this.snp, this.gene){
  plot.mel <- make_IFNeQTL_plot_rankNormCPM_ancestryColored(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel,"MEL")
  plot.krt <- make_IFNeQTL_plot_rankNormCPM_ancestryColored(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt,"KRT")
  plot.frb <- make_IFNeQTL_plot_rankNormCPM_ancestryColored(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb,"FRB")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}

# functions for 3 donor colored scatter plots in CPM
make_PBSeQTL_plot_CPM_3donorColored <- function(this.snp,this.gene,CPM.pbs,PBS_knownVar, celltype) {
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
    ggtitle(paste0("PBS","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM"))
  p1
  return(p1)
}
make_IFNeQTL_plot_CPM_3donorColored <- function(this.snp,this.gene,CPM.ifn,IFN_knownVar, celltype) {
  
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
  this.CPM.IFN <- CPM.ifn[this.gene,] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","IFN")) %>% 
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
  p1<-ggplot(this.CPM.IFN, aes(x=genotype, y=IFN)) +
    geom_exec(geom_point, data = highlightdata0, group = "donor",
              color = "#5A5A5A", size = 3) +
    geom_exec(geom_point, data = highlightdata1, group = "donor",
              color = "blue", size = 4) +
    geom_exec(geom_point, data = highlightdata2, group = "donor",
              color = "brown", size = 4) +
    geom_exec(geom_point, data = highlightdata3, group = "donor",
              color = "black", size = 4) +
    ggtitle(paste0("IFN","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," CPM"))
  p1
  return(p1)
}
make_log2FC_plot_3donorColored <- function(this.snp,this.gene,CPM.pbs,CPM.ifn, celltype) {
  
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
  this.CPM <- this.CPM %>% dplyr::select(-c("PBS","IFN"))
  
  # pivot to long format and then plot
  this.CPM.long <- this.CPM %>% pivot_longer(.,cols=c(log2FC))
  
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
    theme(legend.position="none") +
    ggtitle(paste0("log2FC","\n",celltype,"\n",this.snp,":",this.gene)) +
    ylab("log2(CPM+10/CPM+10)")
  return(p1)
}
make_PBSeQTL_plot_CPM_3donorColored_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_PBSeQTL_plot_CPM_3donorColored(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel,"MEL") + theme(legend.position="none")
  plot.krt <- make_PBSeQTL_plot_CPM_3donorColored(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt,"KRT") + theme(legend.position="none")
  plot.frb <- make_PBSeQTL_plot_CPM_3donorColored(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb,"FRB") + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}
make_IFNeQTL_plot_CPM_3donorColored_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_IFNeQTL_plot_CPM_3donorColored(this.snp,this.gene,CPM.ifn.mel,IFN_knownVar.mel,"MEL") + theme(legend.position="none")
  plot.krt <- make_IFNeQTL_plot_CPM_3donorColored(this.snp,this.gene,CPM.ifn.krt,IFN_knownVar.krt,"KRT") + theme(legend.position="none")
  plot.frb <- make_IFNeQTL_plot_CPM_3donorColored(this.snp,this.gene,CPM.ifn.frb,IFN_knownVar.frb,"FRB") + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}
make_log2FC_plot_CPM_3donorColored_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_log2FC_plot_3donorColored(this.snp,this.gene,CPM.pbs.mel,CPM.ifn.mel,"MEL") + theme(legend.position="none")
  plot.krt <- make_log2FC_plot_3donorColored(this.snp,this.gene,CPM.pbs.krt,CPM.ifn.krt,"KRT") + theme(legend.position="none")
  plot.frb <- make_log2FC_plot_3donorColored(this.snp,this.gene,CPM.pbs.frb,CPM.ifn.frb,"FRB") + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}

# functions for 3 donor colored scatter plots in rank normalized CPM
make_PBSeQTL_plot_rankNormCPM_3donorColored <- function(this.snp,this.gene,CPM.pbs,PBS_knownVar, celltype) {
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
  this.CPM.PBS <- CPM.pbs[this.gene,] %>% rankNorm() %>% as.data.frame() %>%
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
    ggtitle(paste0("PBS","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," rank normalized value"))
  p1
  return(p1)
}
make_IFNeQTL_plot_rankNormCPM_3donorColored <- function(this.snp,this.gene,CPM.ifn,IFN_knownVar, celltype) {
  
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
  this.CPM.IFN <- CPM.ifn[this.gene,] %>% rankNorm() %>% as.data.frame() %>%
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
    ggtitle(paste0("IFN","\n",celltype,"\n",this.snp,":",this.gene)) +
    xlab("genotype") + ylab(paste0(this.gene," rank normalized value"))
  p1
  return(p1)
}
make_log2FC_plot_rankNormlog2FC_3donorColored <- function(this.snp,this.gene,CPM.pbs,CPM.ifn, celltype) {
  
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
  this.CPM$log2FC <- rankNorm( log2((this.CPM$IFN + 10) / (this.CPM$PBS + 10)) )
  this.CPM <- this.CPM %>% dplyr::select(-c("PBS","IFN"))
  
  # pivot to long format and then plot
  this.CPM.long <- this.CPM %>% pivot_longer(.,cols=c(log2FC))
  
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
    theme(legend.position="none") +
    ggtitle(paste0("log2FC","\n",celltype,"\n",this.snp,":",this.gene)) +
    ylab("rank normalized log2(CPM+10/CPM+10)")
  return(p1)
}
make_PBSeQTL_plot_rankNormCPM_3donorColored_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_PBSeQTL_plot_rankNormCPM_3donorColored(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel,"MEL") + theme(legend.position="none")
  plot.krt <- make_PBSeQTL_plot_rankNormCPM_3donorColored(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt,"KRT") + theme(legend.position="none")
  plot.frb <- make_PBSeQTL_plot_rankNormCPM_3donorColored(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb,"FRB") + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}
make_IFNeQTL_plot_rankNormCPM_3donorColored_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_IFNeQTL_plot_rankNormCPM_3donorColored(this.snp,this.gene,CPM.ifn.mel,IFN_knownVar.mel,"MEL") + theme(legend.position="none")
  plot.krt <- make_IFNeQTL_plot_rankNormCPM_3donorColored(this.snp,this.gene,CPM.ifn.krt,IFN_knownVar.krt,"KRT") + theme(legend.position="none")
  plot.frb <- make_IFNeQTL_plot_rankNormCPM_3donorColored(this.snp,this.gene,CPM.ifn.frb,IFN_knownVar.frb,"FRB") + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}
make_log2FC_plot_rankNormlog2FC_3donorColored_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_log2FC_plot_rankNormlog2FC_3donorColored(this.snp,this.gene,CPM.pbs.mel,CPM.ifn.mel,"MEL") + theme(legend.position="none")
  plot.krt <- make_log2FC_plot_rankNormlog2FC_3donorColored(this.snp,this.gene,CPM.pbs.krt,CPM.ifn.krt,"KRT") + theme(legend.position="none")
  plot.frb <- make_log2FC_plot_rankNormlog2FC_3donorColored(this.snp,this.gene,CPM.pbs.frb,CPM.ifn.frb,"FRB") + theme(legend.position="none")
  
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
}

# functions for paired reQTL plot
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
    ggtitle(paste0(celltype,"\n",this.snp,":",this.gene)) +
    theme(legend.position="none") +
    facet_wrap(~ genotype)
  
  return(p)
}
make_reQTL_plot_CPM_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.mel, CPM.ifn.mel,"MEL")
  plot.krt <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.krt, CPM.ifn.krt,"KRT")
  plot.frb <- make_reQTL_plot_CPM(this.snp,this.gene,CPM.pbs.frb, CPM.ifn.frb,"FRB")
  p = plot_grid(plot.mel, plot.krt, plot.frb, ncol = 3, nrow = 1)
  p
  return(p)
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


###### make plots ######
make_all_the_plots <- function(snp, gene, outDir) {

  # ancestry colored scatter plots
  scatterPlot_ancestryColored_PBSeQTL_CPM <- make_PBSeQTL_plot_CPM_ancestryColored_3cts(snp, gene)
  scatterPlot_ancestryColored_IFNeQTL_CPM <- make_IFNeQTL_plot_CPM_ancestryColored_3cts(snp, gene)
  
  scatterPlot_ancestryColored_PBSeQTL_rankNormCPM <- make_PBSeQTL_plot_rankNormCPM_ancestryColored_3cts(snp, gene)
  scatterPlot_ancestryColored_IFNeQTL_rankNormCPM <- make_IFNeQTL_plot_rankNormCPM_ancestryColored_3cts(snp, gene)
  
  # donor colored scatter plots
  scatterPlot_3donorColored_PBSeQTL_CPM <- make_PBSeQTL_plot_CPM_3donorColored_3cts(snp, gene)
  scatterPlot_3donorColored_IFNeQTL_CPM <- make_IFNeQTL_plot_CPM_3donorColored_3cts(snp, gene)
  scatterPlot_3donorColored_log2FC_CPM <- make_log2FC_plot_CPM_3donorColored_3cts(snp, gene)
  
  scatterPlot_3donorColored_PBSeQTL_rankNormCPM <- make_PBSeQTL_plot_rankNormCPM_3donorColored_3cts(snp, gene)
  scatterPlot_3donorColored_IFNeQTL_rankNormCPM <- make_IFNeQTL_plot_rankNormCPM_3donorColored_3cts(snp, gene)
  scatterPlot_3donorColored_log2FC_rankNormCPM <- make_log2FC_plot_rankNormlog2FC_3donorColored_3cts(snp, gene)
  
  # paired plots
  pairedPlot_3donorColored_CPM <- make_reQTL_plot_CPM_3cts(snp, gene)
  
  ###### save the plots in cluster ###### 
  options(scipen = 0,digits=3) # turn on scientific notation with 3 digits
  # ancestry colored scatter plots
  png(paste0(outDir,"/scatterPlot_ancestryColored_PBSeQTL_CPM.png"),width=1200,height=600)
  print(scatterPlot_ancestryColored_PBSeQTL_CPM); dev.off()
  png(paste0(outDir,"/scatterPlot_ancestryColored_IFNeQTL_CPM.png"),width=1200,height=600)
  print(scatterPlot_ancestryColored_IFNeQTL_CPM); dev.off()
  png(paste0(outDir,"/scatterPlot_ancestryColored_PBSeQTL_rankNormCPM.png"),width=1200,height=600)
  print(scatterPlot_ancestryColored_PBSeQTL_rankNormCPM); dev.off()
  png(paste0(outDir,"/scatterPlot_ancestryColored_IFNeQTL_rankNormCPM.png"),width=1200,height=600)
  print(scatterPlot_ancestryColored_IFNeQTL_rankNormCPM); dev.off()
  
  # donor colored scatter plots
  png(paste0(outDir,"/scatterPlot_3donorColored_PBSeQTL_CPM.png"),width=600,height=400)
  print(scatterPlot_3donorColored_PBSeQTL_CPM); dev.off()
  png(paste0(outDir,"/scatterPlot_3donorColored_IFNeQTL_CPM.png"),width=600,height=400)
  print(scatterPlot_3donorColored_IFNeQTL_CPM); dev.off()
  png(paste0(outDir,"/scatterPlot_3donorColored_log2FC_CPM.png"),width=600,height=400)
  print(scatterPlot_3donorColored_log2FC_CPM); dev.off()
  
  png(paste0(outDir,"/scatterPlot_3donorColored_PBSeQTL_rankNormCPM.png"),width=600,height=400)
  print(scatterPlot_3donorColored_PBSeQTL_rankNormCPM); dev.off()
  png(paste0(outDir,"/scatterPlot_3donorColored_IFNeQTL_rankNormCPM.png"),width=600,height=400)
  print(scatterPlot_3donorColored_IFNeQTL_rankNormCPM); dev.off()
  png(paste0(outDir,"/scatterPlot_3donorColored_log2FC_rankNormCPM.png"),width=600,height=400)
  print(scatterPlot_3donorColored_log2FC_rankNormCPM); dev.off()
  
  # paired plots
  png(paste0(outDir,"/pairedPlot_3donorColored_CPM.png"),width=600,height=400)
  print(pairedPlot_3donorColored_CPM); dev.off()
}

###### load snp gene pairs and make plots one pair at a time
snp_gene_pairs_F <- "~/Downloads/nl/human/skin/eQTLs/DREG/ancestry/KRT/PBSeQTLs_whose_ancestry_is_sig_feature.txt"
snp_gene_pairs <- data.table::fread(snp_gene_pairs_F, header=F) %>% 
  tidyr::separate(., col=V1, into=c("snp","gene"),sep="_") %>% 
  arrange(gene)

df <- snp_gene_pairs
for (i in 1:nrow(df)) {
  rm(this_snp, this_gene, this_outDir, genotypeF, genotype)
  this_snp <- df[i,"snp"]
  this_gene <- df[i,"gene"]
  system(paste0("mkdir ~/Downloads/nl/human/skin/eQTLs/DREG/ancestry/KRT/plot/",this_gene))
  system(paste0("mkdir ~/Downloads/nl/human/skin/eQTLs/DREG/ancestry/KRT/plot/",this_gene,"/",this_snp))
  this_outDir=paste0("~/Downloads/nl/human/skin/eQTLs/DREG/ancestry/KRT/plot/",this_gene,"/",this_snp)
  genotype_table="~/Downloads/nl/human/skin/eQTLs/DREG/ancestry/KRT/selected_snps_genotype.bed"
  genotype = data.table::fread(genotype_table,sep="\t")
  
  make_all_the_plots(this_snp, this_gene, this_outDir)
  system(paste0("chmod g+w ",this_outDir))
  system(paste0("chmod g+w ",this_outDir,"/*"))
  print(paste("[",as.character(i),"]",this_snp,this_gene,Sys.time()))
}


##### load modeling result: MEL ##### 
resultDir="~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results"
f.mel = data.table::fread(paste0(resultDir,"/MEL_modeling_stats_a10b10.txt"), fill=TRUE) %>% unique %>%
  set_colnames(result_table_colnames) %>% mutate(pair=paste0(.$SNP,"_",.$GENE)) %>% 
  mutate(celltype="melanocyte")
f.krt = data.table::fread(paste0(resultDir,"/KRT_modeling_stats_a10b10.txt"), fill=TRUE) %>% unique %>%
  set_colnames(result_table_colnames) %>% mutate(pair=paste0(.$SNP,"_",.$GENE)) %>% 
  mutate(celltype="keratinocyte")
f.frb = data.table::fread(paste0(resultDir,"/FRB_modeling_stats_a10b10.txt"), fill=TRUE) %>% unique %>%
  set_colnames(result_table_colnames) %>% mutate(pair=paste0(.$SNP,"_",.$GENE)) %>% 
  mutate(celltype="fibroblast")

f.output = rbind(f.mel, f.krt, f.frb)[,c(51,1:2,52,3:50)]
data.table::fwrite(f.output, file=paste0(resultDir,"/modeling_results_compiled_mastertable.txt"), 
                   sep = "\t", quote = F, col.names = T, row.names = F)

f.temp = f.output %>% 
  dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("rankNormlog2FC_minimalModel|rankNormCPM_minimalModel",colnames(f.output))])) %>%
  set_colnames(gsub("rankNorm|log2FC_|CPM_|","",colnames(.))) %>%
  set_colnames(gsub("minimalModel_","",colnames(.))) %>%
  set_colnames(gsub("_genotype","",colnames(.)))
reQTL = f.temp %>% dplyr::filter(reQTL_pval != "." & reQTL_pval < 0.00001) %>% mutate(tag="reQTL") %>% mutate(QTL=paste0(celltype,"_",tag)) %>% dplyr::select(c("pair","QTL")) # reQTL: (pval 2479, pperm 2936); PBSeQTL: (pval 3015, pperm 5049); IFNeQTL: (pval 3048, pperm 5418)
PBSeQTL = f.temp %>% dplyr::filter(PBSeQTL_pval != "." & PBSeQTL_pval < 0.00001) %>% mutate(tag="PBSeQTL") %>% mutate(QTL=paste0(celltype,"_",tag)) %>% dplyr::select(c("pair","QTL")) # reQTL: (pval 2479, pperm 2936); PBSeQTL: (pval 3015, pperm 5049); IFNeQTL: (pval 3048, pperm 5418)
IFNeQTL = f.temp %>% dplyr::filter(IFNeQTL_pval != "." & IFNeQTL_pval < 0.00001) %>% mutate(tag="IFNeQTL") %>% mutate(QTL=paste0(celltype,"_",tag)) %>% dplyr::select(c("pair","QTL")) # reQTL: (pval 2479, pperm 2936); PBSeQTL: (pval 3015, pperm 5049); IFNeQTL: (pval 3048, pperm 5418)
QTL_df = rbind(reQTL, PBSeQTL, IFNeQTL) %>% group_by(pair) %>%
  summarize(QTL_type=paste(QTL, collapse=", "))

f.minimalCPM = f.output %>% 
  dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("^log2FC_minimalModel|^CPM_minimalModel",colnames(f.output))])) %>%
  set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
  set_colnames(gsub("minimalModel_","",colnames(.))) %>%
  set_colnames(gsub("_genotype","",colnames(.))) %>%
  right_join(. , QTL_df, by = join_by(pair))
f.minimalCPM = f.minimalCPM[,c(1:4,17,5:16)]
data.table::fwrite(f.minimalCPM, file=paste0(resultDir,"/modeling_results_minimalModel_phenotypeCPM.txt"), 
                   sep = "\t", quote = F, col.names = T, row.names = F)

f.minimalrankNorm = f.output %>% 
  dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("rankNormlog2FC_minimalModel|rankNormCPM_minimalModel",colnames(f.output))])) %>%
  set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
  set_colnames(gsub("rankNorm","",colnames(.))) %>%
  set_colnames(gsub("minimalModel_","",colnames(.))) %>%
  set_colnames(gsub("_genotype","",colnames(.))) %>%
  right_join(. , QTL_df, by = join_by(pair))
f.minimalrankNorm = f.minimalrankNorm[,c(1:4,17,5:16)]
data.table::fwrite(f.minimalrankNorm, file=paste0(resultDir,"/modeling_results_minimalModel_phenotypeRankNormCPM.txt"), 
                   sep = "\t", quote = F, col.names = T, row.names = F)

f.featureSelectedCPM = f.output %>% 
  dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("^log2FC_featureSelected|^CPM_featureSelected",colnames(f.output))])) %>%
  set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
  set_colnames(gsub("featureSelected_","",colnames(.))) %>%
  set_colnames(gsub("_genotype","",colnames(.))) %>%
  right_join(. , QTL_df, by = join_by(pair))
f.featureSelectedCPM = f.featureSelectedCPM[,c(1:4,17,5:16)]
data.table::fwrite(f.featureSelectedCPM, file=paste0(resultDir,"/modeling_results_featureSelectedModel_phenotypeCPM.txt"), 
                   sep = "\t", quote = F, col.names = T, row.names = F)

f.featureSelectedRankNorm = f.output %>% 
  dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("rankNormlog2FC_featureSelected|rankNormCPM_featureSelected",colnames(f.output))])) %>%
  set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
  set_colnames(gsub("rankNorm","",colnames(.))) %>%
  set_colnames(gsub("featureSelected_","",colnames(.))) %>%
  set_colnames(gsub("_genotype","",colnames(.))) %>%
  right_join(. , QTL_df, by = join_by(pair))
f.featureSelectedRankNorm = f.featureSelectedRankNorm[,c(1:4,17,5:16)]
data.table::fwrite(f.featureSelectedRankNorm, file=paste0(resultDir,"/modeling_results_featureSelectedModel_phenotypeRankNormCPM.txt"), 
                   sep = "\t", quote = F, col.names = T, row.names = F)
