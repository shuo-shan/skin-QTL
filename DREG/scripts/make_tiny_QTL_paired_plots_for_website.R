#!/usr/bin/env Rscript
# written by Crystal Shan 09/2023
# This script makes 2-3KB sized png plot per SNP:gene pair.

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
#args = commandArgs(trailingOnly=TRUE)
#outDir=args[1]
#gene=args[2]
#snp=args[3]
# 
#gene="ERAP2"
#snp="rs2910686"

#outDir=paste0("~/Downloads/nl/human/skin/eQTLs/website/data/plots","/",gene,"/",snp)
outDir="~/Downloads/nl/human/skin/eQTLs/case_studies/DHX58/rs8081327"
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

make_reQTL_plot_CPM_3cts(snp,gene)
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


#save.image(file="~/Downloads/nl/human/skin/eQTLs/website/data/plotting_data.RData")
###### make plots ######
make_all_the_plots <- function(snp, gene, outDir) {

  # # ancestry colored scatter plots
  # scatterPlot_ancestryColored_PBSeQTL_CPM <- make_PBSeQTL_plot_CPM_ancestryColored_3cts(snp, gene)
  # scatterPlot_ancestryColored_IFNeQTL_CPM <- make_IFNeQTL_plot_CPM_ancestryColored_3cts(snp, gene)
  # 
  # scatterPlot_ancestryColored_PBSeQTL_rankNormCPM <- make_PBSeQTL_plot_rankNormCPM_ancestryColored_3cts(snp, gene)
  # scatterPlot_ancestryColored_IFNeQTL_rankNormCPM <- make_IFNeQTL_plot_rankNormCPM_ancestryColored_3cts(snp, gene)
  # 
  # # donor colored scatter plots
  # scatterPlot_3donorColored_PBSeQTL_CPM <- make_PBSeQTL_plot_CPM_3donorColored_3cts(snp, gene)
  # scatterPlot_3donorColored_IFNeQTL_CPM <- make_IFNeQTL_plot_CPM_3donorColored_3cts(snp, gene)
  # scatterPlot_3donorColored_log2FC_CPM <- make_log2FC_plot_CPM_3donorColored_3cts(snp, gene)
  # 
  # scatterPlot_3donorColored_PBSeQTL_rankNormCPM <- make_PBSeQTL_plot_rankNormCPM_3donorColored_3cts(snp, gene)
  # scatterPlot_3donorColored_IFNeQTL_rankNormCPM <- make_IFNeQTL_plot_rankNormCPM_3donorColored_3cts(snp, gene)
  # scatterPlot_3donorColored_log2FC_rankNormCPM <- make_log2FC_plot_rankNormlog2FC_3donorColored_3cts(snp, gene)
  # 
  # paired plots
  pairedPlot_3donorColored_CPM <- make_reQTL_plot_CPM_3cts(snp, gene)
  
  ###### save the plots in cluster ###### 
  options(scipen = 0,digits=3) # turn on scientific notation with 3 digits
  # # ancestry colored scatter plots
  # png(paste0(outDir,"/scatterPlot_ancestryColored_PBSeQTL_CPM.png"),width=1200,height=600)
  # print(scatterPlot_ancestryColored_PBSeQTL_CPM); dev.off()
  # png(paste0(outDir,"/scatterPlot_ancestryColored_IFNeQTL_CPM.png"),width=1200,height=600)
  # print(scatterPlot_ancestryColored_IFNeQTL_CPM); dev.off()
  # png(paste0(outDir,"/scatterPlot_ancestryColored_PBSeQTL_rankNormCPM.png"),width=1200,height=600)
  # print(scatterPlot_ancestryColored_PBSeQTL_rankNormCPM); dev.off()
  # png(paste0(outDir,"/scatterPlot_ancestryColored_IFNeQTL_rankNormCPM.png"),width=1200,height=600)
  # print(scatterPlot_ancestryColored_IFNeQTL_rankNormCPM); dev.off()
  # 
  # # donor colored scatter plots
  # png(paste0(outDir,"/scatterPlot_3donorColored_PBSeQTL_CPM.png"),width=600,height=400)
  # print(scatterPlot_3donorColored_PBSeQTL_CPM); dev.off()
  # png(paste0(outDir,"/scatterPlot_3donorColored_IFNeQTL_CPM.png"),width=600,height=400)
  # print(scatterPlot_3donorColored_IFNeQTL_CPM); dev.off()
  # png(paste0(outDir,"/scatterPlot_3donorColored_log2FC_CPM.png"),width=600,height=400)
  # print(scatterPlot_3donorColored_log2FC_CPM); dev.off()
  # 
  # png(paste0(outDir,"/scatterPlot_3donorColored_PBSeQTL_rankNormCPM.png"),width=600,height=400)
  # print(scatterPlot_3donorColored_PBSeQTL_rankNormCPM); dev.off()
  # png(paste0(outDir,"/scatterPlot_3donorColored_IFNeQTL_rankNormCPM.png"),width=600,height=400)
  # print(scatterPlot_3donorColored_IFNeQTL_rankNormCPM); dev.off()
  # png(paste0(outDir,"/scatterPlot_3donorColored_log2FC_rankNormCPM.png"),width=600,height=400)
  # print(scatterPlot_3donorColored_log2FC_rankNormCPM); dev.off()
  # 
  # paired plots
  png(paste0(outDir,"/pairedPlot_3donorColored_CPM_",snp,"_",gene,".png"),width=310,height=150,res=90,type="Xlib")
  print(pairedPlot_3donorColored_CPM); dev.off()
}

###### load snp gene pairs and make plots one pair at a time
bigtable.krt <- readRDS("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_krt.rds")
bigtable.mel <- readRDS("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_mel.rds")
bigtable.frb <- readRDS("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_frb.rds")
betacomp.krt <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/KRT/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE)) %>%
  dplyr::filter(p.betaComp10KPermut<0.01)
betacomp.mel <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/MEL/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE)) %>%
  dplyr::filter(p.betaComp10KPermut<0.01)
betacomp.frb <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/FRB/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE)) %>%
  dplyr::filter(p.betaComp10KPermut<0.01)

p.betaComp < 0.05 & p.permute < 0.001

snp_gene_pairs_F <- "~/Downloads/nl/human/skin/eQTLs/chromatin/ANNOVAR/reQTLs_1E-05/temp_snp_gene_pairs.txt"
snp_gene_pairs <- data.table::fread(snp_gene_pairs_F, header=F) %>% 
  tidyr::separate(., col=V1, into=c("snp","gene"),sep="_") %>% 
  arrange(gene)

df <- snp_gene_pairs
for (i in 1:nrow(df)) {
  #rm(this_snp, this_gene, this_outDir, genotypeF, genotype)
  this_snp <- df[i,"snp"]
  this_gene <- df[i,"gene"]
  this_outDir=paste0("~/Downloads/nl/human/skin/eQTLs/website/data/plots","/",this_gene,"/",this_snp)
  genotypeF=paste0("~/Downloads/nl/human/skin/eQTLs/website/data/plots/",this_gene,"/",this_snp,"/genotype.bed")
  genotype = data.table::fread(genotypeF,header=TRUE,sep="\t")
  
  # Wrap the function call in tryCatch to handle potential errors
  tryCatch({
    make_all_the_plots(this_snp, this_gene, this_outDir)
  }, error = function(e) {
    # Handle error
    cat("An error occurred with SNP:", this_snp, "Gene:", this_gene, "\nError message:", e$message, "\n")
  })
  
  system(paste0("chmod g+w ",this_outDir))
  system(paste0("chmod g+w ",this_outDir,"/*"))
  print(paste("[",as.character(i),"]",this_snp,this_gene,Sys.time()))
}


##### load modeling result
f.mel = data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt")
f.mel = f.mel %>% dplyr::select("GENE","SNP","featureSelected_reQTL_genotype_pval","featureSelected_reQTL_genotype_beta",
                                "featureSelected_PBSeQTL_genotype_pval","featureSelected_PBSeQTL_genotype_beta",
                                "featureSelected_IFNeQTL_genotype_pval","featureSelected_IFNeQTL_genotype_beta")
colnames(f.mel) = gsub("featureSelected_","",colnames(f.mel))
colnames(f.mel) = gsub("_genotype","",colnames(f.mel))
f.mel$pair = paste0(f.mel$SNP,"_",f.mel$GENE)

f.mel.reQTL = f.mel %>%
  dplyr::filter(reQTL_pval!=".") %>%
  dplyr::mutate(reQTL_pval = as.numeric(reQTL_pval)) %>%
  dplyr::filter(reQTL_pval< 0.00001) %>%
  dplyr::mutate(key="MEL_reQTL")
f.mel.PBSeQTL = f.mel %>%
  dplyr::filter(PBSeQTL_pval!=".") %>%
  dplyr::mutate(PBSeQTL_pval = as.numeric(PBSeQTL_pval)) %>%
  dplyr::filter(PBSeQTL_pval< 0.00001) %>%
  dplyr::mutate(key="MEL_PBSeQTL")
f.mel.IFNeQTL = f.mel %>%
  dplyr::filter(IFNeQTL_pval!=".") %>%
  dplyr::mutate(IFNeQTL_pval = as.numeric(IFNeQTL_pval)) %>%
  dplyr::filter(IFNeQTL_pval< 0.00001) %>%
  dplyr::mutate(key="MEL_IFNeQTL")
f.mel.QTLs = rbind(f.mel.reQTL, f.mel.PBSeQTL, f.mel.IFNeQTL)
col1 = f.mel.QTLs[,c("pair","key")] %>%
  group_by(pair) %>% 
  summarize(key=paste(key,collapse=',')) 
col2 = f.mel.QTLs[!duplicated(f.mel.QTLs$pair),-"key"]
f.mel.QTLs.annotated = left_join(col1,col2, by="pair") %>% mutate(celltype="MEL")
rm(f.mel, f.mel.reQTL, f.mel.PBSeQTL, f.mel.IFNeQTL, col1, col2, f.mel.QTLs)

# FRB
f.frb = data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/FRB_new/masteroutput_all_with_colnames.txt")
f.frb = f.frb %>% dplyr::select("GENE","SNP","featureSelected_reQTL_genotype_pval","featureSelected_reQTL_genotype_beta",
                                "featureSelected_PBSeQTL_genotype_pval","featureSelected_PBSeQTL_genotype_beta",
                                "featureSelected_IFNeQTL_genotype_pval","featureSelected_IFNeQTL_genotype_beta")
colnames(f.frb) = gsub("featureSelected_","",colnames(f.frb))
colnames(f.frb) = gsub("_genotype","",colnames(f.frb))
f.frb$pair = paste0(f.frb$SNP,"_",f.frb$GENE)

f.frb.reQTL = f.frb %>%
  dplyr::filter(reQTL_pval!=".") %>%
  dplyr::mutate(reQTL_pval = as.numeric(reQTL_pval)) %>%
  dplyr::filter(reQTL_pval< 0.00001) %>%
  dplyr::mutate(key="FRB_reQTL")
f.frb.PBSeQTL = f.frb %>%
  dplyr::filter(PBSeQTL_pval!=".") %>%
  dplyr::mutate(PBSeQTL_pval = as.numeric(PBSeQTL_pval)) %>%
  dplyr::filter(PBSeQTL_pval< 0.00001) %>%
  dplyr::mutate(key="FRB_PBSeQTL")
f.frb.IFNeQTL = f.frb %>%
  dplyr::filter(IFNeQTL_pval!=".") %>%
  dplyr::mutate(IFNeQTL_pval = as.numeric(IFNeQTL_pval)) %>%
  dplyr::filter(IFNeQTL_pval< 0.00001) %>%
  dplyr::mutate(key="FRB_IFNeQTL")
f.frb.QTLs = rbind(f.frb.reQTL, f.frb.PBSeQTL, f.frb.IFNeQTL)
col1 = f.frb.QTLs[,c("pair","key")] %>%
  group_by(pair) %>% 
  summarize(key=paste(key,collapse=',')) 
col2 = f.frb.QTLs[!duplicated(f.frb.QTLs$pair),-"key"]
f.frb.QTLs.annotated = left_join(col1,col2, by="pair") %>% mutate(celltype="FRB")
rm(f.frb, f.frb.reQTL, f.frb.PBSeQTL, f.frb.IFNeQTL, col1, col2, f.frb.QTLs)

#### KRT
f.krt = data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt")
f.krt = f.krt %>% dplyr::select("GENE","SNP","featureSelected_reQTL_genotype_pval","featureSelected_reQTL_genotype_beta",
                                "featureSelected_PBSeQTL_genotype_pval","featureSelected_PBSeQTL_genotype_beta",
                                "featureSelected_IFNeQTL_genotype_pval","featureSelected_IFNeQTL_genotype_beta")
colnames(f.krt) = gsub("featureSelected_","",colnames(f.krt))
colnames(f.krt) = gsub("_genotype","",colnames(f.krt))
f.krt$pair = paste0(f.krt$SNP,"_",f.krt$GENE)

f.krt.reQTL = f.krt %>%
  dplyr::filter(reQTL_pval!=".") %>%
  dplyr::mutate(reQTL_pval = as.numeric(reQTL_pval)) %>%
  dplyr::filter(reQTL_pval< 0.00001) %>%
  dplyr::mutate(key="KRT_reQTL")
f.krt.PBSeQTL = f.krt %>%
  dplyr::filter(PBSeQTL_pval!=".") %>%
  dplyr::mutate(PBSeQTL_pval = as.numeric(PBSeQTL_pval)) %>%
  dplyr::filter(PBSeQTL_pval< 0.00001) %>%
  dplyr::mutate(key="KRT_PBSeQTL")
f.krt.IFNeQTL = f.krt %>%
  dplyr::filter(IFNeQTL_pval!=".") %>%
  dplyr::mutate(IFNeQTL_pval = as.numeric(IFNeQTL_pval)) %>%
  dplyr::filter(IFNeQTL_pval< 0.00001) %>%
  dplyr::mutate(key="KRT_IFNeQTL")
f.krt.QTLs = rbind(f.krt.reQTL, f.krt.PBSeQTL, f.krt.IFNeQTL)
col1 = f.krt.QTLs[,c("pair","key")] %>%
  group_by(pair) %>% 
  summarize(key=paste(key,collapse=',')) 
col2 = f.krt.QTLs[!duplicated(f.krt.QTLs$pair),-"key"]
f.krt.QTLs.annotated = left_join(col1,col2, by="pair") %>% mutate(celltype="KRT")
rm(f.krt, f.krt.reQTL, f.krt.PBSeQTL, f.krt.IFNeQTL, col1, col2, f.krt.QTLs)

## bind results together:
QTLs = rbind(f.frb.QTLs.annotated, f.krt.QTLs.annotated, f.mel.QTLs.annotated) %>% arrange(pair)
write(unique(QTLs$pair),"~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results/SNPGenePairs_QTLs_1E-5PvalCutoff_3cts.txt")
data.table::fwrite(QTLs,file="~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results/modelingResults_featureSelectedRankNormCPM_QTLs_1E-5PvalCutoff_3cts.txt",quote=F,sep="\t",row.names=F,col.names=T)
# 
# ##### load modeling result: MEL ##### 
# resultDir="~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results"
# f.mel = data.table::fread(paste0(resultDir,"/MEL_modeling_stats_a10b10.txt"), fill=TRUE) %>% unique %>%
#   set_colnames(result_table_colnames) %>% mutate(pair=paste0(.$SNP,"_",.$GENE)) %>% 
#   mutate(celltype="melanocyte")
# f.krt = data.table::fread(paste0(resultDir,"/KRT_modeling_stats_a10b10.txt"), fill=TRUE) %>% unique %>%
#   set_colnames(result_table_colnames) %>% mutate(pair=paste0(.$SNP,"_",.$GENE)) %>% 
#   mutate(celltype="keratinocyte")
# f.frb = data.table::fread(paste0(resultDir,"/FRB_modeling_stats_a10b10.txt"), fill=TRUE) %>% unique %>%
#   set_colnames(result_table_colnames) %>% mutate(pair=paste0(.$SNP,"_",.$GENE)) %>% 
#   mutate(celltype="fibroblast")
# 
# f.output = rbind(f.mel, f.krt, f.frb)[,c(51,1:2,52,3:50)]
# data.table::fwrite(f.output, file=paste0(resultDir,"/modeling_results_compiled_mastertable.txt"), 
#                    sep = "\t", quote = F, col.names = T, row.names = F)
# 
# f.temp = f.output %>% 
#   dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("rankNormlog2FC_minimalModel|rankNormCPM_minimalModel",colnames(f.output))])) %>%
#   set_colnames(gsub("rankNorm|log2FC_|CPM_|","",colnames(.))) %>%
#   set_colnames(gsub("minimalModel_","",colnames(.))) %>%
#   set_colnames(gsub("_genotype","",colnames(.)))
# reQTL = f.temp %>% dplyr::filter(reQTL_pval != "." & reQTL_pval < 0.00001) %>% mutate(tag="reQTL") %>% mutate(QTL=paste0(celltype,"_",tag)) %>% dplyr::select(c("pair","QTL")) # reQTL: (pval 2479, pperm 2936); PBSeQTL: (pval 3015, pperm 5049); IFNeQTL: (pval 3048, pperm 5418)
# PBSeQTL = f.temp %>% dplyr::filter(PBSeQTL_pval != "." & PBSeQTL_pval < 0.00001) %>% mutate(tag="PBSeQTL") %>% mutate(QTL=paste0(celltype,"_",tag)) %>% dplyr::select(c("pair","QTL")) # reQTL: (pval 2479, pperm 2936); PBSeQTL: (pval 3015, pperm 5049); IFNeQTL: (pval 3048, pperm 5418)
# IFNeQTL = f.temp %>% dplyr::filter(IFNeQTL_pval != "." & IFNeQTL_pval < 0.00001) %>% mutate(tag="IFNeQTL") %>% mutate(QTL=paste0(celltype,"_",tag)) %>% dplyr::select(c("pair","QTL")) # reQTL: (pval 2479, pperm 2936); PBSeQTL: (pval 3015, pperm 5049); IFNeQTL: (pval 3048, pperm 5418)
# QTL_df = rbind(reQTL, PBSeQTL, IFNeQTL) %>% group_by(pair) %>%
#   summarize(QTL_type=paste(QTL, collapse=", "))
# 
# f.minimalCPM = f.output %>% 
#   dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("^log2FC_minimalModel|^CPM_minimalModel",colnames(f.output))])) %>%
#   set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
#   set_colnames(gsub("minimalModel_","",colnames(.))) %>%
#   set_colnames(gsub("_genotype","",colnames(.))) %>%
#   right_join(. , QTL_df, by = join_by(pair))
# f.minimalCPM = f.minimalCPM[,c(1:4,17,5:16)]
# data.table::fwrite(f.minimalCPM, file=paste0(resultDir,"/modeling_results_minimalModel_phenotypeCPM.txt"), 
#                    sep = "\t", quote = F, col.names = T, row.names = F)
# 
# f.minimalrankNorm = f.output %>% 
#   dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("rankNormlog2FC_minimalModel|rankNormCPM_minimalModel",colnames(f.output))])) %>%
#   set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
#   set_colnames(gsub("rankNorm","",colnames(.))) %>%
#   set_colnames(gsub("minimalModel_","",colnames(.))) %>%
#   set_colnames(gsub("_genotype","",colnames(.))) %>%
#   right_join(. , QTL_df, by = join_by(pair))
# f.minimalrankNorm = f.minimalrankNorm[,c(1:4,17,5:16)]
# data.table::fwrite(f.minimalrankNorm, file=paste0(resultDir,"/modeling_results_minimalModel_phenotypeRankNormCPM.txt"), 
#                    sep = "\t", quote = F, col.names = T, row.names = F)
# 
# f.featureSelectedCPM = f.output %>% 
#   dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("^log2FC_featureSelected|^CPM_featureSelected",colnames(f.output))])) %>%
#   set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
#   set_colnames(gsub("featureSelected_","",colnames(.))) %>%
#   set_colnames(gsub("_genotype","",colnames(.))) %>%
#   right_join(. , QTL_df, by = join_by(pair))
# f.featureSelectedCPM = f.featureSelectedCPM[,c(1:4,17,5:16)]
# data.table::fwrite(f.featureSelectedCPM, file=paste0(resultDir,"/modeling_results_featureSelectedModel_phenotypeCPM.txt"), 
#                    sep = "\t", quote = F, col.names = T, row.names = F)
# 
# f.featureSelectedRankNorm = f.output %>% 
#   dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("rankNormlog2FC_featureSelected|rankNormCPM_featureSelected",colnames(f.output))])) %>%
#   set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
#   set_colnames(gsub("rankNorm","",colnames(.))) %>%
#   set_colnames(gsub("featureSelected_","",colnames(.))) %>%
#   set_colnames(gsub("_genotype","",colnames(.))) %>%
#   right_join(. , QTL_df, by = join_by(pair))
# f.featureSelectedRankNorm = f.featureSelectedRankNorm[,c(1:4,17,5:16)]
# data.table::fwrite(f.featureSelectedRankNorm, file=paste0(resultDir,"/modeling_results_featureSelectedModel_phenotypeRankNormCPM.txt"), 
#                    sep = "\t", quote = F, col.names = T, row.names = F)
