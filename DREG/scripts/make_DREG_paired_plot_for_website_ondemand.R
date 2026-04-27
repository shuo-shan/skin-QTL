#!/usr/bin/env Rscript
# written by Crystal Shan 09/2023

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)

# parsing input arguments
gene=args[2]
snp=args[3]
# 
# gene="DHX58"
# snp="rs739636"

# load genotype
genotype_table="~/Downloads/nl/human/skin/eQTLs/DREG/master_filtered_genotype.bed"
genotype = data.table::fread(genotype_table,header=TRUE,sep="\t")

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
CPM.pbs.mel <- read.table(CPM_table_PBS.mel,sep="\t",header=TRUE) %>% column_to_rownames("X")
CPM.ifn.mel <- read.table(CPM_table_IFN.mel,sep="\t",header=TRUE) %>% column_to_rownames("X")
# load KRT data
CPM_table_PBS.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt"
CPM.pbs.krt <- read.table(CPM_table_PBS.krt,header=TRUE,sep="\t") %>% column_to_rownames("X")
CPM.ifn.krt <- read.table(CPM_table_IFN.krt,header=TRUE,sep="\t") %>% column_to_rownames("X")
# load FRB data
CPM_table_PBS.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt"
CPM.pbs.frb <- read.table(CPM_table_PBS.frb,header=TRUE,sep="\t") %>% column_to_rownames("X")
CPM.ifn.frb <- read.table(CPM_table_IFN.frb,header=TRUE,sep="\t") %>% column_to_rownames("X")

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


###### load snp gene pairs and make plots one pair at a time
make_reQTL_plot_CPM_3cts(snp, gene)


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


f.featureSelectedRankNorm = f.output %>% 
  dplyr::select(c("pair","SNP","GENE","celltype",colnames(f.output)[grepl("rankNormlog2FC_featureSelected|rankNormCPM_featureSelected",colnames(f.output))])) %>%
  set_colnames(gsub("log2FC_|CPM_|","",colnames(.))) %>%
  set_colnames(gsub("rankNorm","",colnames(.))) %>%
  set_colnames(gsub("featureSelected_","",colnames(.))) %>%
  set_colnames(gsub("_genotype","",colnames(.))) %>%
  right_join(. , QTL_df, by = join_by(pair))
f.featureSelectedRankNorm = f.featureSelectedRankNorm[,c(1:4,17,5:16)]
