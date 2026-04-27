library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
args = commandArgs(trailingOnly=TRUE)
dir=args[1]
fastQTLResF=args[2]
edSiteCountF=args[3]
edSiteLevelF=args[4]
genotypeBedF=args[5]

dir="~/Downloads/nl/human/skin/eQTLs/edQTL/output"
fastQTLResF=paste0(dir,"/","foreskin.edMat.10cov.20samps.noXYM.qqnorm.nominal")
edSiteCountF=paste0(dir,"/","foreskin.edMat.10cov.15samps.txt")
edSiteLevelF=paste0(dir,"/","foreskin.edMat.10cov.20samps.noXYM.qqnorm.bed")
genotypeBedF=paste0(dir,"/","genotype.bed")

# fetch the editing-site and edQTL list
fastQTLRes = data.table::fread(fastQTLResF)
dict = fastQTLRes[,c("V1","V6","V7","V9","V11")] %>% arrange(V11)
colnames(dict) = c("edSite","edQTL","distance","slope","padj")
dict$pair = paste0(dict$edQTL,":",dict$edSite)

# fetch the editing-site count of edited reads and total reads
edSiteCount = data.table::fread(edSiteCountF) %>% 
  tidyr::separate(. , col=chrom, into=c("chr","start","end"),sep=":") %>%
  mutate(edSite=paste0(chr,"_",end)) %>%
  dplyr::filter(edSite %in% dict$edSite) %>%
  set_rownames(.$edSite) %>%
  dplyr::select(-c("chr","start","end","edSite"))

# fetch the editing-site editing level
edSiteLevel = data.table::fread(edSiteLevelF) %>% dplyr::filter(name %in% dict$edSite)
rownames(edSiteLevel) = edSiteLevel$name
edSiteLevel = edSiteLevel %>% dplyr::select(-c("#chr","start","end","name"))

# fetch the editing-site closest expressed genes
edSite_gene_lookup=data.table::fread(paste0(dir,"/edsites_closest_expressed_genes.txt")) %>% 
  dplyr::filter(V1 %in% dict$edSite)
# load MEL data
CPM_table_PBS.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/CPM_expressedGenes_IFN.txt"
PBS_knownVar_file.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/metadata_PBS.txt"
IFN_knownVar_file.mel="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/metadata_IFN.txt"
CPM.pbs.mel <- read.table(CPM_table_PBS.mel,sep="\t",header=TRUE) 
CPM.ifn.mel <- read.table(CPM_table_IFN.mel,sep="\t",header=TRUE) 
PBS_knownVar.mel <- read.table(PBS_knownVar_file.mel, header=TRUE,sep="\t") 
IFN_knownVar.mel <- read.table(IFN_knownVar_file.mel, header=TRUE,sep="\t") 
# load KRT data
CPM_table_PBS.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt"
PBS_knownVar_file.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/metadata_PBS.txt"
IFN_knownVar_file.krt="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/metadata_IFN.txt"
CPM.pbs.krt <- read.table(CPM_table_PBS.krt,header=TRUE,sep="\t") 
CPM.ifn.krt <- read.table(CPM_table_IFN.krt,header=TRUE,sep="\t") 
PBS_knownVar.krt <- read.table(PBS_knownVar_file.krt, sep="\t",header=TRUE) 
IFN_knownVar.krt <- read.table(IFN_knownVar_file.krt, sep="\t",header=TRUE) 
# load FRB data
CPM_table_PBS.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/CPM_expressedGenes_PBS.txt"
CPM_table_IFN.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt"
PBS_knownVar_file.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/metadata_PBS.txt"
IFN_knownVar_file.frb="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/metadata_IFN.txt"
CPM.pbs.frb <- read.table(CPM_table_PBS.frb,header=TRUE,sep="\t") 
CPM.ifn.frb <- read.table(CPM_table_IFN.frb,header=TRUE,sep="\t")
PBS_knownVar.frb <- read.table(PBS_knownVar_file.frb, header=TRUE,sep="\t") 
IFN_knownVar.frb <- read.table(IFN_knownVar_file.frb, header=TRUE,sep="\t") 


# fetch the edQTL genotype
genotypeBed = data.table::fread(genotypeBedF)
colnames(genotypeBed) = colnames(genotypeBed) %>% gsub("ID","edQTL",.) %>% gsub(".GT","",.) 
genotypeBed=genotypeBed[!duplicated(genotypeBed$edQTL),]

# functions
make_edQTL_plot_qqnorm <- function(this.edQTL,this.edSite) {
  # pick the edQTL
  this.edQTL.genotype <- genotypeBed %>% filter(edQTL==this.edQTL)
  edQTL.ref=this.edQTL.genotype$REF
  edQTL.alt=this.edQTL.genotype$ALT
  
  this.genotype=this.edQTL.genotype %>% 
    dplyr::select(-c("CHROM","START","END","edQTL","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(edQTL.ref,edQTL.ref), . ) %>%
    gsub("0/1",paste0(edQTL.ref,edQTL.alt), . ) %>%
    gsub("1/1",paste0(edQTL.alt,edQTL.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype = this.genotype %>% na.omit()
  rm(this.edQTL.genotype)
  
  # pick the edSite
  this.edSite.level <- edSiteLevel[which(rownames(edSiteLevel)==this.edSite),] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","editingLevel")) %>% 
    inner_join( . , this.genotype, by="donor")
  this.edSite.level$genotype <- ordered(this.edSite.level$genotype, 
                                   levels <- c(paste0(edQTL.ref,edQTL.ref),
                                               paste0(edQTL.ref,edQTL.alt),
                                               paste0(edQTL.alt,edQTL.alt)))
  this.edSite.level$color <- ifelse(this.edSite.level$donor=="F25","F25",
                               ifelse(this.edSite.level$donor=="F49","F49",
                                      ifelse(this.edSite.level$donor=="F55","F55","-")))
  this.edSite.level$color <- ordered(this.edSite.level$color, levels=c("F25","F49","F55","-"))
  
  # make scatter plots
  p1<-ggplot(this.edSite.level, aes(x=genotype, y=editingLevel, color=color)) +
    geom_jitter(size=3,width=0.1,alpha=0.8) +
    ggtitle(paste0(this.edQTL,":",this.edSite)) +
    scale_color_manual(breaks=c("F25","F49","F55","-"),
                       values = c("blue", "brown", "black","#999999")) +
    xlab("genotype") + ylab(paste0(this.edSite," editing level qqnorm score")) +
    theme(legend.position = "none")
  p1
  return(p1)
}
make_edQTL_plot_ratio <- function(this.edQTL,this.edSite) {
  # pick the edQTL
  this.edQTL.genotype <- genotypeBed %>% filter(edQTL==this.edQTL)
  edQTL.ref=this.edQTL.genotype[1,]$REF
  edQTL.alt=this.edQTL.genotype[1,]$ALT
  
  this.genotype=this.edQTL.genotype[1,] %>% 
    dplyr::select(-c("CHROM","START","END","edQTL","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(edQTL.ref,edQTL.ref), . ) %>%
    gsub("0/1",paste0(edQTL.ref,edQTL.alt), . ) %>%
    gsub("1/1",paste0(edQTL.alt,edQTL.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype = this.genotype %>% na.omit()
  rm(this.edQTL.genotype)
  
  # pick the edSite
  this.edSite.count.temp <- edSiteCount[which(rownames(edSiteCount)==this.edSite),] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","editingCount")) 
  this.edSite.count <- this.edSite.count.temp %>% 
    inner_join( . , this.genotype, by="donor") %>%
    tidyr::separate(. , col=editingCount, into=c("editedReads","totalReads"),sep="/") %>%
    mutate(editingRatio = as.numeric(editedReads)/as.numeric(totalReads)) 
  
  this.edSite.count$genotype <- ordered(this.edSite.count$genotype, 
                                        levels <- c(paste0(edQTL.ref,edQTL.ref),
                                                    paste0(edQTL.ref,edQTL.alt),
                                                    paste0(edQTL.alt,edQTL.alt)))
  this.edSite.count$color <- ifelse(this.edSite.count$donor=="F25","F25",
                                    ifelse(this.edSite.count$donor=="F49","F49",
                                           ifelse(this.edSite.count$donor=="F55","F55","-")))
  this.edSite.count$color <- ordered(this.edSite.count$color, levels=c("F25","F49","F55","-"))
  
  # make scatter plots
  p1<-ggplot(this.edSite.count, aes(x=genotype, y=editingRatio, color=color)) +
    geom_jitter(size=3,width=0.1,alpha=0.8) +
    ggtitle(paste0(this.edQTL,":",this.edSite)) +
    scale_color_manual(breaks=c("F25","F49","F55","-"),
                       values = c("blue", "brown", "black","#999999")) +
    xlab("genotype") + ylab(paste0(this.edSite," editing ratio")) +
    theme(legend.position = "none")
  p1
  return(p1)
}
make_edSite_coverage_histogram <- function(this.edSite) {
  # pick the edSite
  this.edSite.count <- edSiteCount[which(rownames(edSiteCount)==this.edSite),] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","editingCount")) %>% 
    tidyr::separate(. , col=editingCount, into=c("editedReads","totalReads"),sep="/") %>%
    arrange(totalReads) 
  this.edSite.count$donor <- factor(this.edSite.count$donor, levels=this.edSite.count$donor)
  
  this.edSite.count.long <- this.edSite.count %>%
    pivot_longer( . , c("editedReads","totalReads"),names_to="type") %>%
    mutate(value = as.numeric(value))
  
  p1 <- ggplot(this.edSite.count.long, aes(x=value,y=donor,fill=type)) +
    geom_bar(position="stack",stat="identity")
  p1
  return(p1)
}
make_edQTL_plot_coverage_histogram <- function(this.edQTL, this.edSite) {
  # pick the edQTL
  this.edQTL.genotype <- genotypeBed %>% filter(edQTL==this.edQTL)
  edQTL.ref=this.edQTL.genotype$REF
  edQTL.alt=this.edQTL.genotype$ALT
  
  this.genotype=this.edQTL.genotype %>% 
    dplyr::select(-c("CHROM","START","END","edQTL","REF","ALT")) %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
  
  this.genotype$genotype = this.genotype$genotype %>%
    gsub("0/0",paste0(edQTL.ref,edQTL.ref), . ) %>%
    gsub("0/1",paste0(edQTL.ref,edQTL.alt), . ) %>%
    gsub("1/1",paste0(edQTL.alt,edQTL.alt), . ) %>%
    gsub("./.",NA, . )
  this.genotype = this.genotype %>% na.omit()
  rm(this.edQTL.genotype)
  
  # pick the edSite
  this.edSite.count <- edSiteCount[which(rownames(edSiteCount)==this.edSite),] %>% t() %>% as.data.frame() %>%
    rownames_to_column("donor") %>% set_colnames(c("donor","editingCount")) %>% 
    tidyr::separate(. , col=editingCount, into=c("editedReads","totalReads"),sep="/") %>%
    arrange(totalReads) 
  this.edSite.count$editedReads <- as.numeric(this.edSite.count$editedReads)
  this.edSite.count$totalReads <- as.numeric(this.edSite.count$totalReads)
  this.edSite.count$nonEditedReads = this.edSite.count$totalReads - this.edSite.count$editedReads
  this.edSite.count$ratio <- paste0(this.edSite.count$donor," (",this.edSite.count$editedReads,"/",this.edSite.count$totalReads,")")
  this.edSite.count$donor <- factor(this.edSite.count$donor, levels=this.edSite.count$donor)
  this.edSite.count <- this.edSite.count %>% dplyr::select(-totalReads)
  
  
  this.edSite.count.long <- this.edSite.count %>%
    pivot_longer( . , c("editedReads","nonEditedReads"),names_to="type") %>%
    mutate(value = as.numeric(value)) %>%
    left_join( . , this.genotype, by="donor")
  
  this.edSite.count.long$genotype <- factor(this.edSite.count.long$genotype, levels=c(paste0(edQTL.ref,edQTL.ref),
                                                                                      paste0(edQTL.ref,edQTL.alt),
                                                                                      paste0(edQTL.alt,edQTL.alt),
                                                                                      NA))
  
  p1 <- ggplot(this.edSite.count.long, aes(x=value,y=donor,fill=type)) +
    geom_bar(position="stack",stat="identity") +
    scale_y_discrete(breaks=this.edSite.count.long$donor,
                     labels=this.edSite.count.long$ratio) +
    facet_wrap(~genotype)
  p1
  return(p1)
}
make_PBSeQTL_plot_CPM_nocolor <- function(this.snp,this.gene,CPM.pbs,PBS_knownVar) {
  # pick the snp
  this.snp.genotype <- genotypeBed %>% filter(edQTL==this.snp)
  snp.ref=this.snp.genotype$REF
  snp.alt=this.snp.genotype$ALT
  
  this.genotype=this.snp.genotype %>% 
    dplyr::select(-c("CHROM","START","END","edQTL","REF","ALT")) %>% t() %>% as.data.frame() %>%
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
  this.CPM.PBS$color <- ifelse(this.CPM.PBS$donor=="F25","F25",
                               ifelse(this.CPM.PBS$donor=="F49","F49",
                                      ifelse(this.CPM.PBS$donor=="F55","F55","-")))
  this.CPM.PBS$color <- ordered(this.CPM.PBS$color, levels=c("F25","F49","F55","-"))
  
  # make scatter plots
  p1<-ggplot(this.CPM.PBS, aes(x=genotype, y=PBS, color=color)) +
    geom_jitter(size=3,width=0.1) +
    #ggtitle(paste0(this.snp,":",this.gene)) +
    scale_color_manual(values = c("blue", "brown", "black","#999999")) +
    xlab("genotype") + ylab(paste0(this.gene," CPM")) +
    theme(legend.position = "none")
  p1
  return(p1)
}
make_PBSeQTL_plot_CPM_3cts <- function(this.snp,this.gene) {
  plot.mel <- make_PBSeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.pbs.mel,PBS_knownVar.mel) + theme(legend.position="none")
  plot.krt <- make_PBSeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.pbs.krt,PBS_knownVar.krt) + theme(legend.position="none")
  plot.frb <- make_PBSeQTL_plot_CPM_nocolor(this.snp,this.gene,CPM.pbs.frb,PBS_knownVar.frb) + theme(legend.position="none")

  p = plot_grid(plot.mel, plot.krt, plot.frb, labels=c("MEL", "KRT", "FRB"),ncol = 3, nrow = 1)
  p
  return(p)
}
# MAKE THE PLOTS!
this.df <- dict %>% dplyr::filter(padj<0.05) %>% dplyr::filter(abs(slope)<3)
plot_list = list()
options(scipen = 0,digits=3) # turn on scientific notation with 3 digits
for (i in 1:nrow(this.df)) {
  this.edQTL=this.df[i,"edQTL"] %>% as.character
  this.edSite=this.df[i,"edSite"] %>% as.character
  this.gene=edSite_gene_lookup[which(edSite_gene_lookup$V1==this.edSite),2][[1]][1]
  dist_edSite_gene=edSite_gene_lookup[which(edSite_gene_lookup$V1==this.edSite),3] %>% unlist() %>% as.character()
  plot_edQTL_qqnorm <- make_edQTL_plot_qqnorm(this.edQTL, this.edSite)
  plot_edQTL_ratio  <- make_edQTL_plot_ratio(this.edQTL, this.edSite)
  plot_edQTL_coverage_histogram <- make_edQTL_plot_coverage_histogram(this.edQTL, this.edSite)
  temp.plot <- make_PBSeQTL_plot_CPM_3cts(this.edQTL,this.gene) 
  this.title=paste0(this.gene," is ",dist_edSite_gene,"bp away from ",this.edSite)
  title <- ggdraw() + draw_label(this.title, fontface='bold', x=0.01, hjust=0)
  plot_eQTL_3cts_edQTL_edSiteGene <- plot_grid(title,temp.plot, ncol=1,rel_heights=c(0.1,1))
  
  tt <- ttheme_minimal(base_size = 5, core=list(fg_params=list(cex = 2.0)))
  mainidx <- which(dict$edSite==this.edSite & dict$edQTL==this.edQTL)
  tab=dict[mainidx,c("distance","slope","padj")]
  plot_table = plot_grid(tableGrob(tab, cols = colnames(tab), rows=NULL, theme = tt), ncol=1, nrow=1)
  
  p.top = plot_grid(plot_edQTL_qqnorm, plot_edQTL_ratio, ncol=2,nrow=1)
  p.left = plot_grid(plot_eQTL_3cts_edQTL_edSiteGene,plot_table, ncol=1,nrow=2)
  p.right = plot_edQTL_coverage_histogram
  p.bottom = plot_grid(p.left,p.right,ncol=2,nrow=1)
  p = plot_grid(p.top,p.bottom,ncol=1,nrow=2)
  plot_list[[i]] = p
}
pdf(paste0(dir,"/edQTL_plots_significant_results_padj_0.05.pdf"),width=15,height=8)
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}
dev.off()
options(scipen = 100) # turn back to normal format
cat("done plotting edQTL plot! \n")

# check edGenes
newtab <- left_join(this.df, edSite_gene_lookup, by=c("edSite"="V1"))
write.table(newtab$V2,"~/Downloads/nl/human/skin/eQTLs/edQTL/output/edGenes.txt",quote = F,col.names = F,row.names = F)

# check GWAS SNP overlap
gwas_snps=data.table::fread(paste0(dir,"/gwas_snps_and_LD_tags_long.txt"),header=F)
this.df %>% dplyr::filter(this.df$edQTL %in% gwas_snps$V1)
