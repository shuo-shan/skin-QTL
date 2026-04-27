library(dplyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(reactable)
library(DESeq2)
library(edgeR)

############# FUNCTIONS ############
create_table <- function(Dir, file_lst) {
  f=file_lst[1]
  sample=gsub("_sorted_esat\\.txt\\.gene\\.txt","",f)
  df <- read.table(paste0(Dir,"/",f), sep="\t", header=TRUE) %>% 
    tidyr::unite( .,"region",c("Symbol","chr","strand"),sep=";",remove=TRUE) %>% 
    mutate(sum=rowSums(select_if(., is.numeric))) %>%
    dplyr::select(c("region",sum)) %>% 
    magrittr::set_colnames(c("region",sample))
  
  if (length(file_lst) > 1) {
    for (f in file_lst[2:length(file_lst)]) {
      sample=gsub("_sorted_esat\\.txt\\.gene\\.txt","",f)
      df.this <- read.table(paste0(Dir,"/",f), sep="\t", header=TRUE) %>% 
        tidyr::unite( .,"region",c("Symbol","chr","strand"),sep=";",remove=TRUE) %>% 
        mutate(sum=rowSums(select_if(., is.numeric))) %>%
        dplyr::select(c("region",sum)) %>% 
        magrittr::set_colnames(c("region",sample))
      
      df <- dplyr::left_join(df,df.this,by="region")
    }
  }
  df <- df %>% tibble::column_to_rownames( . , "region")
  return(df)
  rm(f,sample,df,df.this)
}

# check which variable correlates with each PC the most.
best_correlation_with_pc <- function(pcs, metadata, pc.n) {
  total.pcc=c()
  pca.this=pcs[,pc.n]
  for (x in colnames(metadata)) {
    metadata.this <- as.numeric(as.factor(metadata[,x]))
    total.pcc=c(total.pcc,cor(pca.this,metadata.this))
  }
  total.df=data.frame(feature=colnames(metadata),pcc=total.pcc,pcc.magnitude=abs(total.pcc))
  return(total.df)
}
########################################################################################################

############# Load raw ESAT counts & metadata
### List of directories, change these for your own Dropbox paths
Dir.MEL <- "/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_01122022/analysis"
Dir.FRB <- "/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_01152022/analysis"
Dir.KRT <- "/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_01142022/analysis"
### Step 1. Change directory for each cell type
Dir <- Dir.MEL
### Step 2. Load data
counts.keep <- read.table(paste0(Dir,"/raw_counts.txt"), sep="\t", header=TRUE)
metadata <- read.table(paste0(Dir,"/metadata.txt"), sep="\t", header=TRUE)
########################################################################################################

############# TMM normalization 
dge <- DGEList(counts.keep)
dge <- calcNormFactors(dge)
dge$samples # check normalizaiton factors
CPM <- cpm(dge, log=FALSE)

metadata.pbs <- metadata %>% dplyr::filter(condition=="PBS")
metadata.ifn <- metadata %>% dplyr::filter(condition=="IFN")
# Getting normalized CPM matrix
CPM.pbs <- as.data.frame(CPM)[,as.character(metadata.pbs$sample)] %>% 
  magrittr::set_colnames(as.character(metadata.pbs$donor)) 
CPM.ifn <- as.data.frame(CPM)[,as.character(metadata.ifn$sample)] %>% 
  magrittr::set_colnames(as.character(metadata.ifn$donor)) 

# Getting log transformed CPM
logCPM <- log2(CPM+1) %>% as.data.frame() 
logCPM.pbs <- as.data.frame(logCPM)[,as.character(metadata.pbs$sample)] %>% 
  magrittr::set_colnames(as.character(metadata.pbs$donor)) 
logCPM.ifn <- as.data.frame(logCPM)[,as.character(metadata.ifn$sample)] %>% 
  magrittr::set_colnames(as.character(metadata.ifn$donor)) 
log2FC <- logCPM.ifn - logCPM.pbs

# pick genes of top row variance
mat.logCPM.lnsd <- apply(logCPM,MARGIN=1,FUN=sd)
mat.logCPM.lncv <- sqrt( (exp(mat.logCPM.lnsd) )^2 -1) #coefficient of variance
logCPM.top <- logCPM[order(mat.logCPM.lncv,decreasing=T),][1:2000,]
########################################################################################################

############# PCA
pca = prcomp(t(logCPM.top),centered= T,scale = T)
loadings=as.data.frame(pca$rotation)
# percent variance in each pc
plot(pca)
pcs = pca$x
pca.pctVar = pca$sdev^2/sum(pca$sdev^2)
# check which variable correlates with each PC the most. the number 1 corresponds to PC1. 
# Change it to 2,3,4,etc and run the command, and view the output "tot" to see which metadata
# factor correlates with each PC the most.
tot=best_correlation_with_pc(pcs,metadata,1)

# plot
# top PC information
d=data.frame(PC1=pca$x[,1],PC2=pca$x[,2], PC3=pca$x[,3], PC4=pca$x[,4], PC5=pca$x[,5],PC6=pca$x[,6],sample=rownames(pcs))
d=inner_join(d,metadata,by="sample") %>% magrittr::set_rownames(d$sample)
colnames(metadata)
#d$group_tags <- cut(d$Unique_Reads_Aligned_STAR, breaks=6)
d$group_tags <- cut(d$total_reads, 
                    breaks=c(min(d$total_reads)-1,
                             10e+06, 15e+06, 20e+06,
                             max(d$total_reads)))
d %>%
  mutate(libprep_seq_dates=paste(d$lib_prep_date,d$seq_date)) %>%
  ggplot(., aes(y=PC1,x=PC2,col=seq_cycle_num,shape=condition,label=biosample)) + 
  geom_point(size=3) + 
  geom_text_repel(show.legend = F,size=2) + 
  ylab(sprintf("PC1 (%2.1f%% variance)",100*pca.pctVar[1])) + 
  xlab(sprintf("PC2 (%2.1f%% variance)",100*pca.pctVar[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
########################################################################################################

############# DE analysis
df <- apply(counts.keep[rowSums(counts.keep)>=10,],2,round,0) %>% as.data.frame() # filter out lowly expressed genes first.

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = metadata,
                              design = ~ batch + condition)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>% as.data.frame()
summary(results(dds,contrast=c("condition", "IFN", "PBS")))

res <- results(dds,contrast=c("condition", "IFN", "PBS")) %>% as.data.frame %>% rownames_to_column("gene")
res.up <- res %>% filter(log2FoldChange > 1.5) %>% filter(padj < 0.05)
res.down <- res %>% filter(log2FoldChange < -1.5) %>% filter(padj < 0.05)
res.nochange <- res %>% filter(abs(log2FoldChange) <= 1.5) %>% filter(padj < 0.05) %>% filter(!is.na(padj))
########################################################################################################

############# Plot paired data in CPM  
g="IFIT2" # <----- change this for each gene of interest

# pick out TMM-normalized CPM paired data & run wilcox test
df.paired.this <- data.frame(PBS=as.data.frame(t(CPM.pbs[g,])),
                             IFN=as.data.frame(t(CPM.ifn[g,]))) %>% set_colnames(c("PBS","IFN"))
df.paired.this.test = wilcox.test(t(CPM.pbs[g,]),t(CPM.ifn[g,]))

# DESeq2 analysis result
df.DEseq2.this <- res %>% dplyr::filter(gene==g)

# Plot
ggpaired(df.paired.this, cond1="PBS", cond2="IFN", 
         color="condition",palette=c("grey","red"),
         title=paste0("TMM-normalized CPM of ",g,"\n",
                      "wilcox test: p=",format(df.paired.this.test$p.value, digits=3,scientific=T),"\n",
                      "DE analysis: log2FC=", round(df.DEseq2.this$log2FoldChange,3)," ; p.adj=", 
                      format(df.DEseq2.this$padj,digits=3,scientific=T),"\n",
                      "donors=",as.character(nrow(df.paired.this))),
         xlab="condition",ylab="TMM-normalized CPM")
########################################################################################################

############# Hierarchical clustering & Heatmap
df <- logCPM

# set of candidate genes for clustering
n=20 # top number of genes with most abs(log2FC)
candidate_genes <- res %>% as.data.frame %>%
  dplyr::filter(padj < 0.05) %>%                   # restrict FDR
  dplyr::filter(abs(log2FoldChange) >= 1.5) %>%    # restrict effect size
  dplyr::top_n(. , n, abs(log2FoldChange)) %>%     # restrict top FC response
  pull(gene) %>%
  unique()

# filter data matrix based on candidate genes
df <- df %>% dplyr::filter(rownames(.) %in% candidate_genes)

# perform k-means-clustering to determine k_row.
wss <- (nrow(scale(df))-1)*sum(apply(scale(df),2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scale(df),centers=i)$withinss)
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares") 

# hierarchical clustering to check number of genes that fall in clusters.
gene_hclust <- hclust(dist(df,method="euclidean"), method = "complete")
cutree(gene_hclust, k = 4) %>% 
  enframe() %>% 
  rename(gene=name,cluster=value) %>%
  tidyr::separate(.,gene,into=c("gene","response"),sep="_") %>%
  group_by(response,cluster) %>%
  summarise(n=n())

# heatmap
heatmaply(df, 
          k_row = 2,
          scale = "row",
          fontsize_row = 5,
          fontsize_col = 5,
          col_side_colors = metadata[,c("condition","batch")],
          main = "row-scaled log2 TMM-normalized CPM")
