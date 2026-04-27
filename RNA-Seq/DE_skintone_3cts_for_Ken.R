### Note: [1] change the dir at line 23 to load the environment at line 27.
###       [2] line 29 - 70 is the code chunk that runs DESeq2 for skintone. 
###           You don't have to run it cuz the results are already stored in the 
###           res.frb, res.krt, and res.mel variables in the environment.
###       [3] to make the plot for a gene, change gene name at line 73, 
###           then run the entire rest of the script. if plot doesn't show up,
###           type "p" in the console to plot.

library(dplyr)
library(tidyverse)
library(DESeq2)
library(magrittr)
library(stringr)
library(gplots)
library(ggrepel)
library(ggfortify)
library(reactable)
library(gridExtra)
library(ggpubr)
library(plotly)
library(cowplot)

Dir="~/Downloads/nl/human/skin/eQTLs/RNA-Seq" # <---- change this to yours
setwd(Dir)

# load data and metadata
load(paste0(Dir,'/myEnvironment_allcts_skintoneDE_forKen.RData'))

### conduct DESeq2 analysis on skintone to look for DE genes
# melanocytes
temp.CPM <- CPM[which(rowSums(CPM[,])>=10),which(grepl("MEL_PBS",colnames(CPM)))]
df <- counts[rownames(temp.CPM),colnames(temp.CPM)]
df.metadata <- metadata[which(metadata$celltype=="MEL" & metadata$condition=="PBS"),]
dds <- DESeqDataSetFromMatrix(countData = df, colData = df.metadata, design = ~ batch + skintone)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>% as.data.frame()
summary(results(dds,contrast=c("skintone", "medium_to_dark", "white")))
res.mel <- results(dds,contrast=c("skintone", "medium_to_dark", "white")) %>% as.data.frame %>% rownames_to_column("gene")
res.up <- res.mel %>% filter(log2FoldChange > 1.5) %>% filter(padj < 0.05)
res.down <- res.mel %>% filter(log2FoldChange < -1.5) %>% filter(padj < 0.05)
skintone.DEG.mel <- c(res.up$gene,res.down$gene)

# keratinocytes
temp.CPM <- CPM[which(rowSums(CPM[,])>=10),which(grepl("KRT_PBS",colnames(CPM)))]
df <- counts[rownames(temp.CPM),colnames(temp.CPM)]
df.metadata <- metadata[which(metadata$celltype=="KRT" & metadata$condition=="PBS"),]
dds <- DESeqDataSetFromMatrix(countData = df, colData = df.metadata, design = ~ batch + skintone)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>% as.data.frame()
summary(results(dds,contrast=c("skintone", "medium_to_dark", "white")))
res.krt <- results(dds,contrast=c("skintone", "medium_to_dark", "white")) %>% as.data.frame %>% rownames_to_column("gene")
res.up <- res.krt %>% filter(log2FoldChange > 1.5) %>% filter(padj < 0.05)
res.down <- res.krt %>% filter(log2FoldChange < -1.5) %>% filter(padj < 0.05)
skintone.DEG.krt <- c(res.up$gene,res.down$gene)

# fibroblasts
temp.CPM <- CPM[which(rowSums(CPM[,])>=10),which(grepl("FRB_PBS",colnames(CPM)))]
df <- counts[rownames(temp.CPM),colnames(temp.CPM)]
df.metadata <- metadata[which(metadata$celltype=="FRB" & metadata$condition=="PBS"),]
dds <- DESeqDataSetFromMatrix(countData = df, colData = df.metadata, design = ~ batch + skintone)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE) %>% as.data.frame()
summary(results(dds,contrast=c("skintone", "medium_to_dark", "white")))
res.frb <- results(dds,contrast=c("skintone", "medium_to_dark", "white")) %>% as.data.frame %>% rownames_to_column("gene")
res.up <- res.frb %>% filter(log2FoldChange > 1.5) %>% filter(padj < 0.05)
res.down <- res.frb %>% filter(log2FoldChange < -1.5) %>% filter(padj < 0.05)
skintone.DEG.frb <- c(res.up$gene,res.down$gene)

########## plot CPM for all donors across all cell types for gene of interest #########
g="GUCA1A" ## <- change the name of the gene you want to plot here

# pick out CPM paired data for each donor
skintone.DEG <- unique(c(skintone.DEG.mel, skintone.DEG.krt, skintone.DEG.frb))
skintone.DEG.CPM <- CPM[skintone.DEG,]
this.CPM <- as.data.frame(CPM[g,]) %>% rownames_to_column("sample") %>% left_join( . , metadata, by="sample")
df.paired.this.long <- data.frame(condition=c(rep("PBS",ncol(CPM.pbs)), rep("IFN",ncol(CPM.ifn))),
                                  CPM=c(t(CPM.pbs[g,]),t(CPM.ifn[g,])),
                                  sample_for_splitting=c(colnames(CPM.pbs),colnames(CPM.ifn))) %>%
  mutate(biosample=sample_for_splitting) %>%
  tidyr::separate(., col=sample_for_splitting, into=c("donor","celltype","condition"),sep="_") %>%
  left_join( . , metadata[,c("biosample","skintone")], by="biosample")
df.paired.this.long$condition <- factor(df.paired.this.long$condition, levels=c("PBS","IFN"))
df.paired.this.long$skintone <- factor(df.paired.this.long$skintone, levels=c("white","medium_to_dark"))

# DESeq2 analysis result for all 3 cell types
df.DEseq2.this <- rbind(res.mel[which(res.mel$gene==g),] %>% mutate(celltype="MEL"),
                        res.krt[which(res.krt$gene==g),] %>% mutate(celltype="KRT"),
                        res.frb[which(res.frb$gene==g),] %>% mutate(celltype="FRB"))

# Plot melanocyte data
plt.MEL <- df.paired.this.long %>% 
  dplyr::filter(celltype=="MEL") %>%
  ggplot(.,aes(x=condition, y=CPM)) +
  geom_line(aes(group=donor)) +
  geom_point(aes(color=donor)) +
  labs(title=paste0(g,"\n",
                    "skintone DE analysis:","\n"," log2FC=", round(df.DEseq2.this[1,]$log2FoldChange,3),"\n"," p.adj=", 
                    format(df.DEseq2.this[1,]$padj,digits=3,scientific=T),"\n")) +
  theme(title =element_text(size=5, face='bold'), legend.position="none") +
  facet_wrap(~skintone)

# Plot keratinocyte data
plt.KRT <- df.paired.this.long %>% 
  dplyr::filter(celltype=="KRT") %>%
  ggplot(.,aes(x=condition, y=CPM)) +
  geom_line(aes(group=donor)) +
  geom_point(aes(color=donor)) +
  labs(title=paste0(g,"\n",
                    "skintone DE analysis:","\n"," log2FC=", round(df.DEseq2.this[2,]$log2FoldChange,3),"\n"," p.adj=", 
                    format(df.DEseq2.this[2,]$padj,digits=3,scientific=T),"\n")) +
  theme(title =element_text(size=5, face='bold'), legend.position="none") +
  facet_wrap(~skintone)

# Plot fibroblast data
plt.FRB <- df.paired.this.long %>% 
  dplyr::filter(celltype=="FRB") %>%
  ggplot(.,aes(x=condition, y=CPM)) +
  geom_line(aes(group=donor)) +
  geom_point(aes(color=donor)) +
  labs(title=paste0(g,"\n",
                    "skintone DE analysis:","\n"," log2FC=", round(df.DEseq2.this[3,]$log2FoldChange,3),"\n"," p.adj=", 
                    format(df.DEseq2.this[3,]$padj,digits=3,scientific=T),"\n")) +
  theme(title =element_text(size=5, face='bold'), legend.position="none") +
  facet_wrap(~skintone)

# integrate three plots side by side
p <- plot_grid(plt.MEL, plt.KRT, plt.FRB, ncol = 3, nrow = 1)
# show plot in Plots tab
p
