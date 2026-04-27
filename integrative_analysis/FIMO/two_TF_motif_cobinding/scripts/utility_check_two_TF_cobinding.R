#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
TF1=args[1]
TF2=args[2]
Dir=paste0("/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_cobinding/results","/",TF1)
TF1="GATA2"
TF2="RUNX1"
Dir=paste0("~/Downloads/nl/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_cobinding/results","/",TF1)

######### central hypothesis ##########
# TF might regulate highly-TF-correlating genes, compared to non-correlated genes. 
# when TF1 highly correlate with genes, it suggests that TF1 might regulate these genes,
# compared to the poorly correlated genes.
# And for TF2 that correlates the same genes as TF1, it might regulate the same genes as TF1.
# Therefore if I take these TF-correlating genes, 
# I should see both TF1 and TF2 binding their promoters, compared to TF-noncorrelating genes.

# Hypothesis #1: There is an association between gene correlation with TF1 and 
#                the presence of TF1 motifs in hiCor.promoters.
# Rows: Highly correlating genes vs. Poorly correlating genes.
# Columns: Promoters with TF1 motif vs. Promoters without TF1 motif.


# Hypothesis #2: There is an association between the presence of TF1 motifs and the
#                presence of TF2 motifs in promoters of highly correlating genes.
# Rows: Promoters with TF1 motif vs. Promoters without TF1 motif.
# Columns: Promoters with TF2 motif vs. Promoters without TF2 motif.


######### TF1-highly correlating genes ###########
### How many promoters are there?
this.dir <- paste0(Dir,"/",TF1,"_",TF1,"HighCorrelationGenesPromoters")
hiCor.promoters.total <- data.table::fread(paste0(this.dir,"/regions.bed"))
hiCor.promoters.total.N <- length(unique(hiCor.promoters.total$V4))
### How many promoters have TF1 motif?
hiCor.promoters.hasTF1 <- data.table::fread(paste0(this.dir,"/fimo_sig_",TF1,".bed"))
hiCor.promoters.hasTF1.N <- length(unique(hiCor.promoters.hasTF1$V4))
### How many promoters don't have TF1 motif?
hiCor.promoters.noTF1.N <- hiCor.promoters.total.N - hiCor.promoters.hasTF1.N 
### How many promoters that have TF1 motif also have TF2 motif?
this.dir <- paste0(Dir,"/",TF2,"/",TF1,"TopMatching",TF1,"HighCorrelationGenesPromoters")
tryCatch({
  hiCor.promoters.hasTF1TF2 <- data.table::fread(paste0(this.dir,"/fimo_sig_",TF2,".bed"))
  hiCor.promoters.hasTF1TF2.N <- nrow(hiCor.promoters.hasTF1TF2)
}, error = function(e) {
  hiCor.promoters.hasTF1TF2.N <- 0
})

### How many promoters that have TF1 motif do not have TF2 motif?
hiCor.promoters.hasTF1noTF2.N <- hiCor.promoters.hasTF1.N - hiCor.promoters.hasTF1TF2.N


######### TF1-poorly correlating genes ###########
### How many promoters are there?
this.dir <- paste0(Dir,"/",TF1,"_",TF1,"LowCorrelationGenesPromoters")
lowCor.promoters.total <- data.table::fread(paste0(this.dir,"/regions.bed"))
lowCor.promoters.total.N <- length(unique(lowCor.promoters.total$V4))
### How many promoters have TF1 motif?
tryCatch ({
  lowCor.promoters.hasTF1 <- data.table::fread(paste0(this.dir,"/fimo_sig_",TF1,".bed"))
  lowCor.promoters.hasTF1.N <- length(unique(lowCor.promoters.hasTF1$V4))
}, error = function(e) {
  lowCor.promoters.hasTF1.N <- 0
})

### How many promoters don't have TF1 motif?
lowCor.promoters.noTF1.N <- lowCor.promoters.total.N - lowCor.promoters.hasTF1.N 
### How many promoters that have TF1 motif also have TF2 motif?
this.dir <- paste0(Dir,"/",TF2,"/",TF1,"TopMatching",TF1,"LowCorrelationGenesPromoters")
tryCatch ({
  lowCor.promoters.hasTF1TF2 <- data.table::fread(paste0(this.dir,"/fimo_sig_",TF2,".bed"))
  lowCor.promoters.hasTF1TF2.N <- length(unique(lowCor.promoters.hasTF1TF2$V4))
}, error = function(e) {
  lowCor.promoters.hasTF1TF2.N <- 0
})
### How many promoters that have TF1 motif do not have TF2 motif?
lowCor.promoters.hasTF1noTF2.N <- lowCor.promoters.hasTF1.N - lowCor.promoters.hasTF1TF2.N

######### HIGH-Cor vs. LOW-Cor ######### 
### construct a contingency table
lowCor.promoters.total.N
lowCor.promoters.hasTF1.N
lowCor.promoters.hasTF1TF2.N
lowCor.promoters.hasTF1noTF2.N 

hiCor.promoters.total.N
hiCor.promoters.hasTF1.N
hiCor.promoters.hasTF1TF2.N
hiCor.promoters.hasTF1noTF2.N 

# compare the prescence of TF2 among promoters that have TF1
contingency_table <- matrix(c(lowCor.promoters.hasTF1TF2.N, hiCor.promoters.hasTF1TF2.N,  # Presence of TF2 in lowCor and hiCor
                              lowCor.promoters.hasTF1noTF2.N , hiCor.promoters.hasTF1noTF2.N ), # Absence of TF2 in lowCor and hiCor
                            nrow = 2, byrow = TRUE,
                            dimnames = list(c("TF2_Present", "TF2_Absent"),
                                            c("lowCor", "hiCor")))

# Perform Fisher's exact test
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)


# anchor to TF1 motif center, calculate how far away TF2 motif distance is relative to TF1 center
TF1_positive_motifCenter <-
  data.table::fread(paste0(Dir,"/",TF1,"_",TF1,"HighCorrelationGenesEnhancers/fimo_topRanked_",TF1,"_motif_center.bed")) %>%
  set_colnames(c("TF1_center_chr","TF1_center_start","TF1_center_end","name","score","strand"))

TF1_positive_motifFlank <-
  data.table::fread(paste0(Dir,"/",TF1,"_",TF1,"HighCorrelationGenesEnhancers/fimo_topRanked_",TF1,"_motif_center_flanking_1kb.bed")) %>%
  set_colnames(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name","score","strand")) %>%
  dplyr::select(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name"))

TF1_positive_df <- left_join(TF1_positive_motifCenter, TF1_positive_motifFlank, by="name") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF1_center_chr","TF1_center_start","TF1_center_end","TF1_flank"))

# TF1 motif center in poorly correlated gene enhancers
TF1_negative_motifCenter <-
  data.table::fread(paste0(Dir,"/",TF1,"_",TF1,"LowCorrelationGenesEnhancers/fimo_topRanked_",TF1,"_motif_center.bed")) %>%
  set_colnames(c("TF1_center_chr","TF1_center_start","TF1_center_end","name","score","strand"))

TF1_negative_motifFlank <-
  data.table::fread(paste0(Dir,"/",TF1,"_",TF1,"LowCorrelationGenesEnhancers/fimo_topRanked_",TF1,"_motif_center_flanking_1kb.bed")) %>%
  set_colnames(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name","score","strand")) %>%
  dplyr::select(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name"))

TF1_negative_df <- left_join(TF1_negative_motifCenter, TF1_negative_motifFlank, by="name") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF1_center_chr","TF1_center_start","TF1_center_end","TF1_flank"))

# TF2
Dir.TF2 <- paste0(Dir,"/",TF2)
# TF2
TF2_motifCenter_highCorr <-
  data.table::fread(paste0(Dir.TF2,"/",list.files(Dir.TF2, pattern="High"))) %>%
  set_colnames(c("TF2_center_chr","TF2_center_start","TF2_center_end","name","score","strand")) %>%
  tidyr::separate(., col=name, into=c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","TF2_start","TF2_end"),sep="_") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF2_center_chr","TF2_center_start","TF2_center_end","TF1_flank")) %>%
  left_join( . , TF1_positive_df, by="TF1_flank") %>%
  mutate(distance=TF1_center_start - TF2_center_start) %>%
  mutate(tag=paste0(TF1,"HighCorrelationGenesEnhancers"))

TF2_motifCenter_lowCorr <-
  data.table::fread(paste0(Dir.TF2,"/",list.files(Dir.TF2, pattern="Low"))) %>%
  set_colnames(c("TF2_center_chr","TF2_center_start","TF2_center_end","name","score","strand")) %>%
  tidyr::separate(., col=name, into=c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","TF2_start","TF2_end"),sep="_") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF2_center_chr","TF2_center_start","TF2_center_end","TF1_flank")) %>%
  left_join( . , TF1_negative_df, by="TF1_flank") %>%
  mutate(distance=TF1_center_start - TF2_center_start) %>%
  mutate(tag=paste0(TF1,"LowCorrelationGenesEnhancers"))

# plot
outDir <- Dir
pdf(paste0(outDir,"/histograms_",TF1,"_",TF2,".pdf"))
df <- rbind(TF2_motifCenter_highCorr , TF2_motifCenter_lowCorr )
ggplot(df, aes(x=distance,fill=tag)) +
  geom_histogram(color="white",position="identity",alpha=0.7) +
  theme(legend.position="top") +
  xlab("motif center distance (bp)") +
  ggtitle(paste("motif center distance between",TF1,"and",TF2,sep=" ")) +
  facet_wrap(~tag)
dev.off()

# wilcoxon Rank Sum test
x <- TF2_motifCenter_highCorr$distance %>% na.omit()
y <- TF2_motifCenter_lowCorr$distance %>% na.omit()
res <- data.frame(x=TF1, y=TF2, pval=unlist(wilcox.test(x,y))[2])
data.table::fwrite(res, file=paste0(outDir,"/wilcoxon_test_pval_",TF1,"_",TF2,".txt"),quote=F,col.names=F,sep="\t")
