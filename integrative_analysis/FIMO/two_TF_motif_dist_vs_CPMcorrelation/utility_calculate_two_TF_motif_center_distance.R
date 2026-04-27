#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
TF1=args[1]
TF2=args[2]
Dir=paste0("/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO","/",TF1)
# TF1="GATA2"
# TF2="RUNX1"
# Dir=paste0("~/Downloads/nl/human/skin/eQTLs/integrative_analysis/FIMO","/",TF1)

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
