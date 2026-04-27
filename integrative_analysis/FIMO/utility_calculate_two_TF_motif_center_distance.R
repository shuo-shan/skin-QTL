library(dplyr)
library(magrittr)
library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)

Dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/FIMO"
setwd(Dir)

# anchor to TF1 motif center, calculate how far away TF2 motif distance is relative to TF1 center
TF1="STAT1"
TF2="IRF1"
TF1_positive_motifCenter <- 
  data.table::fread(paste0(Dir,"/STAT1_STAT1HighCorrelationGenesEnhancers/fimo_topRanked_STAT1_motif_center.bed")) %>%
  set_colnames(c("TF1_center_chr","TF1_center_start","TF1_center_end","name","score","strand")) 
  
TF1_positive_motifFlank <- 
  data.table::fread(paste0(Dir,"/STAT1_STAT1HighCorrelationGenesEnhancers/fimo_topRanked_STAT1_motif_center_flanking_1kb.bed")) %>%
  set_colnames(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name","score","strand")) %>%
  dplyr::select(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name"))

TF1_positive_df <- left_join(TF1_positive_motifCenter, TF1_positive_motifFlank, by="name") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF1_center_chr","TF1_center_start","TF1_center_end","TF1_flank"))
# TF2
TF2_motifCenter <- 
  data.table::fread(paste0(Dir,"/IRF1_STAT1TopMatchingSTAT1HighCorrelationGenesEnhancers/fimo_topRanked_IRF1_motif_center.bed")) %>%
  set_colnames(c("TF2_center_chr","TF2_center_start","TF2_center_end","name","score","strand")) %>%
  tidyr::separate(., col=name, into=c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","TF2_start","TF2_end"),sep="_") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF2_center_chr","TF2_center_start","TF2_center_end","TF1_flank")) %>%
  left_join( . , TF1_positive_df, by="TF1_flank") %>%
  mutate(distance=TF1_center_start - TF2_center_start) %>%
  mutate(tag="STAT1HighCorrelationGenesEnhancers")


######### negative control
# anchor to TF1 motif center, calculate how far away TF2 motif distance is relative to TF1 center
TF1_negative_motifCenter <- 
  data.table::fread(paste0(Dir,"/STAT1_STAT1LowCorrelationGenesEnhancers/fimo_topRanked_STAT1_motif_center.bed")) %>%
  set_colnames(c("TF1_center_chr","TF1_center_start","TF1_center_end","name","score","strand")) 

TF1_negative_motifFlank <- 
  data.table::fread(paste0(Dir,"/STAT1_STAT1LowCorrelationGenesEnhancers/fimo_topRanked_STAT1_motif_center_flanking_1kb.bed")) %>%
  set_colnames(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name","score","strand")) %>%
  dplyr::select(c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","name"))

TF1_negative_df <- left_join(TF1_negative_motifCenter, TF1_negative_motifFlank, by="name") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF1_center_chr","TF1_center_start","TF1_center_end","TF1_flank"))

# TF2
TF2_motifCenter_negative <- 
  data.table::fread(paste0(Dir,"/IRF1_STAT1TopMatchingSTAT1LowCorrelationGenesEnhancers/fimo_topRanked_IRF1_motif_center.bed")) %>%
  set_colnames(c("TF2_center_chr","TF2_center_start","TF2_center_end","name","score","strand")) %>%
  tidyr::separate(., col=name, into=c("TF1_flank_chr","TF1_flank_start","TF1_flank_end","TF2_start","TF2_end"),sep="_") %>%
  mutate(TF1_flank=paste(TF1_flank_chr,TF1_flank_start,TF1_flank_end,sep="_")) %>%
  dplyr::select(c("TF2_center_chr","TF2_center_start","TF2_center_end","TF1_flank")) %>%
  left_join( . , TF1_negative_df, by="TF1_flank") %>%
  mutate(distance=TF1_center_start - TF2_center_start) %>%
  mutate(tag="STAT1LowCorrelationGenesEnhancers")

# plot!
df <- rbind(TF2_motifCenter, TF2_motifCenter_negative)
ggplot(df, aes(x=distance,fill=tag)) +
  geom_histogram(color="white",position="identity",alpha=0.7) +
  theme(legend.position="top") +
  xlab("motif center distance (bp)") +
  ggtitle(paste("motif center distance between",TF1,"and",TF2,sep=" ")) +
  facet_wrap(~tag)
