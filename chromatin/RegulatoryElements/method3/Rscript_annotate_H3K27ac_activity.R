library(dplyr)
library(tidyverse)
library(magrittr)
library(stringr)

Dir="~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3"
setwd(Dir)
load(paste0(Dir,'/myEnvironment_annotate_H3K27ac_activity_for_atac.RData'))
#save.image(paste0(Dir,"/myEnvironment_annotate_H3K27ac_activity_for_atac.RData"))

# load tables
atac <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_300bp_flanking_window.bed")
colnames(atac) <- c("chr","start","end","peak","score")

atac.in.frb.ifn <- data.table::fread(paste0(Dir,"/temp_atac_active_in_frb_ifn.txt"),header=F)
atac.in.frb.pbs <- data.table::fread(paste0(Dir,"/temp_atac_active_in_frb_pbs.txt"),header=F)
atac.in.krt.ifn <- data.table::fread(paste0(Dir,"/temp_atac_active_in_krt_ifn.txt"),header=F)
atac.in.krt.pbs <- data.table::fread(paste0(Dir,"/temp_atac_active_in_krt_pbs.txt"),header=F)
atac.in.mel.ifn <- data.table::fread(paste0(Dir,"/temp_atac_active_in_mel_ifn.txt"),header=F)
atac.in.mel.pbs <- data.table::fread(paste0(Dir,"/temp_atac_active_in_mel_pbs.txt"),header=F)

df <- data.frame(peak=atac$peak)
df$KRT_PBS <- lapply(df$peak, function(x) ifelse(x %in% atac.in.krt.pbs$V1, "active_KRT_PBS", "inactive_KRT_PBS")) %>% unlist()
df$KRT_IFN <- lapply(df$peak, function(x) ifelse(x %in% atac.in.krt.ifn$V1, "active_KRT_IFN", "inactive_KRT_IFN")) %>% unlist()
df$MEL_PBS <- lapply(df$peak, function(x) ifelse(x %in% atac.in.mel.pbs$V1, "active_MEL_PBS", "inactive_MEL_PBS")) %>% unlist()
df$MEL_IFN <- lapply(df$peak, function(x) ifelse(x %in% atac.in.mel.ifn$V1, "active_MEL_IFN", "inactive_MEL_IFN")) %>% unlist()
df$FRB_PBS <- lapply(df$peak, function(x) ifelse(x %in% atac.in.frb.pbs$V1, "active_FRB_PBS", "inactive_FRB_PBS")) %>% unlist()
df$FRB_IFN <- lapply(df$peak, function(x) ifelse(x %in% atac.in.frb.ifn$V1, "active_FRB_IFN", "inactive_FRB_IFN")) %>% unlist()

# count number of inactive peaks in the 6 conditions
df$inactive_count <- apply(df,1,function(x) sum(grepl("inactive",x)))

# remove atac peaks that are inactive in all 6 cts x conditions
df2 <- df[which(df$inactive_count != 6),]

# write to file
data.table::fwrite(df2, file=paste0(Dir,"/dictionary_atac_peaks_annotated_H3K27ac_activity.txt"), quote=FALSE, sep="\t",row.names=F, col.names=F)

