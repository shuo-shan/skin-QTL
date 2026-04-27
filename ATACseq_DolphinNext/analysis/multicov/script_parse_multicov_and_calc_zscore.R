#!/usr/bin/env Rscript
# goal: parse coverage data of ATACseq signal on cross-celltype, cross-condition common set of peaks for each group of interest
# goal: calculate zscore for each sample
# written by Crystal Shan 02/2022
cat("heya!!!! =D, loading libraries now... \n")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))

df=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/multicov/multicov.merged.txt",header=TRUE,sep="\t")
atac.dir="/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/analysis/multicov"
df.info=df %>% dplyr::select(c("chr","start","end")) %>% tidyr::unite(.,"peak",c("chr","start","end"),sep="_",remove=FALSE)
df.info=df.info[,c("chr","start","end","peak")]
df.krt_PBS=df %>% dplyr::select(contains(c("F25K_PBS_S2","F49K_PBS_S1","F55K_PBS_S1"))) %>% set_colnames(c("F25","F49","F55")) %>% cbind(df.info,.)
df.frb_PBS=df %>% dplyr::select(contains(c("F25F_PBS_S2","F49F_PBS_S1","F55F_PBS_S1"))) %>% set_colnames(c("F25","F49","F55")) %>% cbind(df.info,.)
df.mel_PBS=df %>% dplyr::select(contains(c("F25M_PBS_S2","F49M_PBS_S1","F55M_PBS_S2"))) %>% set_colnames(c("F25","F49","F55")) %>% cbind(df.info,.)

df.krt_PBS.zscore=apply(df.krt_PBS[,c("F25","F49","F55")],2,scale) %>% as.data.frame() %>% cbind(df.info,.)
df.frb_PBS.zscore=apply(df.frb_PBS[,c("F25","F49","F55")],2,scale) %>% as.data.frame() %>% cbind(df.info,.)
df.mel_PBS.zscore=apply(df.mel_PBS[,c("F25","F49","F55")],2,scale) %>% as.data.frame() %>% cbind(df.info,.)

write.table(df.krt_PBS,paste0(atac.dir,"/atac_cov_KRT_PBS.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(df.frb_PBS,paste0(atac.dir,"/atac_cov_FRB_PBS.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(df.mel_PBS,paste0(atac.dir,"/atac_cov_MEL_PBS.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(df.krt_PBS.zscore,paste0(atac.dir,"/atac_zscore_KRT_PBS.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(df.frb_PBS.zscore,paste0(atac.dir,"/atac_zscore_FRB_PBS.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(df.mel_PBS.zscore,paste0(atac.dir,"/atac_zscore_MEL_PBS.bed"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
