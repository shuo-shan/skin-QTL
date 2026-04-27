library(tidyverse)
library(dplyr)
library(magrittr)

####### load tables
# beta-comparison model results
f1.krt <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/KRT/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE))
f1.mel <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/MEL/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE))
f1.frb <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/FRB/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE))
# allele-specific
f2=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/allele_specific/QTLs_1E-02/QTL_1E-02_annotated_AlleleSpecificMark.bed") %>%
  mutate(tag=paste(chr,start,end,ID,REF,ALT,sep="_"))
# TF_peak
f3=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/TF_peak_overlapping/QTL_1E-02/QTL_overlapping_TF_peaks_collapsed.bed")%>%
  mutate(tag=paste(chr,start,end,ID,REF,ALT,sep="_")) %>%
  select(-c(chr,start,end,ID,REF,ALT))
# TF motif
f4=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/TF_motif_overlapping/QTL_1E-02/QTL_1E-02_overlapping_TFmotif.bed") %>%
  mutate(tag=paste(chr,start,end,ID,REF,ALT,sep="_")) %>%
  select(-c(chr,start,end,ID,REF,ALT))
# cRE and dynamics
f5=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/promoter_enhancer_overlapping/QTLs_1E-02/QTL_overlapping_CRE.bed") %>%
  mutate(tag=paste(chr,start,end,ID,REF,ALT,sep="_")) %>%
  select(-c(chr,start,end,ID,REF,ALT))
# LD trait
tempf6=data.table::fread("~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results/QTLs_LDtrait_association_ALLpopulation_r20.1.txt")
tempf6.1 <- tempf6 %>% dplyr::filter(R2 >= 0.6) %>% dplyr::select(c(Query,GWAS_Trait)) %>% group_by(Query) %>% summarize(GWAS_Trait=paste(GWAS_Trait, collapse=','))
tempf6.2 <- tempf6 %>% dplyr::filter(R2 >= 0.6) %>% 
  arrange(desc(R2)) %>%
  mutate(RS_Number_annotated=paste0(RS_Number,"(R2=",round(R2,2),")")) %>%
  dplyr::select(c(Query,RS_Number_annotated)) %>% 
  group_by(Query) %>% 
  summarize(high_LD_GWAS_SNP_and_R2=paste(RS_Number_annotated, collapse=','))
f6 <- left_join(tempf6.1, tempf6.2, by="Query") %>% set_colnames(c("ID","GWAS_Trait","high_LD_GWAS_SNP"))
# ANNOVAR
f7=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/ANNOVAR/QTLs_1E-02/QTL.hg38_multianno.txt") %>%
  mutate(start=Start-1) %>%
  mutate(tag=paste(Chr,start,End,avsnp147,Ref,Alt,sep="_")) %>%
  select(-c(Chr,start,Start,End,avsnp147,Ref,Alt))

# TSS
f8=data.table::fread("~/Downloads/nl/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_genelevel_transcription_start_sites.bed",header=F) %>%
  set_colnames(c("gene_chr","gene_start","gene_end","gene_strand","gene"))

####### join tables
res23 <- left_join(f2,f3,by="tag")
res234 <- left_join(res23, f4, by="tag")
res2345 <- left_join(res234, f5, by="tag")
res23456 <- left_join(res2345, f6, by="ID")
res234567 <- left_join(res23456, f7, by="tag")

res234567.MEL <- res234567 %>% dplyr::select(-contains("KRT")) %>% dplyr::select(-contains("FRB"))
res234567.FRB <- res234567 %>% dplyr::select(-contains("KRT")) %>% dplyr::select(-contains("MEL"))
res234567.KRT <- res234567 %>% dplyr::select(-contains("MEL")) %>% dplyr::select(-contains("FRB"))

rm(f2,f3,f4,f5,tempf6,tempf6.1,tempf6.2,f6,f7,res23,res234,res2345,res23456,res234567)

####### load modeling result
resF.mel=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/MEL_QTL_modeling_result.txt")
bigtable.mel <- resF.mel %>% mutate(SNPGENE=paste0(QTL,"_",gene)) %>%
  left_join( . , f1.mel, by=c("SNPGENE"="tag")) %>%
  left_join( . , res234567.MEL, by=c("QTL"="ID")) %>%
  left_join( . , f8, by="gene") %>%
  mutate(tss_pos = ifelse(gene_strand == "+", gene_start, gene_end)) %>%
  mutate(dist_to_tss = abs(start - tss_pos)) %>%
  dplyr::select(-c(SNPGENE,tag,SNP,GENE,tss_pos,gene_chr,gene_start,gene_end,gene_strand))
bigtable.mel$z.betaComp <- as.numeric(bigtable.mel$z.betaComp)
bigtable.mel$p.betaCompPnorm <- as.numeric(bigtable.mel$p.betaCompPnorm)
bigtable.mel$p.betaComp10KPermut <- as.numeric(bigtable.mel$p.betaComp10KPermut)
bigtable.mel$reQTL_pval <- as.numeric(bigtable.mel$reQTL_pval)
bigtable.mel$reQTL_beta <- as.numeric(bigtable.mel$reQTL_beta)
bigtable.mel$reQTL_se <- as.numeric(bigtable.mel$reQTL_se)
bigtable.mel$PBSeQTL_pval <- as.numeric(bigtable.mel$PBSeQTL_pval)
bigtable.mel$PBSeQTL_beta <- as.numeric(bigtable.mel$PBSeQTL_beta)
bigtable.mel$PBSeQTL_se <- as.numeric(bigtable.mel$PBSeQTL_se)
bigtable.mel$IFNeQTL_pval <- as.numeric(bigtable.mel$IFNeQTL_pval)
bigtable.mel$IFNeQTL_beta <- as.numeric(bigtable.mel$IFNeQTL_beta)
bigtable.mel$IFNeQTL_se <- as.numeric(bigtable.mel$IFNeQTL_se)
saveRDS(bigtable.mel, "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_mel.rds")
data.table::fwrite(bigtable.mel, "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_MEL.txt",
                   quote=F, sep="\t", row.names=F, col.names=T)
bigtable.mel %>%
  dplyr::filter((!is.na(reQTL_pval) & reQTL_pval < 0.00001) | 
                  (!is.na(p.betaComp10KPermut) & p.betaComp10KPermut < 0.001)) %>%
  data.table::fwrite( . , "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_MEL_reQTL.txt",
                      quote=F, sep="\t", row.names=F, col.names=T)

rm(f1.mel,res234567.MEL)

resF.krt=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/KRT_QTL_modeling_result.txt")
bigtable.krt <- resF.krt %>% mutate(SNPGENE=paste0(QTL,"_",gene)) %>%
  left_join( . , f1.krt, by=c("SNPGENE"="tag")) %>%
  left_join( . , res234567.KRT, by=c("QTL"="ID")) %>%
  left_join( . , f8, by="gene") %>%
  mutate(tss_pos = ifelse(gene_strand == "+", gene_start, gene_end)) %>%
  mutate(dist_to_tss = abs(start - tss_pos)) %>%
  dplyr::select(-c(SNPGENE,tag,SNP,GENE,tss_pos,gene_chr,gene_start,gene_end,gene_strand)) 
bigtable.krt$z.betaComp <- as.numeric(bigtable.krt$z.betaComp)
bigtable.krt$p.betaCompPnorm <- as.numeric(bigtable.krt$p.betaCompPnorm)
bigtable.krt$p.betaComp10KPermut <- as.numeric(bigtable.krt$p.betaComp10KPermut)
bigtable.krt$reQTL_pval <- as.numeric(bigtable.krt$reQTL_pval)
bigtable.krt$reQTL_beta <- as.numeric(bigtable.krt$reQTL_beta)
bigtable.krt$reQTL_se <- as.numeric(bigtable.krt$reQTL_se)
bigtable.krt$PBSeQTL_pval <- as.numeric(bigtable.krt$PBSeQTL_pval)
bigtable.krt$PBSeQTL_beta <- as.numeric(bigtable.krt$PBSeQTL_beta)
bigtable.krt$PBSeQTL_se <- as.numeric(bigtable.krt$PBSeQTL_se)
bigtable.krt$IFNeQTL_pval <- as.numeric(bigtable.krt$IFNeQTL_pval)
bigtable.krt$IFNeQTL_beta <- as.numeric(bigtable.krt$IFNeQTL_beta)
bigtable.krt$IFNeQTL_se <- as.numeric(bigtable.krt$IFNeQTL_se)
saveRDS(bigtable.krt, "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_krt.rds")
#data.table::fwrite(bigtable.krt, "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_KRT.txt",
#                   quote=F, sep="\t", row.names=F, col.names=T)

bigtable.krt %>%
  dplyr::filter((!is.na(reQTL_pval) & reQTL_pval < 0.00001) | 
                  (!is.na(p.betaComp10KPermut) & p.betaComp10KPermut < 0.001)) %>%
  data.table::fwrite( . , "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_KRT_reQTL.txt",
                      quote=F, sep="\t", row.names=F, col.names=T)

rm(f1.krt,res234567.KRT)

resF.frb=data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/FRB_QTL_modeling_result.txt")
bigtable.frb <- resF.frb %>% mutate(SNPGENE=paste0(QTL,"_",gene)) %>%
  left_join( . , f1.frb, by=c("SNPGENE"="tag")) %>%
  left_join( . , res234567.FRB, by=c("QTL"="ID")) %>%
  left_join( . , f8, by="gene") %>%
  mutate(tss_pos = ifelse(gene_strand == "+", gene_start, gene_end)) %>%
  mutate(dist_to_tss = abs(start - tss_pos)) %>%
  dplyr::select(-c(SNPGENE,tag,SNP,GENE,tss_pos,gene_chr,gene_start,gene_end,gene_strand))
bigtable.frb$z.betaComp <- as.numeric(bigtable.frb$z.betaComp)
bigtable.frb$p.betaCompPnorm <- as.numeric(bigtable.frb$p.betaCompPnorm)
bigtable.frb$p.betaComp10KPermut <- as.numeric(bigtable.frb$p.betaComp10KPermut)
bigtable.frb$reQTL_pval <- as.numeric(bigtable.frb$reQTL_pval)
bigtable.frb$reQTL_beta <- as.numeric(bigtable.frb$reQTL_beta)
bigtable.frb$reQTL_se <- as.numeric(bigtable.frb$reQTL_se)
bigtable.frb$PBSeQTL_pval <- as.numeric(bigtable.frb$PBSeQTL_pval)
bigtable.frb$PBSeQTL_beta <- as.numeric(bigtable.frb$PBSeQTL_beta)
bigtable.frb$PBSeQTL_se <- as.numeric(bigtable.frb$PBSeQTL_se)
bigtable.frb$IFNeQTL_pval <- as.numeric(bigtable.frb$IFNeQTL_pval)
bigtable.frb$IFNeQTL_beta <- as.numeric(bigtable.frb$IFNeQTL_beta)
bigtable.frb$IFNeQTL_se <- as.numeric(bigtable.frb$IFNeQTL_se)
saveRDS(bigtable.frb, "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_frb.rds")
#data.table::fwrite(bigtable.frb, "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_FRB.txt",
#                   quote=F, sep="\t", row.names=F, col.names=T)
bigtable.frb %>%
  dplyr::filter((!is.na(reQTL_pval) & reQTL_pval < 0.00001) | 
                  (!is.na(p.betaComp10KPermut) & p.betaComp10KPermut < 0.001)) %>%
  data.table::fwrite( . , "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_FRB_reQTL.txt",
                      quote=F, sep="\t", row.names=F, col.names=T)
rm(f1.frb,res234567.FRB)
rm(resF.frb, resF.krt, resF.mel)