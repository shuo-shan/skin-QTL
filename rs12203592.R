# Rscript 06/12/2022 Crystal Shan

library(dplyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggpubr)

# 1. load genotype file
genotype=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/rs12203592_genotype_hg19.bed", header=TRUE)
genotype$ALT="T" # R automatically identified alternative allele, "T" as "true".

# 2. make metadata
metadata=genotype[7:50] %>% t() %>% as.data.frame() %>% rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
metadata$genotype=metadata$genotype %>%
  gsub("0/0","CC", . ) %>%
  gsub("0/1","CT", . ) %>%
  gsub("1/1","TT", . )

# 3. load RNA-seq TMM-normalized CPM data
cpm.mel=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/IRF4_cpm_mel.txt", header=TRUE) %>% t() %>% as.data.frame()
cpm.krt=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/IRF4_cpm_krt.txt", header=TRUE) %>% t() %>% as.data.frame()
cpm.frb=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/IRF4_cpm_frb.txt", header=TRUE) %>% t() %>% as.data.frame()

# 4. make the plots 
cpm=cpm.mel
ct="MEL"
rownames(cpm)=rownames(cpm) %>%
  gsub("_S.","", .) %>%
  gsub(paste0("_",ct),"", .)
cpm=cpm %>% rownames_to_column("sample") %>% separate( . , sample, c("donor","condition"), sep="_", remove=TRUE) 
df=left_join(cpm, metadata, by="donor") %>% drop_na(.)
df$IRF4=log2(df$IRF4+1)
df.paired=pivot_wider(df[,c("donor","genotype","IRF4","condition")], 
                              names_from="condition", values_from="IRF4")
ggplot(df.paired , aes(x=genotype, y=PBS, color=genotype)) +
  geom_boxplot(width=0.5) +
  geom_jitter(size=3, alpha=0.4,width=0.1) +
  ggtitle("IRF4 expression in PBS") +
  xlab("rs12203592 genotype") + ylab("log2 CPM") + 
  theme_bw() +
  theme(legend.position="none", 
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=14,face="bold"))

ggpaired(df.paired,cond1="PBS",cond2="IFN",facet.by="genotype",
         color="condition",palette="aaas",line.color = "gray", line.size = 0.4,
         title=paste0("IRF4 expression"),
         width=0,
         xlab="condition",ylab="log2 CPM")








