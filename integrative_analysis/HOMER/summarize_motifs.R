library(dplyr)
library(magrittr)
library(tidyverse)
library(gplots)
library(ggfortify)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
library(ComplexHeatmap)

Dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/HOMER"
setwd(Dir)

##############
df.promoter = data.table::fread(paste0(Dir,"/motifSig_promoter_compiled.txt"),header=F) %>% set_colnames(c("motif","annotation"))
df.promoter$annotation = df.promoter$annotation %>% gsub("_promoter","",.)
df.promoter.collapsed = df.promoter %>% distinct_all() %>% group_by(motif) %>% 
  summarize(annotation=paste(annotation, collapse=','))
temp = df.promoter %>% distinct_all() %>% group_by(motif) %>% mutate(count=n()) %>% 
  dplyr::select(c("motif","count")) %>% arrange(motif) %>%
  distinct_at("motif", .keep_all=TRUE)
df.promoter.collapsed = left_join(df.promoter.collapsed, temp, by="motif")
write.table(df.promoter.collapsed, file=paste0(Dir,"/motifSig_promoter_compiled_collapsed.txt"),
            quote=FALSE, sep="\t", row.names=F, col.names=F)


##############
df.enhancer = data.table::fread(paste0(Dir,"/motifSig_enhancer_compiled.txt"),header=F) %>% set_colnames(c("motif","annotation"))
df.enhancer$annotation = df.enhancer$annotation %>% gsub("_enhancer","",.)
df.enhancer.collapsed = df.enhancer %>% distinct_all() %>% group_by(motif) %>% 
  summarize(annotation=paste(annotation, collapse=','))
temp = df.enhancer %>% distinct_all() %>% group_by(motif) %>% mutate(count=n()) %>% 
  dplyr::select(c("motif","count")) %>% arrange(motif) %>%
  distinct_at("motif", .keep_all=TRUE)
df.enhancer.collapsed = left_join(df.enhancer.collapsed, temp, by="motif")
write.table(df.enhancer.collapsed, file=paste0(Dir,"/motifSig_enhancer_compiled_collapsed.txt"),
            quote=FALSE, sep="\t", row.names=F, col.names=F)

