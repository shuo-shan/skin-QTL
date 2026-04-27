library(tidyverse)
library(ggplot2)
Dir="/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink/"
setwd(Dir)
save.image(file=paste0(Dir,'/rs706779_LD_explore.RData'))
#load(paste0(Dir,'/rs706779_LD_explore.RData')

# load data
ld.table <- read.table(paste0(Dir,"rs706779.600SNPs.flanking.LDr2.txt"),sep="\t", header=TRUE) 
dict <- read.table(paste0(Dir,"37donors.chr10.txt"),sep=" ", header=FALSE) %>% set_colnames(c("RS_number","chr","position"))

# compile
df <- inner_join(ld.table,dict,by="RS_number") %>% dplyr::select(-"chr") %>% set_colnames(c("RS_number","r2","position"))
df$dist <- df$position-6056861 # position of rs706779 is 6056861

# plot: visualize the LD r2 value 
p=ggplot(df,aes(x=dist,y=r2)) +
  geom_line() +
  xlab("distance to rs706779 (bp)")
p
