# explore total library sequencing reads for each biosample
library(dplyr)
library(tidyverse)

Dir="/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq"
setwd(Dir)
lib_info <- read.table(file="/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/lib_size_12162021.txt",header=T,sep="\t")

biosample_info = lib_info %>% 
  dplyr::group_by(biosample) %>%
  dplyr::summarize(biosample_reads = sum(total_reads_sample))

biosample_info = biosample_info %>% tidyr::separate(.,biosample,into=c("donor","celltype","condition"),sep="_",remove=FALSE)
