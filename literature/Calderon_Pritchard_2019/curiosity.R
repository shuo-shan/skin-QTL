library(dplyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(ggrepel)
library(ggfortify)
library(reactable)

# load data
rna1 <- read.table("/nl/umw_manuel_garber/human/skin/eQTLs/literature/Calderon_Pritchard_2019/metadata_RNAseq.txt",sep=",",header=TRUE,row.names=1)
rna2 <- read.table("/nl/umw_manuel_garber/human/skin/eQTLs/literature/Calderon_Pritchard_2019/rna.txt",sep="\t",header=F) %>% magrittr::set_colnames(c("GEO_Accession..exp.","spots"))
rna <- inner_join(rna2,rna1,by="GEO_Accession..exp.")
rna$spots <- rna$spots %>% gsub("M","",.) %>% as.numeric
rm(rna1,rna2)

atac1 <- read.table("/nl/umw_manuel_garber/human/skin/eQTLs/literature/Calderon_Pritchard_2019/metadata_ATACseq.txt",sep=",",header=TRUE,row.names=1)
atac2 <- read.table("/nl/umw_manuel_garber/human/skin/eQTLs/literature/Calderon_Pritchard_2019/atac.txt",sep="\t",header=F) %>% magrittr::set_colnames(c("GEO_Accession..exp.","spots"))
atac <- inner_join(atac2,atac1,by="GEO_Accession..exp.")
atac$spots <- atac$spots %>% gsub("M","",.) %>% as.numeric
rm(atac1,atac2)


# summary statistics of read number (in millions)
plot(hist(rna$spots))
summary(rna$spots)
summary(rna$AvgSpotLen)

plot(hist(atac$spots))
summary(atac$spots)
summary(atac$AvgSpotLen)
