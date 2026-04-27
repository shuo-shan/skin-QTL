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

Dir="~/Downloads/nl/human/skin/eQTLs/edQTL/output"
setwd(Dir)
#load(paste0(Dir,"/"))

# load in editing level table
header <- read.table(paste0(Dir,"/temp.header")) %>% as.character()
edLevel <- read.table(paste0(Dir,"/sig.edSites.edLevels.txt"))
colnames(edLevel) <- header

# for every donor, calculate the editing ratio for each editing site
donors <- header[2:length(header)]
edRatio <- matrix(nrow=nrow(edLevel), ncol=length(donors))
for ( j in 1:length(donors)) {
  this.donor <- donors[j]
  this.donor.edLevels <- data.frame(edLevel=edLevel[,this.donor]) %>% set_rownames(edLevel$chrom)
  this.donor.edLevels$edReads <- apply(this.donor.edLevels, 2, function(x) gsub("/.*","", x)) %>% as.numeric()
  this.donor.edLevels$totalReads <- lapply(this.donor.edLevels$edLevel, function(x) gsub(".*/","", x)) %>% as.numeric()
  this.donor.edLevels$edRatio <- this.donor.edLevels$edReads / this.donor.edLevels$totalReads
  this.donor.edLevels$edRatio[is.nan(this.donor.edLevels$edRatio)] <- NA
  edRatio[,j] <- this.donor.edLevels$edRatio
}

edRatio <- as.data.frame(edRatio) %>% set_colnames(donors) %>% set_rownames(edLevel$chrom)

# sum up the editing ratio for each donor
edRatio.sum <- data.frame(donor=donors,
                          edRatioSum=colSums(edRatio, na.rm=TRUE))

# write to file
write.table(edRatio.sum, file=paste0(Dir,"/sum_of_editingRatio_across_19_convincing_edSites.txt"), quote=F, row.names=F, col.names=T)

