library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)
library(gridExtra)
library(AddOns)
#devtools::install_github("Yuqing66/AddOns")

# Yuqing inDrop data
srt1 <- readRDS("~/Dropbox (UMass Medical School)/skin_eQTL/RNA-seq/singleCellData/merge6_pub4.rds")
unique(srt1$CellType)
levels(unique(srt1$orig.ident))

srt1_sub <- srt1[,srt1$orig.ident %in% c("CB043","VB065","VB076","VB077","VB096","VB150","VB173")]
VlnPlot.ssc(srt1_sub, c("DHX58"), split.by = "CellType", group.by = "orig.ident",text_sizes = c(15, 12, 10, 12, 10, 8, 2.5),group.order = c("VB065","VB076","VB077","VB096"),colors = c("black","grey","blue","pink","red"))

marker.genes.list<-c("DHX58","STAT1","STAT3")
VlnPlot.compact(srt_sub, assay = "RNA", genes = unique(unlist(marker.genes.list)), group.by = "orig.ident", 
                split.by = c("CellType","gene"), split.scale = "free", 
                split.label.pos = c("left", "bottom"), split.label.textonly = T, 
                split.label.rotate = "both", axis.hide = "both", flip = F, legend.hide = F, log_scale = T,
                violin.linewidth = 0.05) 

# Erica 10x data
srt2 <- readRDS("~/Dropbox (UMass Medical School)/skin_eQTL/RNA-seq/singleCellData/Erica_10x/filtered.RDS")
unique(srt2$CellType)
