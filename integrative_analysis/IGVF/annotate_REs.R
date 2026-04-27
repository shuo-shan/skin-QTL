library(dplyr)
library(magrittr)
library(tidyverse)
library(gplots)
library(ggfortify)
library(gridExtra)
library(grid)
library(cowplot)
library(ggpubr)
library(DESeq2)
library(ComplexHeatmap)

Dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis"

#res.atac.dds.frb = readRDS("~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_FRB.Rdata")
#res.atac.dds.krt = readRDS("~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_KRT.Rdata")
#res.atac.dds.mel = readRDS("~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_MEL.Rdata")
#res.chip.dds.frb = readRDS("~/Downloads/nl/human/skin/eQTLs/ChIP-seq/DESeq2/DESeq2_results_H3K27acChIPseq_peaks_1kb_FRB.Rdata")
#res.chip.dds.krt = readRDS("~/Downloads/nl/human/skin/eQTLs/ChIP-seq/DESeq2/DESeq2_results_H3K27acChIPseq_peaks_1kb_KRT.Rdata")
#res.chip.dds.mel = readRDS("~/Downloads/nl/human/skin/eQTLs/ChIP-seq/DESeq2/DESeq2_results_H3K27acChIPseq_peaks_1kb_MEL.Rdata")
#link_gene_promoter = readRDS(paste0(Dir,"/heatmap_RNA_round2/link_gene_promoter_within_300kbp_of_TSS_allcts.Rdata"))
#link_gene_enhancer = readRDS(paste0(Dir,"/heatmap_RNA_round2/link_gene_enhancer_within_300kbp_of_TSS_allcts.Rdata"))
FRB_cRE_dynamics = readRDS(paste0(Dir,"/heatmap_RNA_round2/FRB_gene_and_cRE_RNA_ATAC_and_H3K27ac_dynamics.Rdata"))
MEL_cRE_dynamics = readRDS(paste0(Dir,"/heatmap_RNA_round2/MEL_gene_and_cRE_RNA_ATAC_and_H3K27ac_dynamics.Rdata"))
KRT_cRE_dynamics = readRDS(paste0(Dir,"/heatmap_RNA_round2/KRT_gene_and_cRE_RNA_ATAC_and_H3K27ac_dynamics.Rdata"))

outDir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/IGVF"

# rename
FRB_cRE_dynamics <- FRB_cRE_dynamics %>% 
  dplyr::select(-c(gene_cluster, gene_in_cluster_order)) %>%
  rename_with(~ paste0( . , "_FRB"), -region)
KRT_cRE_dynamics <- KRT_cRE_dynamics %>% 
  dplyr::select(-c(gene_cluster, gene_in_cluster_order)) %>%
  rename_with(~ paste0( . , "_KRT"), -region)
MEL_cRE_dynamics <- MEL_cRE_dynamics %>% 
  dplyr::select(-c(gene_cluster, gene_in_cluster_order)) %>%
  rename_with(~ paste0( . , "_MEL"), -region)
cRE_dynamics_merged <- FRB_cRE_dynamics %>%
  full_join(KRT_cRE_dynamics, by = "region") %>%
  full_join(MEL_cRE_dynamics, by = "region")

# load ATAC coverage table
cov_table <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/merged_coverage_ATACseq_allcts_allregions_1kbp.bed")
cov_matrix <- cov_table %>% dplyr::select(-c("chr","start","end","name"))
rownames(cov_matrix) <- cov_table$name
metadata <- data.frame(sample=colnames(cov_matrix)) %>%
  tidyr::separate(., col="sample",sep="_",into=c("donor","celltype","condition"), remove=FALSE)

ggplot(cov_table, aes(x = F55_FRB_PBS)) +
  geom_density() +
  scale_x_log10() +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
  labs(title = "Density of Coverage per Region")

### regions that are open at baseline in each celltype in at least 2 donors
# Out of 426,620,
# 115,257 commonly open regions in baseline
# 180,324 cell type specific open regions in baseline
threshold <- 50
open_matrix <- (cov_matrix > threshold) * 1
donors_open_in_baseline <- data.frame(donors_open_in_FRB = rowSums(open_matrix[ , grep("FRB_PBS", colnames(open_matrix))]),
                                      donors_open_in_MEL = rowSums(open_matrix[ , grep("MEL_PBS", colnames(open_matrix))]),
                                      donors_open_in_KRT = rowSums(open_matrix[ , grep("KRT_PBS", colnames(open_matrix))]))
open_in_baseline <- ((donors_open_in_baseline >= 2) * 1) %>% as.data.frame()
regions_cts_open <- rownames(open_in_baseline)[which(rowSums(open_in_baseline) == 1)]
regions_pairshared_open <- rownames(open_in_baseline)[which(rowSums(open_in_baseline) == 2)]
regions_common_open <- rownames(open_in_baseline)[which(rowSums(open_in_baseline) == 3)]
regions_cts_open_KRT <- intersect(rownames(donors_open_in_baseline[which(donors_open_in_baseline$donors_open_in_KRT>=2),]),
                                  regions_cts_open)
regions_cts_open_MEL <- intersect(rownames(donors_open_in_baseline[which(donors_open_in_baseline$donors_open_in_MEL>=2),]),
                                  regions_cts_open)
regions_cts_open_FRB <- intersect(rownames(donors_open_in_baseline[which(donors_open_in_baseline$donors_open_in_FRB>=2),]),
                                  regions_cts_open)
write(regions_cts_open, paste0(outDir,"/regions_cts_open.txt"))
write(regions_common_open, paste0(outDir,"/regions_common_open.txt"))
write(regions_pairshared_open, paste0(outDir,"/regions_pairshared_open.txt"))
write(c(regions_cts_open,regions_pairshared_open,regions_common_open), paste0(outDir,"/regions_open_baseline.txt"))
### regions that are commonly / cts inducible
#res.atac.dds.frb = readRDS("~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_FRB.Rdata")
#res.atac.dds.krt = readRDS("~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_KRT.Rdata")
#res.atac.dds.mel = readRDS("~/Downloads/nl/human/skin/eQTLs/ATACseq_DolphinNext/DESeq2/DESeq2_results_ATACseq_peaks_1kb_MEL.Rdata")
dynamic_atac_matrix = data.frame(dynamic_FRB = res.atac.dds.frb$dynamic,
                                 dynamic_KRT = res.atac.dds.krt$dynamic,
                                 dynamic_MEL = res.atac.dds.mel$dynamic)
rownames(dynamic_atac_matrix) = res.atac.dds.frb$region
rm(res.atac.dds.frb, res.atac.dds.krt, res.atac.dds.mel)
induced_atac_matrix = (dynamic_atac_matrix == "open") * 1
induced_atac_matrix[is.na(induced_atac_matrix)] <- 0
regions_induced <- rownames(induced_atac_matrix)[which(rowSums(induced_atac_matrix) >= 1)]
regions_cts_induced <- rownames(induced_atac_matrix)[which(rowSums(induced_atac_matrix) == 1)]
regions_pairshared_induced <- rownames(induced_atac_matrix)[which(rowSums(induced_atac_matrix) == 2)]
regions_common_induced <- rownames(induced_atac_matrix)[which(rowSums(induced_atac_matrix) == 3)]

regions_cts_induced_KRT <- intersect(regions_cts_induced,
                                     rownames(dynamic_atac_matrix[which(dynamic_atac_matrix$dynamic_KRT=="open"),]))

regions_cts_induced_MEL <- intersect(regions_cts_induced,
                                     rownames(dynamic_atac_matrix[which(dynamic_atac_matrix$dynamic_MEL=="open"),]))
regions_cts_induced_FRB <- intersect(regions_cts_induced,
                                     rownames(dynamic_atac_matrix[which(dynamic_atac_matrix$dynamic_FRB=="open"),]))

write(regions_cts_induced , paste0(outDir,"/regions_cts_induced.txt"))
write(regions_common_induced, paste0(outDir,"/regions_common_induced.txt"))
write(regions_induced , paste0(outDir,"/regions_induced.txt"))
write(regions_pairshared_induced , paste0(outDir,"/regions_pairshared_induced.txt"))



