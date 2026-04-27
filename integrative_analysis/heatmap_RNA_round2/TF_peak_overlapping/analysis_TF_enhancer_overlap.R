# written by shuo.shan@umassmed.edu 05/2024
library(tidyverse)
library(magrittr)
library(ComplexHeatmap)

dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/TF_peak_overlapping"
link_gene_enhancer = readRDS(paste0(dir,"/../link_gene_enhancer_within_300kbp_of_TSS_allcts.Rdata"))
load(paste0(dir,"/../myEnvironment_heatmap_RNAseq_DEgenes_14kmm-hc_04182023_padj0.05_log2FC1.5_avgCPM10.RData"))

#### compile enhancer TF overlapping table ####
# enhancer was defined to be linked to both the nearest upstream and downstream expressed gene for that celltype
f = data.table::fread(paste0(dir,"/enhancer_overlapping_TF_peaks_cleaned.txt"))
# aggregate data
aggregated_data <- f %>%
  group_by(gene, celltype, TF) %>%
  summarise(overlap_count = n(), .groups = 'drop')
# Pivot the data to create a matrix format
tf_matrix <- aggregated_data %>%
  pivot_wider(names_from = TF, values_from = overlap_count, values_fill = list(overlap_count = 0))

#### Heatmap of TF enhancer Overlaps: #1 cell-type-agnostic enhancers ####
# Specify split instructions
aggregated_data_ctagnostic <- f %>%
  group_by(gene, TF) %>%
  summarise(overlap_count = n(), .groups = 'drop')
tf_matrix_ctagnostic <- aggregated_data_ctagnostic %>%
  pivot_wider(names_from = TF, values_from = overlap_count, values_fill = list(overlap_count = 0)) %>%
  dplyr::filter(gene %in% intersect(gene, rownames(km_res.rna.new))) %>%
  column_to_rownames("gene") %>%
  as.matrix()

rowidx <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% unique(rownames( tf_matrix_ctagnostic ))) %>% 
  dplyr::select(c(newclass, within_cluster_order))
split_instructions.rna <- as.character( rowidx$newclass )

heatmap <- Heatmap(tf_matrix_ctagnostic, 
                   name = "Overlap Count", 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE,
                   show_row_names = FALSE, 
                   show_column_names = TRUE,
                   column_names_gp = grid::gpar(fontsize = 8),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_rot = 45,
                   row_gap = unit(1, "mm"),
                   border = TRUE,
                   row_title = "Genes", 
                   column_title = "Number of Overlap of TF and enhancers of DE genes, cell-type-agnostic",
                   column_dend_height = unit(1, "cm"),
                   row_split = split_instructions.rna)  # Use row_split to divide rows based on split_instructions

pdf(paste0(dir,"/heatmap_enhancer_TF_overlap_celltype_agnostic.pdf"),width=24,height=30)
draw(heatmap) # Adjust padding to ensure labels fit
dev.off()


# binary data matrix (as long as TF overlap any enhancer of the gene, count as 1)
#tf_matrix_ctagnostic[tf_matrix_ctagnostic>0] <- 1
#rowidx <- km_res.rna.new %>% 
#  dplyr::filter(rownames(.) %in% unique(rownames( tf_matrix_ctagnostic ))) %>% 
#  dplyr::select(c(newclass, within_cluster_order))
#split_instructions.rna <- as.character( rowidx$newclass )
heatmap <- Heatmap(tf_matrix_ctagnostic, 
                   name = "Overlap Count", 
                   cluster_rows = TRUE, 
                   cluster_columns = TRUE,
                   show_row_names = FALSE, 
                   show_column_names = TRUE,
                   column_names_gp = grid::gpar(fontsize = 8),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_rot = 45,
                   row_gap = unit(1, "mm"),
                   border = TRUE,
                   row_title = "Genes", 
                   column_title = "Number of Overlap of TF and enhancers of DE genes, cell-type-agnostic",
                   column_dend_height = unit(1, "cm"))  # Use row_split to divide rows based on split_instructions

pdf(paste0(dir,"/heatmap_enhancer_TF_overlap_celltype_agnostic_binary_rowclustered.pdf"),width=24,height=30)
draw(heatmap) # Adjust padding to ensure labels fit
dev.off()


#### Heatmap of TF enhancer Overlaps: #2 MEL enhancers ####
# Specify split instructions
aggregated_data_MEL <- f %>%
  group_by(gene, celltype, TF) %>%
  summarise(overlap_count = n(), .groups = 'drop') %>%
  dplyr::filter(celltype=="MEL") %>%
  dplyr::select(-celltype)

tf_matrix_MEL <- aggregated_data_MEL %>%
  pivot_wider(names_from = TF, values_from = overlap_count, values_fill = list(overlap_count = 0)) %>%
  dplyr::filter(gene %in% intersect(gene, rownames(km_res.rna.new))) %>%
  column_to_rownames("gene") %>%
  as.matrix()

rowidx <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% unique(rownames( tf_matrix_MEL ))) %>% 
  dplyr::select(c(newclass, within_cluster_order))
split_instructions.rna <- as.character( rowidx$newclass )

heatmap <- Heatmap(tf_matrix_MEL, 
                   name = "Overlap Count", 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE,
                   show_row_names = FALSE, 
                   show_column_names = TRUE,
                   column_names_gp = grid::gpar(fontsize = 8),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_rot = 45,
                   row_gap = unit(1, "mm"),
                   border = TRUE,
                   row_title = "Genes", 
                   column_title = "Number of Overlap of TF and enhancers of DE genes, MEL enhancers",
                   column_dend_height = unit(1, "cm"),
                   row_split = split_instructions.rna)  # Use row_split to divide rows based on split_instructions

# Draw the heatmap
pdf(paste0(dir,"/heatmap_enhancer_TF_overlap_MEL.pdf"),width=24,height=30)
draw(heatmap) # Adjust padding to ensure labels fit
dev.off()

# draw binary
tf_matrix_MEL[tf_matrix_MEL>0] <- 1
rowidx <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% unique(rownames( tf_matrix_MEL ))) %>% 
  dplyr::select(c(newclass, within_cluster_order))
split_instructions.rna <- as.character( rowidx$newclass )

heatmap <- Heatmap(tf_matrix_MEL, 
                   name = "Overlap Count", 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE,
                   show_row_names = FALSE, 
                   show_column_names = TRUE,
                   column_names_gp = grid::gpar(fontsize = 8),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_rot = 45,
                   row_gap = unit(1, "mm"),
                   border = TRUE,
                   row_title = "Genes", 
                   column_title = "Number of Overlap of TF and enhancers of DE genes, MEL enhancers",
                   column_dend_height = unit(1, "cm"),
                   row_split = split_instructions.rna)  # Use row_split to divide rows based on split_instructions

# Draw the heatmap
pdf(paste0(dir,"/heatmap_enhancer_TF_overlap_MEL_binary.pdf"),width=24,height=30)
draw(heatmap) # Adjust padding to ensure labels fit
dev.off()
# CREB1, RFX5, NFYA, NFYB seem to be enriched in cluster9.check out their expression levels.
make_paired_plot("NFYB",CPM_rna["NFYB",])
# which genes' enhancers overlap each TF?
gene_list <- unique( f[f$TF=="CREB1",]$gene )
gene_summary <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% gene_list) %>% 
  dplyr::select(c(newclass)) %>%
  group_by(newclass) %>%
  summarise(count_thisgene = n(), .groups = 'drop')
allgenes_summary <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% unique(f$gene)) %>% 
  dplyr::select(c(newclass)) %>%
  group_by(newclass) %>%
  summarise(count_allgenes = n(), .groups = 'drop')
df <- left_join(gene_summary, allgenes_summary)


#### Heatmap of TF enhancer Overlaps: #2 FRB enhancers ####
# Specify split instructions
aggregated_data_FRB <- f %>%
  group_by(gene, celltype, TF) %>%
  summarise(overlap_count = n(), .groups = 'drop') %>%
  dplyr::filter(celltype=="FRB") %>%
  dplyr::select(-celltype)

tf_matrix_FRB <- aggregated_data_FRB %>%
  pivot_wider(names_from = TF, values_from = overlap_count, values_fill = list(overlap_count = 0)) %>%
  dplyr::filter(gene %in% intersect(gene, rownames(km_res.rna.new))) %>%
  column_to_rownames("gene") %>%
  as.matrix()

rowidx <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% unique(rownames( tf_matrix_FRB ))) %>% 
  dplyr::select(c(newclass, within_cluster_order))
split_instructions.rna <- as.character( rowidx$newclass )

heatmap <- Heatmap(tf_matrix_FRB, 
                   name = "Overlap Count", 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE,
                   show_row_names = FALSE, 
                   show_column_names = TRUE,
                   column_names_gp = grid::gpar(fontsize = 8),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_rot = 45,
                   row_gap = unit(1, "mm"),
                   border = TRUE,
                   row_title = "Genes", 
                   column_title = "Number of Overlap of TF and enhancers of DE genes, FRB enhancers",
                   column_dend_height = unit(1, "cm"),
                   row_split = split_instructions.rna)  # Use row_split to divide rows based on split_instructions

# Draw the heatmap
pdf(paste0(dir,"/heatmap_enhancer_TF_overlap_FRB.pdf"),width=24,height=30)
draw(heatmap) # Adjust padding to ensure labels fit
dev.off()

#### Heatmap of TF enhancer Overlaps: #2 KRT enhancers ####
# Specify split instructions
aggregated_data_KRT <- f %>%
  group_by(gene, celltype, TF) %>%
  summarise(overlap_count = n(), .groups = 'drop') %>%
  dplyr::filter(celltype=="KRT") %>%
  dplyr::select(-celltype)

tf_matrix_KRT <- aggregated_data_KRT %>%
  pivot_wider(names_from = TF, values_from = overlap_count, values_fill = list(overlap_count = 0)) %>%
  dplyr::filter(gene %in% intersect(gene, rownames(km_res.rna.new))) %>%
  column_to_rownames("gene") %>%
  as.matrix()

rowidx <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% unique(rownames( tf_matrix_KRT ))) %>% 
  dplyr::select(c(newclass, within_cluster_order))
split_instructions.rna <- as.character( rowidx$newclass )

heatmap <- Heatmap(tf_matrix_KRT, 
                   name = "Overlap Count", 
                   cluster_rows = FALSE, 
                   cluster_columns = TRUE,
                   show_row_names = FALSE, 
                   show_column_names = TRUE,
                   column_names_gp = grid::gpar(fontsize = 8),
                   row_names_side = "left",
                   column_names_side = "top",
                   column_names_rot = 45,
                   row_gap = unit(1, "mm"),
                   border = TRUE,
                   row_title = "Genes", 
                   column_title = "Number of Overlap of TF and enhancers of DE genes, KRT enhancers",
                   column_dend_height = unit(1, "cm"),
                   row_split = split_instructions.rna)  # Use row_split to divide rows based on split_instructions

# Draw the heatmap
pdf(paste0(dir,"/heatmap_enhancer_TF_overlap_KRT.pdf"),width=24,height=30)
draw(heatmap) # Adjust padding to ensure labels fit
dev.off()

#### Differential Enrichment Analysis: commonly induced vs. MEL specific ####
# To identify TFs that are significantly associated with commonly-induced genes vs. cell-type-specific genes,
# Perform Fisher's Exact Test to compare the frequency of each TF's overlap with enhancers in commonly-induced genes vs. cell-type-specific genes
#### PCA ####
# PCA of cell-type-agnostic enhancers
tf_matrix_ctagnostic

#### Venn Diagram ####
# To visualize unique vs. shared TFs across commonly induced genes and cell-type-specific genes
gene_cluster_info <- data.frame(gene=rownames(km_res.rna.new), cluster=km_res.rna.new[,"newclass"])
f2 <- left_join(f, gene_cluster_info, by="gene") %>% na.omit()  

commonly_induced_genes <- f2 %>% dplyr::filter(cluster %in% c(1,2,3))
MEL_specific_genes <- f2 %>% dplyr::filter(cluster %in% c(8,9,10))

TF.union <- union(commonly_induced_genes$TF, MEL_specific_genes$TF)
TF.shared <- intersect(commonly_induced_genes$TF, MEL_specific_genes$TF)
TF.MEL_specific <- TF.union[which((!TF.union %in% TF.shared) & (TF.union %in% MEL_specific_genes$TF))]
TF.MEL_specific <- TF.union[which((!TF.union %in% TF.shared) & (TF.union %in% MEL_specific_genes$TF))]

tf_matrix_MEL <- aggregated_data_MEL %>%
  pivot_wider(names_from = TF, values_from = overlap_count, values_fill = list(overlap_count = 0)) %>%
  dplyr::filter(gene %in% intersect(gene, rownames(km_res.rna.new))) %>%
  column_to_rownames("gene") %>%
  as.matrix()

rowidx <- km_res.rna.new %>% 
  dplyr::filter(rownames(.) %in% unique(rownames( tf_matrix_MEL ))) %>% 
  dplyr::select(c(newclass, within_cluster_order))
split_instructions.rna <- as.character( rowidx$newclass )


