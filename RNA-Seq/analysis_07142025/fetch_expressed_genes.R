fetch_expressed_genes <- function(ct, this_condition) {
  
  # # for debugging only
  # ct <- "FRB"
  # this_condition <- "IFNG"
  
  dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
  CPM_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
  META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
  
  # load
  CPM <- fread(CPM_FILE) %>%
    dplyr::filter(final_gene!="ARMCX5-GPRASP2") %>%
    column_to_rownames("final_gene") %>%
    dplyr::select(-c(gene,name))
  meta <- fread(META_FILE)
  
  # subset
  this_samples <- meta %>%
    dplyr::filter(celltype==ct & condition==this_condition) %>%
    pull(sample)
  
  CPM_subset <- CPM %>%
    dplyr::select(all_of(this_samples))
  
  # filter by expression
  n_min_donor <- ceiling(0.20 * ncol(CPM_subset))
  keep_genes <- (rowSums(CPM_subset >= 1) >= n_min_donor) & 
    (rowSums(CPM_subset >= 10) >= 5)
  CPM_keep <- CPM_subset[keep_genes,]
  
  return(rownames(CPM_keep))
}
