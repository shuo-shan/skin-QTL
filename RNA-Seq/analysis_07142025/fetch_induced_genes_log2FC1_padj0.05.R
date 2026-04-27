metagene_dict <- data.table::fread("/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt") %>% 
  dplyr::select(-id) %>% distinct()
metagene_dict_collapsed <- metagene_dict %>% group_by(meta) %>% 
  summarize(source=paste(name,collapse=','))

convert_meta_genes <- function(gene_vector) {
  # Identify meta genes in the input vector
  meta_genes <- gene_vector[grepl("meta", gene_vector)]
  non_meta_genes <- gene_vector[!grepl("meta", gene_vector)]
  
  # Convert meta genes to actual gene names using metagene_dict
  meta_to_gene <- metagene_dict$name[metagene_dict$meta %in% meta_genes]
  
  # Combine non-meta genes and converted meta genes, ensuring uniqueness
  converted_genes <- unique(c(non_meta_genes, meta_to_gene))
  
  return(converted_genes)
}

fetch_induced_genes <- function(ct, this_condition) {
  
  # # for debugging only
  # ct <- "FRB"
  # this_condition <- "IFNG"
  
  load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData")
  vars <- ls(pattern = "^deg")
  
  target_var <- vars[grepl(tolower(ct), vars, ignore.case = TRUE) & 
                       grepl(tolower(this_condition), vars, ignore.case = TRUE)]
  
  deg_obj <- get(target_var)
  
  gene_list <- deg_obj %>%
    dplyr::filter(log2FoldChange > 1) %>%
    pull(gene) %>%
    convert_meta_genes()
  
  return(gene_list)
}
