# How many transcription factors have a QTL?
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(ggtext)
  library(scales)
})



# =======1. PATHS ===================================
dir="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models"
tf_file="/pi/manuel.garber-umw/human/skin/eQTLs/literature/Lambert_2018_human_TFs.txt"
deg_file    <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData"
metagene_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt"
out_dir     <- paste0(dir,"/data")

tf_list <- unique(fread(tf_file, header=F)$V1)
load(deg_file)
metagene_dict <- data.table::fread(metagene_file) %>%
  dplyr::select(-id) %>%
  dplyr::distinct()

# =======2. UTILITY =============================
metagene_dict <- data.table::fread(metagene_file) %>%
  dplyr::select(-id) %>%
  dplyr::distinct()

convert_meta_genes <- function(gene_vector) {
  this_metagenes_idx <- which(grepl("meta", gene_vector))
  temp1 <- data.frame(
    meta = gene_vector[-this_metagenes_idx],
    name = gene_vector[-this_metagenes_idx]
  )
  temp2 <- data.frame(meta = gene_vector[this_metagenes_idx]) %>%
    dplyr::left_join(metagene_dict, by = "meta")
  rbind(temp1, temp2) %>%
    magrittr::set_colnames(c("meta", "gene"))
}

# Helper to get DE genes for a specific ct + stim
get_deg <- function(ct_lower, stim_lower, direction = c("both", "up", "down")) {
  direction <- match.arg(direction)
  var <- ls(pattern = paste0("^deg\\.", ct_lower, "\\.", stim_lower),
            envir = .GlobalEnv)
  if (length(var) == 0) {
    cat(sprintf("  No DE object found for deg.%s.%s\n", ct_lower, stim_lower))
    return(NULL)
  }
  deg_raw <- get(var)
  colnames(deg_raw)[1] <- "meta"
  
  # Convert meta-gene IDs to gene names
  meta_dict <- convert_meta_genes(deg_raw$meta)
  deg_converted <- dplyr::left_join(meta_dict, deg_raw, by = "meta")
  
  # Filter by direction
  if (direction == "up")   deg_converted <- dplyr::filter(deg_converted, log2FoldChange > 0)
  if (direction == "down")  deg_converted <- dplyr::filter(deg_converted, log2FoldChange < 0)
  
  dplyr::select(deg_converted, gene, log2FoldChange, padj) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
}

# ======== 3. Get QTL result ===========================
get_QTL_result <- function(celltype, stim, QTLtype) {
  # total TFs tested 
  total_gene_lst <- unique(fread(paste0(dir,"/",celltype,"/eigenMT/results/",celltype,"_",stim,"_",QTLtype,".eigenMT.txt"))$gene)
  total_tf_tested <- intersect(total_gene_lst, tf_list)
  
  # TFs with QTL
  this_res_file <- paste0(dir,"/",celltype,"/eigenMT/results/",celltype,"_",stim,"_",QTLtype,"_gene_fdr05_table.txt")
  this_qtl_res <- fread(this_res_file, header=T)
  this_qtl_res_tf <- this_qtl_res %>% dplyr::filter(gene %in% tf_list)
  
  # TFs with QTL that are DE genes
  ct_lower <- tolower(celltype)
  cytokines <- c("IFNB", "IFNG", "TNF")
  
  deg_flags <- purrr::map_dfr(tolower(cytokines), function(s) {
    deg <- get_deg(ct_lower, s)
    if (is.null(deg)) return(NULL)
    tibble::tibble(
      gene        = deg$gene,
      cytokine    = toupper(s),
      is_DE       = TRUE,
      DE_direction = dplyr::case_when(
        deg$log2FoldChange > 0 ~ "up",
        deg$log2FoldChange < 0 ~ "down",
        TRUE                   ~ "unchanged"
      ),
      log2FC      = deg$log2FoldChange,
      DE_padj     = deg$padj
    )
  }) %>%
    dplyr::distinct(gene, cytokine, .keep_all = TRUE) %>%
    dplyr::filter(cytokine==stim) %>%
    dplyr::select(gene, is_DE, DE_direction)
  
  if (nrow(this_qtl_res_tf)>0){
    this_qtl_res_tf2 <- left_join(this_qtl_res_tf, deg_flags, by="gene") %>% dplyr::filter(!is.na(is_DE))
    n_tf_qtl_DE <- nrow(this_qtl_res_tf2)
    n_tf_qtl_DE_up <- nrow(this_qtl_res_tf2[which(this_qtl_res_tf2$DE_direction=="up"),])
    n_tf_qtl_DE_down <- nrow(this_qtl_res_tf2[which(this_qtl_res_tf2$DE_direction=="down"),])
    
  } else {
    n_tf_qtl_DE = 0
    n_tf_qtl_DE_up = 0
    n_tf_qtl_DE_down = 0
  }
  
  # summarize count
  n_total_genes_tested <- length(total_gene_lst)
  n_total_TFs_tested <- length(total_tf_tested)
  n_this_qtl_res_tf <- length(unique(this_qtl_res_tf$gene))

  return(
    data.frame(
      celltype             = celltype,
      stim                 = stim,
      QTLtype              = QTLtype,
      n_total_genes_tested = n_total_genes_tested,
      n_total_TFs_tested   = n_total_TFs_tested,
      n_this_qtl_res_tf    = n_this_qtl_res_tf,
      n_tf_qtl_DE          = n_tf_qtl_DE,
      n_tf_qtl_DE_up       = n_tf_qtl_DE_up,
      n_tf_qtl_DE_down     = n_tf_qtl_DE_down
    )
  )
}

all_res <- data.frame()

for (celltype in c("FRB", "MEL", "KRT")) {
  for (stim in c("PBS", "IFNG", "IFNB", "TNF")) {
    QTLtype <- "eQTL"
    this_res <- get_QTL_result(celltype, stim, QTLtype)
    all_res <- rbind(all_res, this_res)
  }
}
for (celltype in c("FRB", "MEL", "KRT")) {
  for (stim in c("IFNG", "IFNB", "TNF")) {
    QTLtype <- "reQTL"
    print(celltype)
    print(stim)
    this_res <- get_QTL_result(celltype, stim, QTLtype)
    all_res <- rbind(all_res, this_res)
  }
}

# summarize ----
all_res$ratio_tf_with_qtl_over_tf_tested <- round(all_res$n_this_qtl_res_tf/all_res$n_total_TFs_tested, 3)

# write ----
fwrite(all_res, file=paste0(out_dir,"/fig2_TF_with_QTL.txt"), sep="\t", quote=F)
