#!/usr/bin/env Rscript
# find all eGene and lead QTL pairs (primary and secondary signals), and permute for transQTL mapping

Sys.setenv(TZ = "America/New_York")
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(future.apply)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)

ct <- args[1] # MEL, KRT, FRB
this_condition <- args[2]
this_QTLtype <- args[3]

# # for debugging only
# ct <- "FRB"
# this_condition <- "IFNG"
# this_QTLtype <- "eQTL"

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
output_file_trans_pairs <- paste0(dir,"/transQTL/eGene_QTL_pairs/",ct,"_",this_condition,"_",this_QTLtype,"_trans_pairs.txt")
output_file_QTLtags <- paste0(dir,"/transQTL/QTL_tags/",ct,"_",this_condition,"_",this_QTLtype,"_QTLtags.txt")

CPM_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc

GENE_CHUNK_DICT_FILE <- paste0(dir,"/data/gene_chunk_dict.txt")
MAIN_QTL_RES_FILE <- paste0(dir,"/eigenMT/results/",ct,"_",this_condition,"_",this_QTLtype,"_gene_fdr05_table.txt")
SECONDARY_QTL_RES_FILE <- paste0(dir,"/conditional_analysis_round1/eigenMT/results/",ct,"_",this_condition,"_",this_QTLtype,"_gene_fdr05_table.txt")

# ---- fetch a table of eGenes and their lead SNPs (primary + secondary) ----
# ---- select QTLs with high effect size ----
chunk_dict <- fread(GENE_CHUNK_DICT_FILE, header=T)
main_QTL_res <- fread(MAIN_QTL_RES_FILE, header=T) %>%
  dplyr::filter(q_gene < 0.005) %>%
  left_join(chunk_dict, by="gene") %>%
  dplyr::mutate(SNPtag=paste0("primary_signal_of_",gene,"_as_",ct,"_",this_condition,"_",this_QTLtype))

secondary_QTL_res <- fread(SECONDARY_QTL_RES_FILE, header=T) %>%
  dplyr::filter(q_gene < 0.005) %>%
  left_join(chunk_dict, by="gene") %>%
  dplyr::mutate(SNPtag=paste0("secondary_signal_of_",gene,"_as_",ct,"_",this_condition,"_",this_QTLtype))

QTL_res <- rbind(main_QTL_res, secondary_QTL_res) %>%
  arrange(gene, SNPtag) %>%
  dplyr::select(c(gene, lead_snp, chunk))

SNP_tag_dict <- rbind(main_QTL_res, secondary_QTL_res) %>%
  arrange(gene, SNPtag) %>%
  dplyr::select(c(lead_snp, SNPtag)) %>%
  group_by(lead_snp) %>% 
  summarize(SNPtag=paste(SNPtag,collapse=',')) %>%
  set_colnames(c("snp","SNPtag"))

# ---- fetch genes that are expressed genes and DE genes (if cytokine) ----
trans_acting_genes_table <- fread("/pi/manuel.garber-umw/human/skin/eQTLs/literature/trans_acting_genes/compiled_trans_acting_candidate_genes_and_category.txt", header=T)

path.Rscript="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_expressed_genes.R"
source(path.Rscript)
expressed_genes <- fetch_expressed_genes(ct, this_condition)
kept_genes <- expressed_genes

if (this_condition != "PBS") {
  path.Rscript="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_DE_genes_log2FC1_padj0.05.R"
  source(path.Rscript)
  DE_genes <- unique(fetch_DE_genes(ct, this_condition))
  kept_genes <- intersect(expressed_genes, DE_genes)
  # 
  # path.Rscript="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_induced_genes_log2FC1_padj0.05.R"
  # source(path.Rscript)
  # induced_genes <- fetch_induced_genes(ct, this_condition)
  # kept_genes <- intersect(expressed_genes, induced_genes)
}

#expressed_regulators <- intersect(trans_acting_genes_table$gene, kept_genes)
#gene_category_table <- trans_acting_genes_table %>%
#  dplyr::filter(gene %in% expressed_regulators) %>%
#  set_colnames(c("gene","gene_category"))

cat("Total expressed/DE genes:", length(kept_genes), "\n")
#cat("Total expressed/DE trans-regulator genes:", length(expressed_regulators), "\n")

# ---- compile genome position annotation for genes and QTLs ----
# genes
gene_list <- kept_genes
gene_chunk_file <- fread(paste0(dir,"/data/gene_chunk_dict.txt")) %>%
  dplyr::filter(gene %in% gene_list)
chunk_list <- sort(unique(gene_chunk_file$chunk))
genes_info_all <- data.frame()
for (i in 1:length(chunk_list)) {
  this_chunk <- chunk_list[i]
  chunk_id <- sprintf("%03d", this_chunk)
  message(paste0("processing ", chunk_id))
  
  pairs_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
  pairs <- fread(pairs_file, header = TRUE)
  colnames(pairs)[4] <- "gene"
  colnames(pairs)[11] <- "snp"
  
  genes_info <- pairs %>% 
    dplyr::select(c(gene_chr, gene_start, gene_end, gene)) %>%
    dplyr::filter(gene %in% gene_list) %>%
    dplyr::mutate(gene_chunk=chunk_id) %>%
    dplyr::distinct() 
  
  genes_info_all <- rbind(genes_info_all, genes_info)
  rm(pairs)
}

# QTLs
QTL_list <- unique(QTL_res$lead_snp)
chunk_list <- sort(unique(QTL_res$chunk))
QTL_info_all <- data.frame()
for (i in 1:length(chunk_list)){

  this_chunk <- chunk_list[i]
  chunk_id <- sprintf("%03d", this_chunk)
  message(paste0("processing ", chunk_id))
  
  pairs_file <- paste0(dir,"/chunks/pairs_chunk_",chunk_id,".tsv")
  pairs <- fread(pairs_file, header = TRUE)
  colnames(pairs)[4] <- "gene"
  colnames(pairs)[11] <- "snp"
  
  QTL_info <- pairs %>%
    dplyr::select(c(SNP_chr, SNP_start, SNP_end, snp)) %>%
    dplyr::filter(snp %in% QTL_list) %>%
    dplyr::mutate(SNP_chunk=chunk_id) %>%
    dplyr::distinct() %>%
    left_join( . , SNP_tag_dict, by="snp")
  
  QTL_info_all <- rbind(QTL_info_all, QTL_info) 
  
  rm(pairs)
}
# keep the first row for each SNP and drop the rest. just to grab the chunk id for SNP for loading gentoype table down the line.
QTL_info_all <- QTL_info_all %>% dplyr::distinct(snp, .keep_all = TRUE)

# ---- compile transQTL test pairs by iterating through the gene table and QTL table ----
gene_QTL_trans_pairs_table <- tidyr::crossing(
  genes_info_all,
  QTL_info_all
) %>%
  mutate(
    cis_trans_category = case_when(
      gene_chr != SNP_chr ~ "different_chr",
      gene_chr == SNP_chr & abs(gene_start - SNP_start) >= 500000 ~ "same_chr_outside_500kb_window",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cis_trans_category))

# ---- chunk pairs and write to file ----
threshold <- 100000

# count rows per gene
gene_summary <- gene_QTL_trans_pairs_table %>%
  count(gene, name = "n_rows") %>%
  arrange(gene)

# assign chunks
assign_chunks <- function(gene_summary, threshold) {
  
  if (sum(gene_summary$n_rows) <= threshold) {
    gene_summary %>%
      mutate(chunk_id = sprintf("%03d", 1))
    
  } else {
    
    chunk_list <- list()
    df_remain <- gene_summary
    i <- 1
    
    while (nrow(df_remain) > 0) {
      
      cs <- cumsum(df_remain$n_rows)
      j_exceed <- which(cs > threshold)[1]
      
      if (is.na(j_exceed)) {
        idx <- seq_len(nrow(df_remain))
      } else {
        idx <- seq_len(j_exceed - 1)
        if (length(idx) == 0) idx <- 1
      }
      
      chunk_genes <- df_remain$gene[idx]
      chunk_name <- sprintf("%03d", i)
      
      chunk_list[[i]] <- tibble(
        gene = chunk_genes,
        chunk_id = chunk_name
      )
      
      df_remain <- df_remain[-idx, , drop = FALSE]
      i <- i + 1
    }
    
    bind_rows(chunk_list)
  }
}

gene_chunk_dict <- assign_chunks(gene_summary, threshold)

# join back to full table
gene_QTL_trans_pairs_table_chunked <- gene_QTL_trans_pairs_table %>%
  left_join(gene_chunk_dict, by = "gene") %>%
  relocate(chunk_id, .before = everything())

# split into a list of data frames by chunk_id
chunk_list <- gene_QTL_trans_pairs_table_chunked %>%
  group_split(chunk_id)

chunk_names <- gene_QTL_trans_pairs_table_chunked %>%
  distinct(chunk_id) %>%
  arrange(chunk_id) %>%
  pull(chunk_id)

# write each chunk table
for (i in seq_along(chunk_list)) {
  this_chunk <- chunk_list[[i]]
  this_name <- chunk_names[i]
  output_file_trans_pairs_chunk <- paste0(dir, "/transQTL/chunks/",this_condition,"/", this_QTLtype,"/gene_QTL_pairs_chunk_", this_name, ".tsv")
  
  write_tsv(this_chunk,output_file_trans_pairs_chunk)
  message(paste0("Writen output to ", output_file_trans_pairs_chunk))
}

# write gene chunk dict
output_file_chunk_dict <- paste0(dir, "/transQTL/data/",ct,"_",this_condition,"_", this_QTLtype,"_gene_transQTL_chunk_dict.txt")
fwrite(gene_chunk_dict, output_file_chunk_dict, quote=F, sep="\t")
message(paste0("Writen output to ", output_file_chunk_dict))

# write full trans pairs table
fwrite(gene_QTL_trans_pairs_table, file=output_file_trans_pairs, sep="\t", quote=F)
message(paste0("Written output to ", output_file_trans_pairs))

# write SNP annotations to which eGenes it associates with
fwrite(SNP_tag_dict, file = output_file_QTLtags, sep="\t", quote=F)
message(paste0("Written output to ", output_file_QTLtags))
