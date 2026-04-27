#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(magrittr)
  library(readr)
  library(patchwork)
  library(tibble)
})


# ---------------------- Set-up --------------------- ####
ct <- "MEL"
basedir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
CPM_FILE="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"

# get eGenes and reGenes for each condition.
this_dir <- paste0(basedir, "/eigenMT/results")

QTL_genes <- list()

QTL_genes$PBS_eQTL  <- fread(paste0(this_dir, "/MEL_PBS_eQTL_gene_fdr05_genelist.txt"), header = FALSE)$V1
QTL_genes$IFNG_eQTL <- fread(paste0(this_dir, "/MEL_IFNG_eQTL_gene_fdr05_genelist.txt"), header = FALSE)$V1
QTL_genes$IFNB_eQTL<- fread(paste0(this_dir, "/MEL_IFNB_eQTL_gene_fdr05_genelist.txt"), header = FALSE)$V1
QTL_genes$TNF_eQTL  <- fread(paste0(this_dir, "/MEL_TNF_eQTL_gene_fdr05_genelist.txt"), header = FALSE)$V1

QTL_genes$IFNG_reQTL<- fread(paste0(this_dir, "/MEL_IFNG_reQTL_gene_fdr05_genelist.txt"), header = FALSE)$V1
QTL_genes$IFNB_reQTL<-fread(paste0(this_dir, "/MEL_IFNB_reQTL_gene_fdr05_genelist.txt"), header = FALSE)$V1
QTL_genes$TNF_reQTL <- fread(paste0(this_dir, "/MEL_TNF_reQTL_gene_fdr05_genelist.txt"), header = FALSE)$V1

# find all expressed genes


# find the TFs in each gene list ####
TF_list <- fread("/pi/manuel.garber-umw/human/skin/eQTLs/literature/Lambert_2018_human_TFs.txt", header = F)$V1
TF_list <- c(TF_list, "ZNF330")
QTL_genes_TFs <- lapply(QTL_genes, intersect, TF_list)
QTL_genes_TFs_sorted <- lapply(QTL_genes_TFs, function(x) sort(unique(x)))
QTL_gene_counts <- sapply(QTL_genes, length)
tf_counts <- sapply(QTL_genes_TFs, length)
tf_frac <- mapply(function(allg, tfg) length(tfg)/length(allg),
                  QTL_genes, QTL_genes_TFs)
QTL_gene_counts
tf_counts
tf_frac

# write TF to file
out <- file.path(basedir, "/data/","QTL_genes_TFs_sorted.txt")

con <- file(out, open = "wt")
for (nm in names(QTL_genes_TFs_sorted)) {
  writeLines(paste0("## ", nm), con)
  writeLines(QTL_genes_TFs_sorted[[nm]], con)
  writeLines("", con)
}
close(con)



# who are the TFs? ####
TF_table <- data.table(
  set = names(QTL_genes_TFs),
  n_TF = sapply(QTL_genes_TFs, length),
  TFs = sapply(QTL_genes_TFs, function(x) paste(sort(unique(x)), collapse = ", "))
)
TF_long <- rbindlist(lapply(names(QTL_genes_TFs), function(nm) {
  data.table(set = nm, TF = QTL_genes_TFs[[nm]])
}))
TF_long

# TFs shared across many sets (likely core regulators)
tf_freq <- sort(table(TF_long$TF), decreasing = TRUE)
head(tf_freq, 30)

# TFs unique to a specific set (likely cytokine-specific program)
tf_by_set <- split(TF_long$TF, TF_long$set)
tf_unique <- lapply(names(tf_by_set), function(s) {
  setdiff(tf_by_set[[s]], unique(unlist(tf_by_set[names(tf_by_set) != s])))
})
names(tf_unique) <- names(tf_by_set)
lapply(tf_unique, head, 30)

# cytokines already up in PBS and do they future go up in cytokines?



