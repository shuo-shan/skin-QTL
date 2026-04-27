suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
  library(coloc)   # coloc 5.2.3
  library(ggplot2)
  library(patchwork)
  library(scales)
})

ct="FRB"
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)
this_gene="IRF3"
chunk_id_lookup <- data.table::fread(paste0(dir, "/data/gene_chunk_dict.txt"))
chunk_id <- unique(chunk_id_lookup[chunk_id_lookup$gene == this_gene, ]$chunk)
chunk_id <- sprintf("%03d", chunk_id)

pair_file <- paste0(dir, "/chunks/pairs_chunk_", chunk_id, ".tsv")
geno_file <- paste0(dir, "/chunks/genotype_pairs_chunk_", chunk_id, ".tsv")
vst_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
cpm_file  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"

modelstats_file <- paste0(dir, "/results/result_", chunk_id, ".tsv")
meta_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
rm(chunk_id_lookup)

# ------------------------------
# Load pairs / genotype / meta / CPM / VST
# ------------------------------
pairs <- fread(pair_file, header = TRUE) %>%
  dplyr::filter(gene_name == this_gene) %>%
  dplyr::mutate(
    SNP = stringr::str_to_lower(SNP_ID),
    key = paste0(gene_name, "_", SNP_ID)
  )

genotype_all <- fread(geno_file) %>%
  dplyr::mutate(ID = stringr::str_to_lower(ID)) %>%
  dplyr::filter(ID %in% pairs$SNP_ID)

meta_all <- readr::read_tsv(meta_file, show_col_types = FALSE)

VST_all <- fread(vst_file) %>%
  column_to_rownames("gene")

VST_subset <- VST_all %>%
  select(matches("FRB") & matches("PBS|IFNG"))

meta_subset <- meta_all %>%
  dplyr::filter(sample %in% colnames(VST_subset))

paired_donors <- table(meta_subset$donor) %>%
  as.data.frame() %>%
  dplyr::filter(Freq==2) %>%
  pull(Var1)

meta_subset_paired <- meta_subset %>%
  dplyr::filter(donor %in% paired_donors)

VST_subset_paired <- VST_subset %>%
  dplyr::select(meta_subset_paired$sample)

VST.PBS <- VST_subset_paired %>%
  dplyr::select(matches("PBS")) 
colnames(VST.PBS) <- gsub(".*_","", gsub("_FRB_PBS_3ctk_S1mod","", colnames(VST.PBS)))

VST.IFNG <- VST_subset_paired %>%
  dplyr::select(matches("IFNG"))
colnames(VST.IFNG) <- gsub(".*_","", gsub("_FRB_IFNG_3ctk_S1mod","", colnames(VST.PBS)))

colnames(VST.PBS)==colnames(VST.IFNG)

VST.delta <- VST.IFNG - VST.PBS
rownames(VST.delta) <- rownames(VST.PBS)

VST.delta.IRF3 <- as.data.frame(VST.delta)["IRF3",] %>% as.numeric()


# how are delta (IFNG - PBS) levels correlated between IRF3 and another gene
tag_gene <- "CGAS"
VST.delta.taggene <- as.data.frame(VST.delta)[tag_gene,] %>% as.numeric()
cor(VST.delta.IRF3, VST.delta.taggene)

# how are IFNG levels correlated
x <- VST.IFNG["IRF3",] %>% as.numeric()
y <- VST.IFNG[tag_gene,] %>% as.numeric()
cor(x, y, method="pearson")

# ------------------------------
# Fetch genome locus range
# ------------------------------
chr <- unique(pairs$gene_chr)
ciswindow_left <- unique(pairs$gene_start) - 500000
ciswindow_right <- unique(pairs$gene_start) + 500000
this_range <- paste(chr, ciswindow_left, ciswindow_right, sep="_")

