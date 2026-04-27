library(tidyverse)
library(DESeq2)
library(magrittr)
library(stringr)
library(rtracklayer)
library(GenomicRanges)


# in cluster:
# cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/gencode
# curl -sS "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz" | zcat > "gencode.v34.primary_assembly.annotation.gtf"


DEGs <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/2025models/data/all_DEGs.txt", header=F)$V1
#gencode_genes <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/literature/gencode/gencode.v34.primary_assembly.annotation.gtf", header=F)
gtf.gencodev34 <- import("~/Downloads/nl/human/skin/eQTLs/literature/gencode/gencode.v34.primary_assembly.annotation.gtf")
gtf.gencodev34.genes <- unique(gtf.gencodev34[gtf.gencodev34$type == "gene"]$gene_name)

vars_deg <- ls(pattern = "^deg")
gene_lists <- lapply(vars_deg, function(var) {
  df <- get(var)
  # keep only significant, upregulated genes
  df %>%
    filter(log2FoldChange > 1, padj < 0.05) %>%
    pull(gene)   # or `rownames(df)` if genes are rownames
})
DEGs <- unique(convert_meta_genes(unlist(gene_lists)))

write(gene_list, "~/Downloads/nl/human/skin/eQTLs/DREG/2025models/data/all_DEGs.txt")
missing <- setdiff(gene_list, gtf.gencodev34.genes)

genes <- gtf.gencodev34[gtf.gencodev34$type == "gene"]

# TSS position (1-based genomic coordinate)
strd    <- as.character(strand(genes))
tss_pos <- ifelse(strd == "+", start(genes), end(genes))

# build BED6 (+ an extra gene_id column)
bed <- data.frame(
  chrom  = as.character(seqnames(genes)),
  start  = pmax(as.integer(tss_pos) - 1L, 0L),  # 0-based
  end    = as.integer(tss_pos),                 # half-open
  name   = if (!is.null(genes$gene_name)) as.character(genes$gene_name) else as.character(genes$gene_id),
  score  = 0,
  strand = strd,
  gene_id = if (!is.null(genes$gene_id)) as.character(genes$gene_id) else NA_character_
)

# sort for readability
bed <- bed[order(bed$chrom, bed$start), ]
bed <- bed[!duplicated(bed$name), , drop = FALSE] # One TSS per gene_name (keep first)

# write as plain BED (no header)
write.table(bed, file = "~/Downloads/nl/human/skin/eQTLs/literature/gencode/gencode_v34_gene_TSS.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# write BED file of DEGs (padj 0.01, absLog2FC 1)
bed.subset <- bed %>% dplyr::filter(name %in% DEGs)
length(unique(bed.subset$name))
length(DEGs)

write.table(bed.subset, file = "~/Downloads/nl/human/skin/eQTLs/literature/gencode/gencode_v34_allDEGs_padj0.01_absLog2FC1_TSS.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# write BED file of expressed genes
bed.subset <- bed %>% dplyr::filter(name %in% unique(convert_meta_genes(rownames(counts.expressed))))
length(unique(bed.subset$name))
length(unique(convert_meta_genes(rownames(counts.expressed))))

write.table(bed.subset, file = "~/Downloads/nl/human/skin/eQTLs/literature/gencode/gencode_v34_allExpressedGenes_TSS.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
