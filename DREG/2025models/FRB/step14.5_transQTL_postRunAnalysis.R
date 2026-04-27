# performs analysis after transQTL mapping and multiple testing correction
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
  library(igraph)
  library(ggraph)
  library(tidyverse)
  library(gridExtra)
  library(grid)
  library(GenomicRanges)
})

# ------ Load data
ct <- "FRB"
condition <- "TNF"
QTLtype <- "eQTL"

dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
CPM_FILE      <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
VST_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt"
META_FILE     <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"   # columns: sample, donor, condition, etc
PEER_FILE     <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/peer/peer_factors/peer_factors_",ct,"_PBS-IFNG-IFNB-TNF.tsv")             # columns: sample, PEER1, PEER2, ...
GENO_PCS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt"      # columns: donor, PC1, PC2, ...
trans_acting_genes_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/literature/trans_acting_genes/compiled_trans_acting_candidate_genes_and_category.txt"
TRANSQTL_RES_FILE <- paste0(dir,"/transQTL/resultsBHcorrected/",ct,"_",condition,"_",QTLtype,"_gene_fdr05_table.txt")
tf_file="/pi/manuel.garber-umw/human/skin/eQTLs/literature/Lambert_2018_human_TFs.txt"
gene_bed_file = paste0(dir,"/data/all_genes.bed")
GENOTYPE_FILE <- paste0(dir,"/transQTL/QTL_tags/",ct,"_",condition,"_",QTLtype,"_SNPs_genotype.txt")

MAIN_CISEQTL_RES_FILE <- paste0(dir,"/eigenMT/results/",ct,"_",condition,"_",QTLtype,"_gene_fdr05_table.txt")
SECONDARY_CISEQTL_RES_FILE <- paste0(dir,"/conditional_analysis_round1/eigenMT/results/",ct,"_",condition,"_",QTLtype,"_gene_fdr05_table.txt")

# ------ Load significant transQTL table ------
res <- fread(TRANSQTL_RES_FILE) %>% 
  tidyr::separate_rows(SNPtag,sep=",")
res$cisGene <-  sub("_.*", "", sub(".*of_", "", res$SNPtag))
idx <- which(colnames(res)=="gene")
colnames(res)[idx] <- "transGene"

# add MAF info ----
genotype_bed <- fread(GENOTYPE_FILE)
gene_bed <- fread(gene_bed_file)
sample_cols <- names(genotype_bed)[6:ncol(genotype_bed)]
genotype_bed <- genotype_bed %>%
  mutate(MAF = {
    alt_count <- rowSums(pick(all_of(sample_cols)), na.rm = TRUE)
    n_alleles  <- 2 * rowSums(!is.na(pick(all_of(sample_cols))))
    alt_freq   <- alt_count / n_alleles
    pmin(alt_freq, 1 - alt_freq)
  })
res <- left_join( res, genotype_bed[,c("ID","MAF")], by=c("snp"="ID"))

# add lead cis-eQTL info ----
ciseQTL_res <- fread(MAIN_CISEQTL_RES_FILE) %>%
  dplyr::select(c(gene, q_gene)) %>%
  set_colnames(c("gene","q_gene_lead_ciseQTL"))
res <- left_join(res, ciseQTL_res, by=c("cisGene"="gene"))

# add SNP number of transGene info ----
df1 <- res[, c("snp","transGene")] %>%
  distinct()
df2 <- as.data.frame(table(df1$snp))
colnames(df2) <- c("snp","n_transGene")
res <- left_join(res, df2)
rm(df1, df2)

# cis-Gene annotation ----
trans_acting_genes_table <- fread(trans_acting_genes_file) %>%
  set_colnames(c("gene","cisGene_category"))
res <- left_join(res, trans_acting_genes_table, by=c("cisGene"="gene"))

# trans-Gene annotation ----
trans_acting_genes_table <- fread(trans_acting_genes_file) %>%
  set_colnames(c("gene","transGene_category"))
res <- left_join(res, trans_acting_genes_table, by=c("transGene"="gene"))

# add GWAS coloc information ----
dir_coloc_genes <- paste0(dir,"/coloc/summary")
coloc_files <- list.files(dir_coloc_genes, 
                          pattern = paste0("^coloc_.*_", condition, "_", QTLtype, "\\.txt$"),
                          full.names = TRUE)
coloc_summary <- coloc_files %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE)) %>%
  dplyr::filter(GWAS_trait != "height") %>%
  dplyr::filter(PP.H4 > 0.7)

coloc_summary_filtered <- coloc_summary %>%
  right_join(res, . , by=c("cisGene"="gene")) %>%
  dplyr::filter(!is.na(transGene))

df <- coloc_summary_filtered %>%
  dplyr::select(c(snp, cisGene, GWAS_trait, PP.H4)) %>%
  distinct()
res <- left_join(res, df, by=c("snp","cisGene"))

# add ATACseq peak overlap information ----
# convert SNP table to GRanges
my_snp_table <- genotype_bed %>%
  dplyr::filter(ID %in% unique(res$snp)) %>%
  dplyr::mutate(start=POS-1) %>%
  dplyr::select(c(CHROM,start,POS,ID)) %>%
  set_colnames(c("chr","start","end","snp"))

snp_gr <- my_snp_table %>%
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

# turn ATACseq peaks into GRanges
atac_table_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline/masterPeaks/ATACseq_peak_annotation.tsv"
atac_table <- fread(atac_table_file)

atac_gr <- atac_table %>%
  makeGRangesFromDataFrame(
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

# findOverlaps = bedtools intersect
hits <- findOverlaps(snp_gr, atac_gr)

# merge result
atac_result <- bind_cols(
  my_snp_table[queryHits(hits), ],
  atac_table[subjectHits(hits), ]
) %>%
  dplyr::select(!starts_with(c("chr", "start", "end")))

res <- left_join( res, atac_result, by="snp")

# examine results
df <- res %>%
  dplyr::select(c(snp,peak_name,open_FRB_PBS,open_FRB_IFNG,peakDynamic_FRB_IFNG)) %>%
  distinct() %>%
  na.omit(peak_name)
table(df$peakDynamic_FRB_IFNG)



# add TFBS information  ----




# ------- POSTER CHILD CANDIDATES #1 ----------
poster_child_candidates <- res %>%
  dplyr::filter(!is.na(GWAS_trait)) %>%          # has coloc hit
  dplyr::filter(PP.H4 > 0.7) %>%              # strong coloc
  dplyr::filter(!is.na(cisGene_category) | !is.na(transGene_category)) %>%  # cis gene is TF
  dplyr::filter(abs(beta) > 0.1)    %>%       # strong effect size - set your threshold
  dplyr::filter(p < 1e-5) %>%                    # strong trans-QTL stats
  dplyr::filter(q_gene_lead_ciseQTL < 1e-5) %>%  # strong cis-QTL stats  
  dplyr::arrange(p, desc(PP.H4))                 # rank by best stats

# is the poster child transGene also a cisGene?
dig_deeper <- res %>%
  dplyr::filter(cisGene %in% unique(poster_child_candidates$transGene))

# ------- POSTER CHILD CANDIDATES #2 ----------
poster_child_candidates <- res %>%
  dplyr::filter(open_FRB_PBS==TRUE) %>%
  dplyr::filter(!is.na(cisGene_category) | !is.na(transGene_category)) %>%  # cis gene is TF
  dplyr::filter(abs(beta) > 0.1)    %>%       # strong effect size - set your threshold
  dplyr::filter(p < 1e-5) %>%                    # strong trans-QTL stats
  dplyr::filter(q_gene_lead_ciseQTL < 1e-5) %>%  # strong cis-QTL stats  
  dplyr::arrange(p, desc(PP.H4))                 # rank by best stats
# "rs74000566" "rs10877020" "rs12154043" "rs76730390" "rs1010167"  "rs12790010"



# ------ 1. all transQTLs: overlap with trans-acting genes ------
trans_acting_genes_table <- fread(trans_acting_genes_file)
# any trans-gene that's known to be trans-acting?
df <- left_join(res, trans_acting_genes_table, by=c("transGene"="gene")) %>% na.omit()
View(df[,c("transGene","categories")] %>% distinct())
message(paste0("there's ",length(unique(df$transGene))," trans-acting gene(s) in trans-gene of this file"))

# any cis-gene that's known to be trans-acting?
df <- left_join(res, trans_acting_genes_table, by=c("cisGene"="gene")) %>% na.omit()
View(df[,c("cisGene","categories")] %>% distinct())
message(paste0("there's ",length(unique(df$cisGene))," trans-acting gene(s) in cis-gene of this file"))

# ------ 2. all transQTLs: overlap with TFs ---------
tf_list <- unique(fread(tf_file, header=F)$V1)

transQTLSNP_eGene_list <- unique(sub("_.*", "", sub(".*of_", "", unique(df$SNPtag))))

df$text <- paste0(df$cisGene,"-",df$lead_snp,"-",df$gene)
df$text

length(tf_list)
length(transQTLGene_list)
length(intersect(tf_list, transQTLGene_list))
length(transQTLSNP_eGene_list)
length(intersect(tf_list, transQTLSNP_eGene_list))

res_tf <- df %>%
  dplyr::filter(cisGene %in% tf_list)

# around 5% of genes with a trans-eQTL are TFs.
# ------ 3. all transQTLs: genes overlap with coloc genes --------
dir_coloc_genes <- paste0(dir,"/coloc/summary")
coloc_files <- list.files(dir_coloc_genes, 
                          pattern = paste0("^coloc_.*_", condition, "_", QTLtype, "\\.txt$"),
                          full.names = TRUE)
coloc_summary <- coloc_files %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE)) %>%
  dplyr::filter(GWAS_trait != "height") %>%
  dplyr::filter(PP.H4 > 0.7)

coloc_summary_filtered <- coloc_summary %>%
  right_join(res, . , by=c("cisGene"="gene")) %>%
  dplyr::filter(!is.na(transGene))

df <- coloc_summary_filtered %>%
  dplyr::select(c(snp, cisGene, GWAS_trait, n_transGene, PP.H4, cisGene_category)) %>%
  distinct()

fwrite(df, paste0(dir,"/transQTL/coloc_table_",ct,"_",condition,"_",QTLtype,".txt"), sep="\t", quote=F, col.names =F)
nrow(df)
length(unique(coloc_summary_filtered$cisGene))

# ------ 4. SNP level analysis #0: look at lead SNPs ----
lead_SNPs <- unique(res$snp)
res %>%
  count(snp) %>%
  ggplot(aes(n)) +
  geom_histogram() +
  labs(x="n trans genes associated with a SNP", y="n SNP") +
  ggtitle(paste0(ct," ",condition," ",QTLtype," transQTLs"))

snp_count <- table(res$snp) %>%
  as.data.frame() %>%
  arrange(desc(Freq))

# ------------------------------------------------------------------------------------------
# ------ SNP level analysis: plot SNP gene pairs and overlap with ATACseq peaks ----
# ------------------------------------------------------------------------------------------
# after step2, in terminal, type this:
# cd ${DIR}
# bash step14.6_transQTL_compile_result_for_snp.sh snp ct
snp_list <- res %>%
  arrange(desc(n_transGene)) %>%
  pull(snp) %>%
  unique 

write(snp_list, paste0(dir,"/transQTL/snps.txt"))

snp_list <- unique(df$snp)
# 3 to 23
for (i in 1:16){
  this_snp=snp_list[i]
  
  message(paste0("processing ",this_snp))
dir.create(paste0(dir, "/transQTL/temp_output/", this_snp), recursive = TRUE, showWarnings = FALSE)

# step1. compiled result table ----
slice <- res %>%
  dplyr::filter(snp==this_snp) %>%
  dplyr::arrange(p)
fwrite(slice, file=paste0(dir,"/transQTL/temp_output/",this_snp,"/compiled_table_",ct,"_",condition,"_",QTLtype,".txt"), quote=F, sep="\t")

slice_pretty <- slice %>% mutate(across(where(is.numeric), ~signif(., 3)))

pdf(paste0(dir,"/transQTL/temp_output/",this_snp,"/compiled_table_",ct,"_",condition,"_",QTLtype,".pdf"),
    width=20, height=nrow(slice)*0.35 + 1.5)
grid.newpage()
grid.table(slice_pretty, rows=NULL,
           theme=ttheme_default(
             base_size=6,
             colhead=list(fg_params=list(fontsize=7, fontface="bold")),
             core=list(fg_params=list(fontsize=6))
           ))
dev.off()

# step2. plot SNP gene pairs and overlap with ATACseq peaks ----
df <- genotype_bed %>%
  dplyr::filter(ID==this_snp) %>%
  select(CHROM, POS, ID) %>%
  mutate(
    start = POS - 1,  # convert to 0-based
    end   = POS
  ) %>%
  select(CHROM, start, end, ID)
fwrite(df, file = paste0(dir,"/transQTL/temp_output/",this_snp,"/SNP_genotype_",ct,"_",condition,"_",QTLtype,".bed"), quote=F, sep="\t")

# step3. SNP genotype correlation with covariates  ----
this_res <- res %>% 
  dplyr::filter(snp==this_snp) %>%
  left_join(. , gene_bed[,c("gene_name","gene_chr")], by=c("transGene"="gene_name"))

this_genotype <- genotype_bed %>% dplyr::filter(ID==this_snp) %>% 
  dplyr::select(-c(CHROM,POS,ID,REF,ALT,MAF)) %>% t() %>%
  as.data.frame() %>%
  rownames_to_column("donor") %>%
  set_colnames(c("donor","genotype"))

# correlate SNP genotype with genotype PC covariates
genotypePC_table <- fread(GENO_PCS_FILE)
df_pc <- left_join(this_genotype, genotypePC_table, by="donor")

# correlate SNP genotype with PEER factors
meta_table <- fread(META_FILE)
peer_table <- fread(PEER_FILE) %>%
  left_join(meta_table, by=c("V1"="sample"))
df_peer <- left_join(this_genotype, peer_table, by="donor") %>% na.omit()

# build correlation table
cor_table <- data.frame(
  covariate = c("genotype_PC1", "genotype_PC2", paste0("PEER", 5:15)),
  correlation = c(
    cor(df_pc$genotype, df_pc$PC1, method="pearson"),
    cor(df_pc$genotype, df_pc$PC2, method="pearson"),
    sapply(paste0("PEER", 5:15), function(p) cor(df_peer$genotype, df_peer[[p]]))
  )
)

print(cor_table)
fwrite(cor_table, file = paste0(dir,"/transQTL/temp_output/",this_snp,"/SNP_genotype_correlation_with_covariates_",ct,"_",condition,"_",QTLtype,".txt"), quote=F, sep="\t")

# step4. make Circos plot per SNP ----
# ── 1. Prepare data ───────────────────────────────────────────────────────────
# Extract cis gene names from SNPtag
cis_genes <- this_res %>%
  distinct(SNPtag) %>%
  mutate(
    cis_gene = str_match(SNPtag,
                         "(?:primary|secondary)_signal_of_([^_]+)_as_")[, 2]
  ) %>%
  filter(!is.na(cis_gene)) %>%
  pull(cis_gene) %>%
  unique()

# Deduplicate: one row per gene (keep primary signal row, drop secondary duplicate)
genes_df <- this_res %>%
  mutate(
    signal = if_else(grepl("primary", cis_trans_category),
                     "primary", "secondary"),
    chr = gene_chr,
    chr_num = as.integer(str_remove(chr, "chr"))
  ) %>%
  dplyr::select(c(snp,transGene,chr,chr_num,q_gene)) %>%
  set_colnames(c("snp","gene","chr","chr_num","q_gene")) %>%
  arrange(chr_num, gene) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(is_cis = gene %in% cis_genes) 

# Add cis genes explicitly onto chr16 if not already in genes_df
# (they may not appear as rows since they are cis, not trans)
snp_chr <- genotype_bed[which(genotype_bed$ID==this_snp),]$CHROM

cis_rows_extra <- tibble(
  gene    = cis_genes,
  chr     = snp_chr,          # chr16
  chr_num = as.numeric(sub("chr","",snp_chr)),
  is_cis  = TRUE,
  q_gene  = NA_real_
) %>%
  filter(!gene %in% genes_df$gene)   # avoid duplicates

genes_df <- bind_rows(genes_df, cis_rows_extra) %>%
  arrange(chr_num, gene)

# Chromosomes present + chr16 even if no trans gene lands there
chrs_in_data <- genes_df %>% pull(chr) %>% unique()
all_chrs_num <- sort(unique(as.integer(str_remove(c(snp_chr,chrs_in_data), "chr"))))
all_chrs     <- paste0("chr", all_chrs_num)

# ── 2. Ring layout ────────────────────────────────────────────────────────────
# Each chromosome gets an arc proportional to gene count (chr16 gets a fixed
# minimum size since it hosts the SNP but may have few/no trans genes)
chr_counts <- genes_df %>%
  count(chr) %>%
  complete(chr = all_chrs, fill = list(n = 0)) %>%
  mutate(
    chr_num   = as.integer(str_remove(chr, "chr")),
    n_display = pmax(n, 1)   # chr16 gets at least 1 unit of arc
  ) %>%
  arrange(chr_num)

total_genes  <- sum(chr_counts$n_display)
gap_deg      <- 3    # degrees of gap between chromosome chunks
n_chrs       <- nrow(chr_counts)
total_gap    <- gap_deg * n_chrs
usable_deg   <- 360 - total_gap

chr_counts <- chr_counts %>%
  mutate(
    arc_deg = n_display / total_genes * usable_deg,
    # cumulative start angle (degrees, 0 = top, clockwise)
    start_deg = cumsum(lag(arc_deg, default = 0)) +
      cumsum(lag(rep(gap_deg, n()), default = 0)),
    end_deg   = start_deg + arc_deg,
    mid_deg   = (start_deg + end_deg) / 2
  )

# ── 3. Helper: degrees → radians, and → x/y on a circle ─────────────────────
deg2rad <- function(d) (d - 90) * pi / 180   # -90 so 0deg = top

ring_xy <- function(deg, r) {
  rad <- deg2rad(deg)
  tibble(x = r * cos(rad), y = r * sin(rad))
}

# ── 4. Build arc polygons for each chromosome ─────────────────────────────────
R_inner <- 5.0
R_outer <- 5.8
R_label <- 6.3   # gene name labels sit here
R_snp   <- 4.3   # SNP node radius (inside the ring)
R_line_start <- 4.95  # lines touch the inner edge of ring

arc_poly <- function(start, end, r_in, r_out, n = 60) {
  angles <- seq(start, end, length.out = n)
  bind_rows(
    map_dfr(angles,         ~ring_xy(.x, r_out)),
    map_dfr(rev(angles),    ~ring_xy(.x, r_in))
  )
}

chr_arcs <- chr_counts %>%
  rowwise() %>%
  mutate(poly = list(arc_poly(start_deg, end_deg, R_inner, R_outer))) %>%
  ungroup()

# ── 5. Chord lines: from SNP (chr16 midpoint) to each trans gene's chr ────────
snp_row   <- chr_counts %>% filter(chr == snp_chr)
snp_angle <- snp_row$mid_deg

# One chord per chromosome (aggregate; could also do one per gene)
chord_df <- genes_df %>%
  left_join(chr_counts %>% dplyr::select(chr, mid_deg, start_deg, end_deg),
            by = "chr") %>%
  dplyr::filter(chr != snp_chr) %>%
  group_by(chr, mid_deg) %>%
  summarise(n_genes = n(), min_q = min(q_gene), .groups = "drop") %>%
  mutate(neg_log_q = -log10(min_q))

# Bezier control point at origin (straight lines through center look clean)
# Use a slight quadratic bezier curving inward
make_bezier <- function(x0, y0, x1, y1, n = 80) {
  # control point pulled toward center (0,0) by 60%
  cx <- (x0 + x1) * 0.15
  cy <- (y0 + y1) * 0.15
  t  <- seq(0, 1, length.out = n)
  tibble(
    x = (1-t)^2 * x0 + 2*(1-t)*t * cx + t^2 * x1,
    y = (1-t)^2 * y0 + 2*(1-t)*t * cy + t^2 * y1
  )
}

snp_xy  <- ring_xy(snp_angle, R_line_start)
chord_paths <- chord_df %>%
  rowwise() %>%
  mutate(
    path = list({
      tgt <- ring_xy(mid_deg, R_line_start)
      make_bezier(snp_xy$x, snp_xy$y, tgt$x, tgt$y) %>%
        mutate(chr = chr, neg_log_q = neg_log_q, n_genes = n_genes)
    })
  ) %>%
  pull(path) %>%
  bind_rows()

# ── 6. Gene labels outside each chromosome chunk ─────────────────────────────
# Helper: wrap a vector of gene names into lines of max n_per_line
# Helper: wrap genes into lines of 5, bolding cis genes
# Uses plotmath-style bold via ggtext <b> tags
wrap_genes_bold <- function(genes, cis_flags, n_per_line = 3) {
  # Tag cis genes with <b>
  tagged <- if_else(cis_flags, paste0("<b>", genes, "</b>"), genes)
  chunks <- split(tagged, ceiling(seq_along(tagged) / n_per_line))
  paste(map_chr(chunks, paste, collapse = ", "), collapse = "<br>")
}

gene_labels <- genes_df %>%
  left_join(chr_counts %>% dplyr::select(chr, mid_deg), by = "chr") %>%
  arrange(chr, desc(is_cis), gene) %>%   # cis genes sort to top of their chr
  group_by(chr, mid_deg) %>%
  summarise(
    label = wrap_genes_bold(gene, is_cis),
    .groups = "drop"
  ) %>%
  mutate(
    rad   = deg2rad(mid_deg),
    x     = R_label * cos(rad),
    y     = R_label * sin(rad),
    angle = mid_deg - 90,
    angle = if_else(mid_deg > 180, angle + 180, angle),
    hjust = if_else(mid_deg <= 180, 0, 1)
  )

# ── 7. Chromosome name labels (just inside outer ring) ────────────────────────
chr_name_labels <- chr_counts %>%
  mutate(
    rad   = deg2rad(mid_deg),
    x     = (R_inner + R_outer) / 2 * cos(rad),
    y     = (R_inner + R_outer) / 2 * sin(rad),
    angle = mid_deg - 90,
    angle = if_else(mid_deg > 180, angle + 180, angle)
  )

# ── 8. Color palette ─────────────────────────────────────────────────────────
chr_palette <- setNames(
  colorRampPalette(c("#185FA5","#0F6E56","#BA7517","#993C1D",
                     "#534AB7","#993556","#3B6D11","#5F5E5A"))(n_chrs),
  chr_counts$chr
)

# ── 9. SNP point ──────────────────────────────────────────────────────────────
snp_point <- ring_xy(snp_angle, R_snp) %>%
  mutate(label = this_snp)

# ── 10. Plot ──────────────────────────────────────────────────────────────────
p <- ggplot() +
  
  # Chromosome arcs
  pmap(list(chr_arcs$poly, chr_arcs$chr), function(poly, chr_id) {
    geom_polygon(data = poly,
                 aes(x = x, y = y),
                 fill  = chr_palette[chr_id],
                 color = "white",
                 linewidth = 0.3,
                 alpha = 0.85)
  }) +
  
  # Chord lines (colored by chromosome, width by n_genes)
  geom_path(
    data = chord_paths,
    aes(x = x, y = y, group = chr,
        color = chr,
        linewidth = n_genes),
    alpha = 0.55
  ) +
  
  # Chromosome name labels (inside arc)
  geom_text(
    data  = chr_name_labels,
    aes(x = x, y = y, label = chr, angle = angle),
    size  = 2.2,
    color = "white",
    fontface = "bold"
  ) +
  
  # Replace the gene label geom with this:
  ggtext::geom_richtext(
    data        = gene_labels,
    aes(x = x, y = y, label = label,
        angle = angle, hjust = hjust,
        color = chr),
    size        = 4,
    lineheight  = 0.9,
    fill        = NA,
    label.size  = 0,
    show.legend = FALSE
  ) +
  
  # SNP node
  geom_point(
    data  = snp_point,
    aes(x = x, y = y),
    shape = 23, size = 5,
    fill  = "#534AB7", color = "white", stroke = 1.2
  ) +
  geom_text(
    data  = snp_point,
    aes(x = x, y = y, label = label),
    vjust = -1, size = 4, fontface = "bold", color = "#534AB7"
  ) +
  
  # Scales
  scale_color_manual(values = chr_palette, guide = "none") +
  scale_linewidth_continuous(
    range  = c(0.5, 3),
    name   = "trans genes\nper chr",
    breaks = c(1, 3, 5, 10)
  ) +
  
  coord_equal(clip = "off") +
  expand_limits(x = c(-9, 9), y = c(-9, 9)) +
  
  labs(
    title    = paste0(this_snp, " trans-eQTL connections by chromosome"),
    subtitle = paste0("SNP on ",snp_chr),
    caption  = "FDR < 0.05 · gene labels grouped by chromosome"
  ) +
  
  theme_void(base_family = "sans") +
  theme(
    plot.title      = element_text(face = "bold", size = 13,
                                   hjust = 3, margin = margin(b = 4)),
    plot.subtitle   = element_text(size = 9, hjust = 3,
                                   color = "grey40", margin = margin(b = 4)),
    plot.caption    = element_text(size = 7, color = "grey60", hjust = 0.5),
    legend.position = "right",
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7),
    plot.margin     = margin(20, 60, 20, 60)   # wide margins for gene labels
  )

ggsave(paste0(dir,"/transQTL/temp_output/",this_snp,"/transQTL_circos_",this_snp,"_",ct,"_",condition,"_",QTLtype,".pdf"),
       p, width = 12, height = 12, device = cairo_pdf)

}
# ------ 6. SNP level analysis #3: pathway analysis --------

# ------ supplemental 1. all tranQTLs: p.value QQ plot and lambda ----------
# Compile all modeling results
dir_results <- paste0(dir, "/transQTL/results/", condition, "/", QTLtype)
files <- list.files(dir_results, pattern = "^result_\\d+\\.tsv$", full.names = TRUE)
modelstats_all <- files %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE))

df <- modelstats_all %>%
  filter(!is.na(p), p > 0, p <= 1)

n <- nrow(df)

df_qq <- df %>%
  arrange(p) %>%
  mutate(
    expected = -log10(ppoints(n)),   # expected under null
    observed = -log10(p)
  )


# ~1.00 → perfect, 1.05–1.1 → mild inflation (often OK in large studies)
# 1.2 → suspicious 🚨
lambda_gc <- median(qchisq(1 - df$p, df = 1)) / qchisq(0.5, df = 1)
lambda_gc

p.qq <- ggplot(df_qq, aes(x = expected, y = observed)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  annotate("text", x = max(df_qq$expected)*0.7, y = max(df_qq$observed)*0.9,
           label = paste0("lambda = ", round(lambda_gc, 3))) +
  labs(
    x = "Expected -log10(p)",
    y = "Observed -log10(p)",
    title = "QQ plot of transQTL p-values"
  ) +
  theme_bw()

png(
  filename = paste0(dir, "/transQTL/results/", condition, "/", QTLtype, "/QQplot_transQTL.png"),
  width = 1800,
  height = 1800,
  res = 300
)
print(p.qq)
dev.off()
rm(df, df_qq, modelstats_all)
gc()

# remove top snps and do it again
top_snps <- res %>%
  count(lead_snp, sort=TRUE) %>%
  pull(lead_snp)

df <- modelstats_all %>%
  filter(!is.na(p), p > 0, p <= 1) %>%
  filter(!snp %in% top_snps)

n <- nrow(df)

lambda_gc <- median(qchisq(1 - df$p, df = 1)) / qchisq(0.5, df = 1)
lambda_gc
