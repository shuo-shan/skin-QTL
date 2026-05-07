#!/usr/bin/env Rscript

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
  library(DESeq2)
})
dir.out = "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/"

#### --------- PART ONE. COMPILE BIG TABLE ------------ ####
# ------- compile all sig cis-QTL SNP table for all celltype x condition ------
sig_snps_all <- data.frame()
for (ct in c("FRB","KRT","MEL")) {
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    for (QTLtype in c("eQTL", "reQTL")) {
      
      if (cond=="PBS" && QTLtype=="reQTL") {
        next
      } 
      
      dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct)
      eigenMT_summary_file <- paste0(dir,"/eigenMT/results/",ct,"_",cond,"_",QTLtype,".eigenMT.txt")
      eigenMT_summary <- fread(eigenMT_summary_file)
      
      # filter by significance threshold
      sig_eigen <- eigenMT_summary %>%
        dplyr::filter(q_gene < 0.05)
      
      # empty check
      if (nrow(sig_eigen) == 0) {
        message(paste0("No significant genes for ", ct, " ", cond, " ", QTLtype))
        next
      }
      
      # sig genes in eigenMT
      sig_genes <- sig_eigen %>%
        dplyr::pull(gene)
      
      # Find max p_gene_eigenMT among significant genes
      cutoff <- sig_eigen %>%
        dplyr::summarise(max_p = max(p_gene_eigenMT)) %>%
        dplyr::pull(max_p)
      
      # sig gene chunks (same for all celltypes)
      if (ct %in% c("KRT", "FRB")) {
        result_chunk_dict <- fread(paste0(dir,"/data/gene_chunk_dict.txt")) %>%
          dplyr::filter(gene %in% sig_genes) %>%
          arrange(chunk)
      } else if (ct == "MEL" ) {
        result_chunk_dict <- fread(paste0(dir,"/data/gene_chunk_dict_old.txt")) %>%
          dplyr::filter(gene %in% sig_genes) %>%
          arrange(chunk)
      }
      chunk_vec <- sort(unique(result_chunk_dict$chunk))
      
      # fetch modeling stats
      sig_snps_across_chunks <- data.frame()
      for (this_chunk in chunk_vec) {
        this_chunk_id <- sprintf("%03d", this_chunk)
        this_chunk_sig_genes <- result_chunk_dict %>%
          dplyr::filter(chunk==this_chunk) %>%
          pull(gene) %>%
          unique()
        message(paste0("processing ", this_chunk_id, " for ", ct," ", cond," ", QTLtype))
        
        modeling_stats_file <- paste0(dir,"/results_QC/modeling_stats_postQC_",ct,"_",cond,"_",QTLtype,"_",this_chunk_id,".txt")
        modeling_stats <- fread(modeling_stats_file) %>%
          dplyr::filter(gene %in% this_chunk_sig_genes)
        
        # fetch significant cis-eQTLs for this chunk
        sig_snps <- modeling_stats %>%
          left_join(eigenMT_summary %>% dplyr::select(gene, Meff), by="gene") %>%
          dplyr::mutate(p_gene_est = p * Meff) %>%
          dplyr::filter(p_gene_est <= cutoff) %>%
          dplyr::mutate(celltype=ct) %>%
          dplyr::mutate(condition=cond) %>%
          dplyr::mutate(QTLtype=QTLtype) %>%
          dplyr::mutate(tag = paste(ct, cond, QTLtype, sep="_")) %>%
          dplyr::select(c(snp, gene, tag, celltype, condition, QTLtype, beta, p)) 
        
        sig_snps_across_chunks <- dplyr::bind_rows(sig_snps_across_chunks, sig_snps)
      }
      sig_snps_all <- dplyr::bind_rows(sig_snps_all, sig_snps_across_chunks)
    }
  }
}

fwrite(sig_snps_all, file = paste0(dir.out,"/sig_snps_all_celltype_conditions.txt"), quote=F, sep="\t")

sig_snps_all %>%
  dplyr::summarise(
    n_snps = dplyr::n(),
    n_genes = dplyr::n_distinct(gene),
    snps_per_gene = dplyr::n() / dplyr::n_distinct(gene),
    .by = c(celltype, condition, QTLtype)
  )


# ------- in bash: create and get bed file for the sig SNPs -------
write(unique(sig_snps_all$snp), file = paste0(dir.out, "/sig_snps_list.txt"))
# in bash, run this:
# workingdir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/
# SNPlist=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/sig_snps_list.txt
# prefix=sig_snps_all_celltype_conditions
# bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/compile_SNP_bed_from_SNP_list.sh ${workingdir} ${SNPlist} ${prefix}
bedFile <- paste0(dir.out,"/QTL_sig_snps_all_celltype_conditions.bed")
sig_snps_bed <- fread(bedFile, sep="\t")
colnames(sig_snps_bed) <- c("chr","start","end","snp","REF","ALT")
length(unique(sig_snps_bed$snp))
length(unique(sig_snps_all$snp))
# awesome, every sig SNP has a row in the bed file

# ------- in bash: create and get ATACseq overlap for the sig SNPs -------
# in bash, run this:
# workingdir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/
# SNPbed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap//QTL_sig_snps_all_celltype_conditions.bed
# prefix=sig_snps_all_celltype_conditions
# bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/annotate_SNP_by_overlapping_ATACseq_peaks.sh ${workingdir} ${SNPbed} ${prefix}
atac_table_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap//sig_snps_all_celltype_conditions_SNP_ATAC_overlap.txt"
atac_table <- fread(atac_table_file) %>%
  dplyr::select(-c(snp_chr, snp_start, snp_end, REF, ALT))
colnames(atac_table)[2:4] <- c("ATACpeak_chr","ATACpeak_start","ATACpeak_end")
sig_snps_all_with_peak <- left_join(sig_snps_all, atac_table, by="snp")


# ------- annotate gene with expression and DE info and potential trans-acting role ------
path.Rscript.fetch_expressed_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_expressed_genes.R"
source(path.Rscript.fetch_expressed_genes)
path.Rscript.fetch_DE_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_DE_genes_log2FC1_padj0.05.R"
source(path.Rscript.fetch_DE_genes)
path.Rscript.fetch_induced_genes="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/fetch_induced_genes_log2FC1_padj0.05.R"
source(path.Rscript.fetch_induced_genes)
trans_acting_genes_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/literature/trans_acting_genes/compiled_trans_acting_candidate_genes_and_category.txt"
trans_acting_genes_table <- fread(trans_acting_genes_file) %>%
  set_colnames(c("gene","gene_category"))

gene_list <- unique(sig_snps_all$gene)

gene_exprs_annotation <- data.frame()
for (ct in c("FRB","KRT","MEL")) {
  for (cond in c("PBS", "IFNB", "IFNG", "TNF")) {
    
    expressed_genes <- unique(fetch_expressed_genes(ct, cond))
    
    if (cond == "PBS") {
      DE_genes <- character(0)
      induced_genes <- character(0)
    } else {
      DE_genes <- unique(fetch_DE_genes(ct, cond))
      induced_genes <- unique(fetch_induced_genes(ct, cond))
    }
    
    this_df <- data.frame(
      gene = gene_list,
      celltype = ct,
      condition = cond,
      is_expressed = gene_list %in% expressed_genes,
      is_DE = if (cond == "PBS") FALSE else gene_list %in% DE_genes,
      is_induced = if (cond == "PBS") FALSE else gene_list %in% induced_genes
    )
    
    gene_exprs_annotation <- dplyr::bind_rows(gene_exprs_annotation, this_df)
  }
}
gene_exprs_annotation <- left_join(gene_exprs_annotation, trans_acting_genes_table, by="gene")
fwrite(gene_exprs_annotation, file = paste0(dir.out,"/sig_eGenes_all_celltype_conditions_expression_DE_annotation_long.txt"), quote=F, sep="\t")

gene_exprs_annotation_wide <- gene_exprs_annotation %>%
  tidyr::pivot_wider(
    names_from = c(celltype, condition),
    values_from = c(is_expressed, is_DE, is_induced)
  )
fwrite(gene_exprs_annotation_wide, file = paste0(dir.out,"/sig_eGenes_all_celltype_conditions_expression_DE_annotation_wide.txt"), quote=F, sep="\t")

# ------- join gene annotation back to sig SNP table with ATACseq peaks ----
sig_snps_big_table <- left_join(sig_snps_all_with_peak, gene_exprs_annotation_wide, by="gene")
fwrite(sig_snps_big_table, file = paste0(dir.out,"/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"), quote=F, sep="\t")

#### --------- PART TWO. How many QTLs have a candidate regulatory region ------------ ####
dir.out = "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap"
df <- read_tsv(paste0(dir.out,"/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))
functional_categories <- unique(df$gene_category)

# ── helper: which "open_celltype_condition" column to use ──────────────────
open_col <- function(ct, cond) paste0("open_", ct, "_", cond)
induced_col <- function(ct) paste0("is_induced_", ct, "_IFNG")  # adjust per condition below
is_expr_col <- function(ct, cond) paste0("is_expressed_", ct, "_", cond)
is_induced_col <- function(ct, cond) paste0("is_induced_", ct, "_", cond)
peak_dynamic_col <- function(ct) paste0("peakDynamic_", ct, "_IFNG")  # gain_accessibility

# ── per-row flags ──────────────────────────────────────────────────────────
# NOTE IMPORTANT: My peakDynamic column currently only includes _IFNG, so for IFNB and TNF, 
#                 the "open peak" and "induced peak" annotation is temporarily using the IFNG dynamic peaks.

# IMPORTANT NOTE 2: I only kept eGenes that are: expressed in PBS, expressed and induced in cytokine
cond_levels <- c("PBS","IFNB","IFNG","TNF")
results <- list()

for (ct in c("FRB","KRT","MEL")) {
  for (cond in c("PBS","IFNG","IFNB","TNF")) {
    for (qt in c("eQTL","reQTL")) {
      
      sub <- df %>% filter(celltype == ct, condition == cond, QTLtype == qt)
      if (nrow(sub) == 0) next
      
      # gene-level filter
      if (cond == "PBS") {
        expr_col <- is_expr_col(ct, "PBS")
        sub <- sub %>% filter(.data[[expr_col]] == TRUE)
      } else {
        ind_col <- is_induced_col(ct, cond)
        if (!ind_col %in% colnames(sub)) next
        sub <- sub %>% filter(.data[[ind_col]] == TRUE)
      }
      if (nrow(sub) == 0) next
      
      # SNP-level flags
      has_peak <- !is.na(sub$peak_name) & sub$peak_name != ""
      
      oc <- open_col(ct, if (cond == "PBS") "PBS" else "IFNG")
      if (!oc %in% colnames(sub)) oc <- open_col(ct, "PBS")
      open_flag <- has_peak & (sub[[oc]] == TRUE)
      
      if (cond != "PBS") {
        pd_col <- paste0("peakDynamic_", ct, "_IFNG")
        if (pd_col %in% colnames(sub)) {
          induced_peak_flag <- has_peak & (sub[[pd_col]] == "gain_accessibility")
        } else {
          induced_peak_flag <- rep(FALSE, nrow(sub))
        }
      } else {
        induced_peak_flag <- rep(FALSE, nrow(sub))
      }
      
      sub <- sub %>%
        mutate(has_peak     = has_peak,
               open_peak    = open_flag,
               induced_peak = induced_peak_flag,
               is_functional = !is.na(gene_category))
      
      gene_sum <- sub %>%
        group_by(gene) %>%
        summarise(
          n_snps              = n(),
          any_open_peak       = any(open_peak),
          any_induced_peak    = any(induced_peak),
          is_functional       = any(!is.na(gene_category)),
          .groups = "drop"
        )
      
      results[[paste(ct, cond, qt, sep="_")]] <- tibble(
        celltype                            = ct,
        condition                           = cond,
        QTLtype                             = qt,
        n_eGenes                            = nrow(gene_sum),
        avg_SNPs_per_gene                   = round(mean(gene_sum$n_snps), 1),
        n_with_open_peak                    = sum(gene_sum$any_open_peak),
        pct_with_open_peak                  = round(sum(gene_sum$any_open_peak)/nrow(gene_sum)*100, 1),
        n_with_induced_peak                 = if (cond != "PBS") sum(gene_sum$any_induced_peak) else NA_integer_,
        pct_with_induced_peak               = if (cond != "PBS") round(sum(gene_sum$any_induced_peak)/nrow(gene_sum)*100,1) else NA_real_,
        n_eGenes_functional                 = sum(gene_sum$is_functional),
        n_eGenes_functional_with_open_peak  = sum(gene_sum$is_functional & gene_sum$any_open_peak),
        n_eGenes_functional_with_induced_peak = if (cond != "PBS") sum(gene_sum$is_functional & gene_sum$any_induced_peak) else NA_integer_
      )
    }
  }
}

result_table <- bind_rows(results)
result_table$condition <- ordered(result_table$condition, levels=c("PBS","IFNB","IFNG","TNF"))
result_table$celltype <- ordered(result_table$celltype, levels=c("KRT","MEL","FRB"))
result_table <- result_table %>% dplyr::arrange(celltype, condition)
fwrite(result_table, paste0(dir.out, "/cis_eQTL_chromatin_summary.txt"), quote=F, sep="\t")

# plot 
plot_df <- result_table %>%
  filter(!is.na(n_eGenes)) %>%
  mutate(
    condition = factor(condition, levels = cond_levels),
    celltype  = factor(celltype, levels = c("FRB","KRT","MEL")),
    bar_label = paste0(celltype, ".", condition, "\n",
                       n_with_open_peak, "/", n_eGenes,
                       "\n(", pct_with_open_peak, "%)")
  )

p <- ggplot(plot_df, aes(x = condition, y = pct_with_open_peak, fill = QTLtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(celltype, ".", condition, "\n",
                               n_with_open_peak, "/", n_eGenes,
                               "\n(", pct_with_open_peak, "%)")),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.5, lineheight = 0.85) +
  facet_wrap(~celltype, scales = "free_x") +
  scale_fill_manual(values = c("eQTL" = "#4472C4", "reQTL" = "#ED7D31")) +
  scale_y_continuous(limits = c(0, 115), expand = c(0,0)) +
  scale_x_discrete(limits = cond_levels) +
  labs(
    title = "eGenes with ≥1 sig SNP overlapping an open ATAC peak \n IFNB and TNF don't have ATACseq info, used PBS openness instead)",
    x = NULL, y = "% eGenes with open peak overlap",
    fill = "QTL type"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text   = element_text(face = "bold"))

ggsave(paste0(dir.out, "/cis_eQTL_chromatin_barplot.pdf"), p, width = 12, height = 6)

# 
plot_df_induced <- result_table %>%
  filter(!is.na(n_with_induced_peak)) %>%  # remove PBS
  mutate(
    condition = factor(condition, levels = cond_levels),
    celltype  = factor(celltype,  levels = c("FRB","KRT","MEL"))
  )

p_induced <- ggplot(plot_df_induced, 
                    aes(x = condition, y = pct_with_induced_peak, fill = QTLtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(celltype, ".", condition, "\n",
                               n_with_induced_peak, "/", n_eGenes,
                               "\n(", pct_with_induced_peak, "%)")),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.5, lineheight = 0.85) +
  facet_wrap(~celltype, scales = "free_x") +
  scale_fill_manual(values = c("eQTL" = "#4472C4", "reQTL" = "#ED7D31")) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  scale_x_discrete(limits = cond_levels) +
  labs(
    title = "eGenes with ≥1 sig SNP overlapping an INDUCED ATAC peak",
    subtitle = "\n IFNB and TNF don't have ATACseq info, used IFNG induction instead \nInduced = gain_accessibility upon cytokine stimulation",
    x    = NULL,
    y    = "% eGenes with induced peak overlap",
    fill = "QTL type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(paste0(dir.out, "/cis_eQTL_chromatin_barplot_induced.pdf"), p_induced, width = 12, height = 6)


#### --------- PART THREE. Do genes with shared QTLs across cell types have a common peak? ---------
df <- read_tsv(paste0(dir.out, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))

# Step 1: for each SNP x gene pair x condition x QTLtype, check in which celltype(s) it's sig eQTL.
snp_celltype <- df %>%
  group_by(snp, gene, condition, QTLtype) %>%
  summarise(
    celltypes_with_QTL = paste(sort(unique(celltype)), collapse=","),
    n_celltypes        = n_distinct(celltype),
    any_open_FRB = any(open_FRB_PBS == TRUE | open_FRB_IFNG == TRUE, na.rm=TRUE),
    any_open_KRT = any(open_KRT_PBS == TRUE | open_KRT_IFNG == TRUE, na.rm=TRUE),
    any_open_MEL = any(open_MEL_PBS == TRUE | open_MEL_IFNG == TRUE, na.rm=TRUE),
    has_peak     = any(!is.na(peak_name)),
    .groups = "drop"
  ) %>%
  mutate(
    is_shared_QTL    = n_celltypes > 1,
    n_open_celltypes = any_open_FRB + any_open_KRT + any_open_MEL,
    has_common_peak  = has_peak & (n_open_celltypes >= n_celltypes)
  )

# Step 2: gene-level summary per condition x QTLtype
gene_summary <- snp_celltype %>%
  group_by(gene, condition, QTLtype) %>%
  summarise(
    n_celltypes_with_QTL = max(n_celltypes),
    is_shared_QTL        = any(is_shared_QTL),
    any_common_peak      = any(has_common_peak),
    .groups = "drop"
  )

# Step 3: summarize per condition x QTLtype
result2 <- gene_summary %>%
  group_by(condition, QTLtype, is_shared_QTL) %>%
  summarise(
    n_genes            = n(),
    n_with_common_peak = sum(any_common_peak),
    pct                = round(n_with_common_peak/n_genes*100, 1),
    .groups = "drop"
  )

# Step 4: Fisher's test per condition x QTLtype
fisher_results <- result2 %>%
  group_by(condition, QTLtype) %>%
  summarise({
    shared     <- cur_data() %>% filter(is_shared_QTL == TRUE)
    specific   <- cur_data() %>% filter(is_shared_QTL == FALSE)
    if (nrow(shared) == 0 | nrow(specific) == 0) return(tibble(OR=NA, p=NA))
    mat <- matrix(c(
      shared$n_with_common_peak,
      shared$n_genes - shared$n_with_common_peak,
      specific$n_with_common_peak,
      specific$n_genes - specific$n_with_common_peak
    ), nrow=2)
    ft <- fisher.test(mat)
    tibble(OR = round(ft$estimate, 1), p = ft$p.value)
  }, .groups = "drop")

# Step 5: plot
cond_levels <- c("PBS","IFNB","IFNG","TNF")

plot_df <- result2 %>%
  mutate(
    condition     = factor(condition, levels = cond_levels),
    QTL_category  = ifelse(is_shared_QTL, "Shared QTL\n(>1 celltype)", "Celltype-specific QTL"),
    QTL_category  = factor(QTL_category, levels = c("Celltype-specific QTL","Shared QTL\n(>1 celltype)"))
  )

annot_df <- fisher_results %>%
  mutate(
    condition = factor(condition, levels = cond_levels),
    label     = paste0("OR=", OR, "\np=", formatC(p, format="e", digits=1))
  )

p2 <- ggplot(plot_df, aes(x = QTL_category, y = pct, fill = QTL_category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(n_with_common_peak, "/", n_genes, "\n(", pct, "%)")),
            vjust = -0.4, size = 3) +
  geom_text(data = annot_df,
            aes(x = 1.5, y = 108, label = label),
            size = 2.8, fontface = "italic", inherit.aes = FALSE) +
  facet_grid(QTLtype ~ condition) +
  scale_fill_manual(values = c(
    "Celltype-specific QTL"     = "#7aaad4",
    "Shared QTL\n(>1 celltype)" = "#1a5ea8"
  )) +
  scale_y_continuous(limits = c(0, 120), expand = c(0,0)) +
  labs(
    title = "Genes with shared QTLs are enriched for common open chromatin peaks",
    subtitle = "Stratified by condition and QTL type",
    x    = NULL,
    y    = "% eGenes with common open peak",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

ggsave(paste0(dir.out, "/question2_shared_QTL_common_peak_by_condition.pdf"), p2, width = 10, height = 7)

#### --------- PART FOUR. IFNB eGenes overlap with IFNG/TNF eGenes? ----
df <- read_tsv(paste0(dir.out, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))

# for each gene x celltype，which in which cytokine conditions there's sig eQTL
# for each celltype x condition, only keep the gene that's induced in that condition.
cytokine_overlap <- df %>%
  filter(condition %in% c("IFNB","IFNG","TNF")) %>%
  filter(
    (condition == "IFNB" & is_induced_FRB_IFNB == TRUE & celltype == "FRB") |
      (condition == "IFNB" & is_induced_KRT_IFNB == TRUE & celltype == "KRT") |
      (condition == "IFNB" & is_induced_MEL_IFNB == TRUE & celltype == "MEL") |
      (condition == "IFNG" & is_induced_FRB_IFNG == TRUE & celltype == "FRB") |
      (condition == "IFNG" & is_induced_KRT_IFNG == TRUE & celltype == "KRT") |
      (condition == "IFNG" & is_induced_MEL_IFNG == TRUE & celltype == "MEL") |
      (condition == "TNF"  & is_induced_FRB_TNF  == TRUE & celltype == "FRB") |
      (condition == "TNF"  & is_induced_KRT_TNF  == TRUE & celltype == "KRT") |
      (condition == "TNF"  & is_induced_MEL_TNF  == TRUE & celltype == "MEL")
  ) %>%
  group_by(gene, celltype) %>%
  summarise(
    has_IFNB = "IFNB" %in% unique(condition),
    has_IFNG = "IFNG" %in% unique(condition),
    has_TNF  = "TNF"  %in% unique(condition),
    .groups  = "drop"
  )

# out of IFNB eGenes, how many of then are IFNG or TNF eGenes?
ifnb_genes <- cytokine_overlap %>%
  filter(has_IFNB == TRUE) %>%
  mutate(
    also_IFNG = has_IFNG,
    also_TNF  = has_TNF,
    IFNB_only = !has_IFNG & !has_TNF
  )

ifnb_genes %>%
  group_by(celltype) %>%
  summarise(
    n_IFNB_eGenes   = n(),
    n_also_IFNG     = sum(also_IFNG),
    n_also_TNF      = sum(also_TNF),
    n_IFNB_only     = sum(IFNB_only),
    pct_IFNB_only   = round(sum(IFNB_only)/n()*100, 1)
  )

# check each cytokine's specificity
# IFNB specificity
ifnb_summary <- cytokine_overlap %>%
  filter(has_IFNB == TRUE) %>%
  group_by(celltype) %>%
  summarise(
    n_eGenes      = n(),
    n_also_IFNG   = sum(has_IFNG),
    n_also_TNF    = sum(has_TNF),
    n_specific    = sum(!has_IFNG & !has_TNF),
    pct_specific  = round(n_specific/n()*100, 1),
    cytokine      = "IFNB"
  )

# IFNG specificity
ifng_summary <- cytokine_overlap %>%
  filter(has_IFNG == TRUE) %>%
  group_by(celltype) %>%
  summarise(
    n_eGenes      = n(),
    n_also_IFNB   = sum(has_IFNB),
    n_also_TNF    = sum(has_TNF),
    n_specific    = sum(!has_IFNB & !has_TNF),
    pct_specific  = round(n_specific/n()*100, 1),
    cytokine      = "IFNG"
  )

# TNF specificity
tnf_summary <- cytokine_overlap %>%
  filter(has_TNF == TRUE) %>%
  group_by(celltype) %>%
  summarise(
    n_eGenes      = n(),
    n_also_IFNB   = sum(has_IFNB),
    n_also_IFNG   = sum(has_IFNG),
    n_specific    = sum(!has_IFNB & !has_IFNG),
    pct_specific  = round(n_specific/n()*100, 1),
    cytokine      = "TNF"
  )

# combined
bind_rows(
  ifnb_summary %>% select(celltype, cytokine, n_eGenes, n_specific, pct_specific),
  ifng_summary %>% select(celltype, cytokine, n_eGenes, n_specific, pct_specific),
  tnf_summary  %>% select(celltype, cytokine, n_eGenes, n_specific, pct_specific)
) %>% arrange(celltype, cytokine)


##### ----- PART FIVE. log2FC info of induced eGenes -----
df <- read_tsv(paste0(dir.out, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))
load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData")
vars <- ls(pattern = "^deg")

# compile all DEGs
deg_all <- bind_rows(
  deg.frb.ifnb %>% mutate(celltype = "FRB", condition = "IFNB"),
  deg.frb.ifng %>% mutate(celltype = "FRB", condition = "IFNG"),
  deg.frb.tnf  %>% mutate(celltype = "FRB", condition = "TNF"),
  deg.krt.ifnb %>% mutate(celltype = "KRT", condition = "IFNB"),
  deg.krt.ifng %>% mutate(celltype = "KRT", condition = "IFNG"),
  deg.krt.tnf  %>% mutate(celltype = "KRT", condition = "TNF"),
  deg.mel.ifnb %>% mutate(celltype = "MEL", condition = "IFNB"),
  deg.mel.ifng %>% mutate(celltype = "MEL", condition = "IFNG"),
  deg.mel.tnf  %>% mutate(celltype = "MEL", condition = "TNF")
)

# eQTL genes per celltype x condition (no PBS)
eqtl_genes <- df %>%
  filter(condition != "PBS") %>%
  distinct(gene, celltype, condition)

# overlap
plot_df <- deg_all %>%
  inner_join(eqtl_genes %>% mutate(has_eQTL = TRUE),
             by = c("gene","celltype","condition")) %>%
  mutate(condition = factor(condition, levels = c("IFNB","IFNG","TNF")))

# plot
p <- ggplot(plot_df, aes(x = condition, y = log2FoldChange, fill = condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~celltype) +
  scale_fill_manual(values = c(
    "IFNB" = "#4472C4",
    "IFNG" = "#ED7D31",
    "TNF"  = "#70AD47"
  )) +
  labs(
    title    = "log2FoldChange of induced DEGs by cytokine",
    subtitle = "All induced DEGs; eQTL genes overlaid",
    x        = NULL,
    y        = "log2FoldChange",
    fill     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )

ggsave(paste0(dir.out, "/IFNB_foldchange_comparison.pdf"), p, width = 10, height = 5)


deg_all %>%
  filter(log2FoldChange > 0) %>%  # 只看induced
  group_by(celltype, condition) %>%
  summarise(n_induced = n_distinct(gene))


###### ---- PART SIX. could it be IFNB eQTLs are lower in QTL-beta compared to other cytokines? -------
beta_df <- df %>%
  filter(condition %in% c("IFNB","IFNG","TNF")) %>%
  distinct(snp, gene, celltype, condition, beta) %>%
  mutate(
    abs_beta  = abs(beta),
    condition = factor(condition, levels = c("IFNB","IFNG","TNF"))
  )
n_egenes <- beta_df %>%
  group_by(celltype, condition) %>%
  summarise(n_eGenes = n_distinct(gene), .groups = "drop") %>%
  mutate(label = paste0("n=", n_eGenes))

p_beta <- ggplot(beta_df, aes(x = condition, y = abs_beta, fill = condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_text(data = n_egenes,
            aes(x = condition, y = -0.05, label = label),
            size = 3, inherit.aes = FALSE) +
  facet_wrap(~celltype) +
  scale_fill_manual(values = c(
    "IFNB" = "#4472C4",
    "IFNG" = "#ED7D31",
    "TNF"  = "#70AD47"
  )) +
  scale_y_continuous(limits = c(-0.1, 3)) +
  labs(
    title = "eQTL effect size (|beta|) by cytokine condition",
    x     = NULL,
    y     = "|beta|",
    fill  = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )


ggsave(paste0(dir.out, "/eQTL_beta_by_condition.pdf"), p_beta, width = 10, height = 5)


###### ---- PART SEVEN. could it be IFNB eQTLs have more selection pressure than other cytokines? ----
# in bash, run this:
# workingdir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/
# SNPbed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap//QTL_sig_snps_all_celltype_conditions.bed
# prefix=sig_snps_all_celltype_conditions
# bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/annotate_SNP_by_overlapping_ATACseq_peaks.sh ${workingdir} ${SNPbed} ${prefix}
phyloP_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/cis_eQTL_ATAC_overlap/phyloP_scores_sig_snps_all_celltype_conditions.txt"
plyloP_table <- fread(phyloP_file, header = F) %>%
  set_colnames(c("snp","size","bp_covered_by_phyloP","phyloP_sum","phyloP_mean_in_snp","phyloP_mean_in_covered")) %>%
  dplyr::select(c(snp, phyloP_mean_in_covered)) %>%
  set_colnames(c("snp", "phyloP"))

sig_snps_big_table_with_phyloP <- left_join(sig_snps_big_table, plyloP_table, by="snp")

# filter: cytokine only, induced genes only
phyloP_plot_df <- sig_snps_big_table_with_phyloP %>%
  filter(condition != "PBS") %>%
  filter(
    (condition == "IFNB" & celltype == "FRB" & is_induced_FRB_IFNB == TRUE) |
      (condition == "IFNB" & celltype == "KRT" & is_induced_KRT_IFNB == TRUE) |
      (condition == "IFNB" & celltype == "MEL" & is_induced_MEL_IFNB == TRUE) |
      (condition == "IFNG" & celltype == "FRB" & is_induced_FRB_IFNG == TRUE) |
      (condition == "IFNG" & celltype == "KRT" & is_induced_KRT_IFNG == TRUE) |
      (condition == "IFNG" & celltype == "MEL" & is_induced_MEL_IFNG == TRUE) |
      (condition == "TNF"  & celltype == "FRB" & is_induced_FRB_TNF  == TRUE) |
      (condition == "TNF"  & celltype == "KRT" & is_induced_KRT_TNF  == TRUE) |
      (condition == "TNF"  & celltype == "MEL" & is_induced_MEL_TNF  == TRUE)
  ) %>%
  filter(!is.na(phyloP)) %>%
  mutate(condition = factor(condition, levels = c("IFNB","IFNG","TNF")),
         celltype  = factor(celltype,  levels = c("FRB","KRT","MEL")))

# n SNPs per celltype x condition for annotation
n_snps <- phyloP_plot_df %>%
  group_by(celltype, condition) %>%
  summarise(n = n_distinct(snp), .groups = "drop") %>%
  mutate(label = paste0("n=", n))

# wilcoxon test IFNB vs IFNG and IFNB vs TNF per celltype
stats_df <- phyloP_plot_df %>%
  group_by(celltype) %>%
  summarise(
    p_IFNB_vs_IFNG = wilcox.test(
      phyloP[condition == "IFNB"],
      phyloP[condition == "IFNG"])$p.value,
    p_IFNB_vs_TNF = wilcox.test(
      phyloP[condition == "IFNB"],
      phyloP[condition == "TNF"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0("IFNB vs IFNG: p=", formatC(p_IFNB_vs_IFNG, format="e", digits=1),
                   "\nIFNB vs TNF: p=",  formatC(p_IFNB_vs_TNF,  format="e", digits=1))
  )

# plot
p_phyloP <- ggplot(phyloP_plot_df, 
                   aes(x = condition, y = phyloP, fill = condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = n_snps,
            aes(x = condition, y = -8, label = label),
            size = 3, inherit.aes = FALSE) +
  geom_text(data = stats_df,
            aes(x = 2, y = 7, label = label),
            size = 2.8, fontface = "italic", inherit.aes = FALSE) +
  facet_wrap(~celltype) +
  scale_fill_manual(values = c(
    "IFNB" = "#4472C4",
    "IFNG" = "#ED7D31",
    "TNF"  = "#70AD47"
  )) +
  scale_y_continuous(limits = c(-9, 8)) +
  labs(
    title    = "phyloP20way conservation score of eQTL SNPs by cytokine",
    subtitle = "Only induced eGenes per condition; higher = more conserved",
    x        = NULL,
    y        = "phyloP20way score",
    fill     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none"
  )

ggsave(paste0(dir.out, "/phyloP_by_condition.pdf"), p_phyloP, width = 10, height = 5)

###### ---- IFNB low eGene issue: log2FC violin plot ####
### load data ----
df <- read_tsv(paste0(dir.out, "/cisQTL_sigSNPs_ATACpeak_exprDEinduced_annotated.txt"))
load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData")
vars <- ls(pattern = "^deg")

# compile all DEGs
deg_all <- bind_rows(
  deg.frb.ifnb %>% mutate(celltype = "FRB", condition = "IFNB"),
  deg.frb.ifng %>% mutate(celltype = "FRB", condition = "IFNG"),
  deg.frb.tnf  %>% mutate(celltype = "FRB", condition = "TNF"),
  deg.krt.ifnb %>% mutate(celltype = "KRT", condition = "IFNB"),
  deg.krt.ifng %>% mutate(celltype = "KRT", condition = "IFNG"),
  deg.krt.tnf  %>% mutate(celltype = "KRT", condition = "TNF"),
  deg.mel.ifnb %>% mutate(celltype = "MEL", condition = "IFNB"),
  deg.mel.ifng %>% mutate(celltype = "MEL", condition = "IFNG"),
  deg.mel.tnf  %>% mutate(celltype = "MEL", condition = "TNF")
)
rm(deg.frb.ifnb, deg.frb.ifng, deg.frb.tnf, 
   deg.krt.ifnb, deg.krt.ifng, deg.krt.tnf,
   deg.mel.ifnb, deg.mel.ifng, deg.mel.tnf)

load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_DESeq2_results_post_outlier_exclusion.RData")
vars <- ls(pattern = "^res.")

# compile all DESeq2 results
deseq_res_all <- bind_rows(
  as.data.frame(res.frb.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype = "FRB", condition = "IFNB"),
  as.data.frame(res.frb.ifng) %>% rownames_to_column("gene") %>% mutate(celltype = "FRB", condition = "IFNG"),
  as.data.frame(res.frb.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype = "FRB", condition = "TNF"),
  as.data.frame(res.krt.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype = "KRT", condition = "IFNB"),
  as.data.frame(res.krt.ifng) %>% rownames_to_column("gene") %>% mutate(celltype = "KRT", condition = "IFNG"),
  as.data.frame(res.krt.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype = "KRT", condition = "TNF"),
  as.data.frame(res.mel.ifnb) %>% rownames_to_column("gene") %>% mutate(celltype = "MEL", condition = "IFNB"),
  as.data.frame(res.mel.ifng) %>% rownames_to_column("gene") %>% mutate(celltype = "MEL", condition = "IFNG"),
  as.data.frame(res.mel.tnf)  %>% rownames_to_column("gene") %>% mutate(celltype = "MEL", condition = "TNF")
)
rm(res.frb.ifnb, res.frb.ifng, res.frb.tnf, 
   res.krt.ifnb, res.krt.ifng, res.krt.tnf,
   res.mel.ifnb, res.mel.ifng, res.mel.tnf)


###  [>>> anchor on IFNG eGenes <<<] ----
# QUESTION 1: Are IFNG eGenes also induced by IFNB? 
celltypes <- c("FRB", "MEL", "KRT")

gene_table_list <- map(celltypes, function(ct) {
  
  # IFNG eGenes for this celltype
  ifng_egenes <- df %>%
    filter(celltype == ct, condition == "IFNG", QTLtype == "eQTL") %>%
    distinct(gene)
  
  # log2FC from deg_all for IFNG and IFNB
  deg_ct <- deg_all %>%
    filter(celltype == ct, condition %in% c("IFNG", "IFNB")) %>%
    select(gene, condition, log2FoldChange)
  
  ifng_egenes %>%
    left_join(deg_ct %>% filter(condition == "IFNG") %>%
                select(gene, log2FC_IFNG = log2FoldChange),
              by = "gene") %>%
    left_join(deg_ct %>% filter(condition == "IFNB") %>%
                select(gene, log2FC_IFNB = log2FoldChange),
              by = "gene") %>%
    mutate(
      celltype = ct,
      induced_in_IFNG = !is.na(log2FC_IFNG),
      induced_in_IFNB = !is.na(log2FC_IFNB)
    )
}) %>%
  set_names(celltypes)

# Quick summary for Q1
map_dfr(celltypes, function(ct) {
  
  ifng_induced <- deg_all %>%
    filter(celltype == ct, condition == "IFNG", padj < 0.05, log2FoldChange > 1.5) %>%
    distinct(gene)
  
  ifng_eg <- df %>%
    filter(celltype == ct, condition == "IFNG", QTLtype == "eQTL") %>%
    distinct(gene)
  
  ifnb_induced <- deg_all %>%
    filter(celltype == ct, condition == "IFNB", padj < 0.05, log2FoldChange > 1.5) %>%
    distinct(gene)
   
  # IFNG eGenes that are also IFNG induced
  ifng_eg_and_induced <- ifng_eg %>%
    filter(gene %in% ifng_induced$gene)
  
  tibble(
    celltype = ct,
    n_IFNG_induced = nrow(ifng_induced),
    n_IFNG_eGenes = nrow(ifng_eg),
    n_IFNG_induced_eGene = nrow(ifng_eg_and_induced),
    # NEW: of those, how many also induced in IFNB?
    n_IFNG_induced_eGene_IFNBinduced = sum(ifng_eg_and_induced$gene %in% ifnb_induced$gene),
    # and as percentage
    pct_also_IFNB = round(100 * sum(ifng_eg_and_induced$gene %in% ifnb_induced$gene) / nrow(ifng_eg_and_induced), 1)
  )
})

# ============================================================
# SETUP: define gene sets and log2FC data
# ============================================================

cytokine_colors <- c(
  "IFNG_sig"     = "#1E8449",   # dark green
  "IFNB_sig"     = "#C0392B",   # dark red
  "IFNB_not_sig" = "#F1948A",   # pink
  "TNF_sig"      = "#1A5276",   # dark blue
  "TNF_not_sig"  = "#85C1E9"    # light blue
)

make_plot_data <- function(anchor_genes_fn) {
  map_dfr(celltypes, function(ct) {
    
    anchor <- anchor_genes_fn(ct)
    
    ifnb_sig_genes <- deg_all %>%
      filter(celltype == ct, condition == "IFNB",
             padj < 0.05, log2FoldChange > 1.5) %>%
      pull(gene)
    
    tnf_sig_genes <- deg_all %>%
      filter(celltype == ct, condition == "TNF",
             padj < 0.05, log2FoldChange > 1.5) %>%
      pull(gene)
    
    deseq_res_all %>%
      filter(celltype == ct,
             condition %in% c("IFNG", "IFNB", "TNF"),  # added IFNG
             gene %in% anchor) %>%
      select(gene, condition, log2FoldChange) %>%
      mutate(
        celltype = ct,
        sig_group = case_when(
          condition == "IFNG"  ~ "IFNG_sig",            # IFNG always sig (it's the anchor)
          condition == "IFNB" & gene %in% ifnb_sig_genes ~ "IFNB_sig",
          condition == "IFNB" & !gene %in% ifnb_sig_genes ~ "IFNB_not_sig",
          condition == "TNF"  & gene %in% tnf_sig_genes  ~ "TNF_sig",
          condition == "TNF"  & !gene %in% tnf_sig_genes ~ "TNF_not_sig"
        )
      )
  })
}

# ============================================================
# ANCHOR 1: All IFNG induced genes
# ============================================================
anchor_ifng_induced <- function(ct) {
  deg_all %>%
    dplyr::filter(celltype == ct, condition == "IFNG",
           padj < 0.05, log2FoldChange > 1.5) %>%
    pull(gene)
}

# ============================================================
# ANCHOR 2: IFNG induced genes that are also IFNG eGenes
# ============================================================
anchor_ifng_induced_egene <- function(ct) {
  ifng_induced <- deg_all %>%
    dplyr::filter(celltype == ct, condition == "IFNG",
           padj < 0.05, log2FoldChange > 1.5) %>%
    pull(gene)
  
  ifng_egenes <- df %>%
    dplyr::filter(celltype == ct, condition == "IFNG", QTLtype == "eQTL") %>%
    pull(gene)
  
  intersect(ifng_induced, ifng_egenes)
}

plot_data_induced      <- make_plot_data(anchor_ifng_induced)
plot_data_induced_egene <- make_plot_data(anchor_ifng_induced_egene)

# ============================================================
# PLOT FUNCTIONS
# ============================================================

make_violin <- function(data, title) {
  
  # calculate n for each group
  n_labels <- data %>%
    ungroup() %>%
    dplyr::count(sig_group) %>%
    mutate(
      sig_group = factor(sig_group, levels = names(cytokine_colors)),
      label = paste0("n=", n)
    )
  
  data %>%
    ungroup() %>%
    mutate(sig_group = factor(sig_group, levels = names(cytokine_colors))) %>%
    ggplot(aes(x = sig_group, y = log2FoldChange, fill = sig_group)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1.5, linetype = "dotted", color = "gray50") +
    geom_text(data = n_labels,
              aes(x = sig_group, label = label),
              y = max(data$log2FoldChange, na.rm = TRUE) * 1.05,
              inherit.aes = FALSE,
              size = 3.5, color = "gray30") +
    scale_fill_manual(values = cytokine_colors) +
    labs(title = title,
         x = NULL, y = "log2FoldChange", fill = NULL) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 20, hjust = 1))
}

make_scatter <- function(data, compare_condition, title) {
  
  scatter_data <- data %>%
    dplyr::ungroup() %>%
    dplyr::filter(condition %in% c("IFNG", compare_condition)) %>%
    dplyr::select(gene, celltype, condition, log2FoldChange, sig_group) %>%
    tidyr::pivot_wider(names_from = condition,
                       values_from = c(log2FoldChange, sig_group)) %>%
    dplyr::rename(
      log2FC_IFNG       = log2FoldChange_IFNG,
      log2FC_compare    = !!sym(paste0("log2FoldChange_", compare_condition)),
      sig_group_compare = !!sym(paste0("sig_group_", compare_condition))
    )
  
  # unified axis limits
  axis_min <- min(c(scatter_data$log2FC_IFNG,
                    scatter_data$log2FC_compare), na.rm = TRUE) * 1.05
  axis_max <- max(abs(c(scatter_data$log2FC_IFNG,
                        scatter_data$log2FC_compare)), na.rm = TRUE) * 1.05
  axis_lim <- c(axis_min, axis_max)
  
  scatter_data %>%
    dplyr::mutate(sig_group_compare = factor(sig_group_compare,
                                             levels = names(cytokine_colors))) %>%
    ggplot(aes(x = log2FC_IFNG, y = log2FC_compare,
               color = sig_group_compare)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 1.5, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 1.5, linetype = "dotted", color = "gray50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = "gray30", alpha = 0.5) +  # y=x reference line
    scale_color_manual(values = cytokine_colors, drop = FALSE) +
    scale_x_continuous(limits = axis_lim) +
    scale_y_continuous(limits = axis_lim) +
    labs(
      title = title,
      x = "log2FC in IFNG",
      y = paste0("log2FC in ", compare_condition),
      color = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# ============================================================
# SAVE PDFs — one celltype per page
# ============================================================

save_pdf <- function(plot_data, plot_fn, filename, title_prefix) {
  pdf(filename, width = 8, height = 6)
  for (ct in celltypes) {
    p <- plot_fn(
      plot_data %>% filter(celltype == ct),
      title = paste0(title_prefix, " — ", ct)
    )
    print(p)
  }
  dev.off()
  cat("Saved:", filename, "\n")
}

# Violin plots
save_pdf(plot_data_induced,
         make_violin,
         paste0(dir.out,"/violin_IFNG_induced_genes.pdf"),
         "IFNG induced genes: log2FC in IFNB and TNF")

save_pdf(plot_data_induced_egene,
         make_violin,
         paste0(dir.out,"violin_IFNG_induced_eGenes.pdf"),
         "IFNG induced eGenes: log2FC in IFNB and TNF")

# scatter
save_scatter_pdf <- function(plot_data, filename, title_prefix) {
  pdf(filename, width = 7, height = 6)
  for (ct in celltypes) {
    for (cond in c("IFNB", "TNF")) {
      p <- make_scatter(
        data = plot_data %>% dplyr::filter(celltype == ct),
        compare_condition = cond,
        title = paste0(title_prefix, " — ", ct, " vs ", cond)
      )
      print(p)
    }
  }
  dev.off()
  cat("Saved:", filename, "\n")
}

# Save scatter PDFs
save_scatter_pdf(plot_data_induced,
                 paste0(dir.out, "scatter_IFNG_induced_genes.pdf"),
                 "IFNG induced genes: IFNG vs")

save_scatter_pdf(plot_data_induced_egene,
                 paste0(dir.out, "/scatter_IFNG_induced_eGenes.pdf"),
                 "IFNG induced eGenes: IFNG vs")
