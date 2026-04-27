#!/usr/bin/env Rscript
# Create overview tables and graphs for reclassify jobs into one table per celltype.
# Run after all BSUB jobs finish.
#

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(future.apply)
  library(magrittr)
  library(purrr)
  library(corpcor)
  library(ggplot2)
  library(qpdf)
})

ct <- "KRT"
dir <- paste0("~/Downloads/nl/human/skin/eQTLs/DREG/2025models/", ct)

# ---- Load reclassify tables ---- 
f.ifnb <- data.table::fread(paste0(dir,"/reclassified/reclassify_summary_IFNB.txt")) 
f.ifng <- data.table::fread(paste0(dir,"/reclassified/reclassify_summary_IFNG.txt")) 
f.tnf <- data.table::fread(paste0(dir,"/reclassified/reclassify_summary_TNF.txt")) 
f <- rbind(f.ifnb, f.ifng, f.tnf) %>%
  mutate(
    abs_delta  = abs(anchorSNP_cytokine_beta - anchorSNP_PBS_beta),
    # anchor condition: whichever has the lower p-value
    anchor_beta = ifelse(anchorSNP_cytokine_p <= anchorSNP_PBS_p, abs(anchorSNP_cytokine_beta), abs(anchorSNP_PBS_beta)),
    anchor_p    = pmin(anchorSNP_cytokine_p, anchorSNP_PBS_p),
    
    priority_score = case_when(
      str_starts(gene_class, "constitutive") ~
        -log10(anchor_p) * anchor_beta * abs_delta,
      
      str_starts(gene_class, "emergent") ~
        -log10(anchorSNP_cytokine_p) * abs(anchorSNP_cytokine_beta) * abs_delta,
      
      str_starts(gene_class, "vanishing") ~
        -log10(anchorSNP_PBS_p) * abs(anchorSNP_PBS_beta) * abs_delta,
      
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(-abs_delta, -anchor_beta, -anchor_p)
cytokines <- c("IFNB","IFNG","TNF")

# ---- Load DE genes ---- 
load("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_degs_abslog2FC1_padj0.05_post_outlier_exclusion.RData")

# prep data table ####
vars_deg <- ls(pattern = "^deg")

# functions ####
# load meta genes to gene name dictionary
metagene_dict <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/literature/metaidname.txt") %>% 
  dplyr::select(-id) %>% distinct()
metagene_dict_collapsed <- metagene_dict %>% group_by(meta) %>% 
  summarize(source=paste(name,collapse=','))

convert_meta_genes <- function(gene_vector) {
  this_metagenes_idx <- which(grepl("meta", gene_vector))
  temp1 <- data.frame(meta=gene_vector[-this_metagenes_idx], 
                      name=gene_vector[-this_metagenes_idx])
  temp2 <- data.frame(meta=gene_vector[this_metagenes_idx]) %>%
    left_join(., metagene_dict, join_by(meta))
  temp3 <- rbind(temp1, temp2) %>% set_colnames(c("meta","gene"))
  
  return(temp3)
}


# make one small df per object: gene + logical column
deg_dfs <- lapply(vars_deg, function(var) {
  data.frame(
    gene = get(var)$gene,
    value = TRUE,
    stringsAsFactors = FALSE
  ) %>%
    distinct(gene, .keep_all = TRUE) %>%
    rename(!!paste0("is_", var) := value)
})

# merge into one wide table
deg_wide <- purrr::reduce(deg_dfs, full_join, by = "gene") %>%
  mutate(across(starts_with("is_deg."), ~tidyr::replace_na(., FALSE)))


# ---- Count how many DE genes per celltype x condition ----
ct="krt"
stim="ifnb"
var <- ls(pattern = paste0("^deg.",ct,".",stim))
this_deg <- get(var)
colnames(this_deg)[1] <- "meta"
this_meta_dict <- convert_meta_genes(this_deg$meta)
this_deg_metaconverted <- left_join(this_meta_dict, this_deg, by="meta")

# up-regulated genes
deg.up <- this_deg_metaconverted[this_deg_metaconverted$log2FoldChange > 0,]
length(unique(deg.up$gene))

this_qtlgene <- f_unambiguous_genes_all_cyto %>%
  dplyr::filter(celltype==toupper(ct) & cytokine==toupper(stim)) %>%
  dplyr::filter(gene %in% deg.up$gene) 

length(unique(this_qtlgene$gene))
table(this_qtlgene$gene_class)

# down regulated genes
deg.down <- this_deg_metaconverted[this_deg_metaconverted$log2FoldChange < 0,]
length(unique(deg.down$gene))

this_qtlgene <- f_unambiguous_genes_all_cyto %>%
  dplyr::filter(celltype==toupper(ct) & cytokine==toupper(stim)) %>%
  dplyr::filter(gene %in% deg.down$gene) %>%
  dplyr::filter(anchorSNP_PBS_p<0.001 & anchorSNP_cytokine_p<0.001)

length(unique(this_qtlgene$gene))
table(this_qtlgene$gene_class)

df <- as.data.frame(table(this_qtlgene$gene_class)) %>%
  rename(gene_class = Var1, count = Freq) %>%
  arrange(desc(count))
ggplot(df, aes(x = reorder(gene_class, count), y = count)) +
  geom_col(fill = "#4C78A8") +
  geom_text(aes(label = count), hjust = -0.2) +
  coord_flip() +
  theme_classic() +
  expand_limits(y = max(df$count) * 1.1)

# ---- Data Summary ----
class_counts_all_cyto <- data.frame()
unambiguous_genes_all_cyto <- data.frame()
f_unambiguous_genes_all_cyto <- data.frame()

for (i in 1:3){
  this_cytokine=cytokines[i]
  
  # pick genes with 1 independent SNPs and store
  temp <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    dplyr::filter(anchorSNP_n_independent==1) %>%
    dplyr::mutate(nclass=1)
  
  # collapse gene with 2 independent SNPs to 1 SNP if they belong to the same class. Pick the one with strongest Z score.
  temp2 <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    dplyr::filter(anchorSNP_n_independent==2) 
  
  gene_lst <- unique(temp2$gene)
  collapsed_table <- data.frame()
  for (i in 1:length(gene_lst)) {
    this_slice <- temp2[temp2$gene==gene_lst[i],]
    if (this_slice[1,]$gene_class==this_slice[2,]$gene_class) { 
      this_res <- this_slice %>% arrange(desc(abs(anchorSNP_betaComparison_z))) %>% dplyr::slice(1)
      this_res$nclass <- 1
    } else {
      this_res <- this_slice
      this_res$nclass <- 2
    }
    collapsed_table <- rbind(collapsed_table, this_res)
  }
  
  #
  f_unambiguous_genes_all_cyto = rbind(f_unambiguous_genes_all_cyto, rbind(temp, collapsed_table))
  
  # Genes with 1 SNP — one row each, class is unambiguous
  single_snp_genes <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    dplyr::filter(anchorSNP_n_independent == 1) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%           # already 1 row per gene
    dplyr::select(gene, gene_class)
  
  # Genes with 2 SNPs, same class — already collapsed to 1 representative row
  two_snp_1class_genes <- collapsed_table %>%
    dplyr::filter(nclass == 1) %>%
    dplyr::select(gene, gene_class)
  
  # Combine and count unique genes per class
  unambiguous_genes <- bind_rows(single_snp_genes, two_snp_1class_genes) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>% # shouldn't overlap, but defensive
    dplyr::mutate(cytokine=this_cytokine)
  
  unambiguous_genes_all_cyto <- rbind(unambiguous_genes_all_cyto, unambiguous_genes)
  
  n_constitutive_stable <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_stable"),"gene"]))
  n_constitutive_amplified <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_amplified"),"gene"]))
  n_constitutive_attenuated <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_attenuated"),"gene"]))
  n_constitutive_switched <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_switched"),"gene"]))
  n_emergent_eQTL <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="emergent"),"gene"]))
  n_vanishing_eQTL <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="vanishing"),"gene"]))
  
  # Order bars by count descending
  class_counts <- data.table(gene_class = c("constitutive_stable", "constitutive_amplified",
                                            "constitutive_attenuated", "constitutive_switched",
                                            "emergent_eQTL", "vanishing_eQTL"),
                             n_genes = c(n_constitutive_stable, n_constitutive_amplified,
                                         n_constitutive_attenuated, n_constitutive_switched,
                                         n_emergent_eQTL, n_vanishing_eQTL))
  
  class_counts$gene_class <- factor(class_counts$gene_class,
                                    levels=c("constitutive_stable", "constitutive_amplified",
                                             "constitutive_attenuated", "constitutive_switched",
                                             "emergent_eQTL", "vanishing_eQTL"))
  
  class_counts$cytokine <- this_cytokine
  
  
  # rbind to main table
  class_counts_all_cyto <- rbind(class_counts_all_cyto, class_counts)
}

fwrite(f_unambiguous_genes_all_cyto, file= paste0(dir,"/reclassified/f_unambiguous_genes_all_cytokines.txt"), quote=F, sep="\t")

# ---- Join DE gene and QTL gene tables ----
unambiguous_genes_all_cyto_DEGs <- right_join(deg_wide, f_unambiguous_genes_all_cyto)
# commonly induced IFNG IFNB DE genes, but only have QTL in one condition
temp <- unambiguous_genes_all_cyto_DEGs %>%
  dplyr::filter(is_deg.krt.ifnb==T & is_deg.krt.ifng==T & is_deg.krt.tnf==F) %>%
  dplyr::select(c(gene, cytokine, gene_class, anchorSNP, priority_score)) %>%
  distinct() %>%
  arrange(gene, desc(priority_score))


# ---- Explore Fieller's method of classifying constitutive class ----
# Fieller's theorem CI for ratio beta_cyto / beta_PBS
fieller_ci <- function(b1, b0, se1, se0, t = 1.96) {
  g <- (t * se0 / b0)^2
  if (is.na(g) || g >= 1) return(c(NA_real_, NA_real_))  # denominator too uncertain
  
  disc <- b1^2 * se0^2 + b0^2 * se1^2 - t^2 * se0^2 * se1^2
  if (disc < 0) return(c(NA_real_, NA_real_))
  
  denom  <- b0^2 - t^2 * se0^2
  ci_lo  <- (b1 * b0 - t * sqrt(disc)) / denom
  ci_hi  <- (b1 * b0 + t * sqrt(disc)) / denom
  sort(c(ci_lo, ci_hi))
}

df2 <- f_unambiguous_genes_all_cyto %>%
  mutate(
    beta_PBS  = anchorSNP_PBS_beta,
    se_PBS    = anchorSNP_PBS_SE,
    beta_cyto = anchorSNP_cytokine_beta,
    se_cyto   = anchorSNP_cytokine_SE,
    ratio_hat = beta_cyto / beta_PBS
  ) %>%
  rowwise() %>%
  mutate(
    fieller = list(fieller_ci(beta_cyto, beta_PBS, se_cyto, se_PBS)),
    ratio_lo = fieller[[1]],
    ratio_hi = fieller[[2]]
  ) %>%
  ungroup() %>%
  select(-fieller) %>%
  mutate(
    gene_class_new = case_when(
      gene_class %in% c("emergent", "vanishing") ~ gene_class,
      
      # direction flip — must check before any magnitude-based cutoffs
      ratio_hat < 0                              ~ "constitutive_switched",
      
      # denominator too uncertain for reliable ratio estimate
      is.na(ratio_lo)                            ~ "constitutive_stable",
      
      # CI excludes 1 on the low side: significantly attenuated
      ratio_hi < 1                               ~ "constitutive_attenuated",
      
      # CI excludes 1 on the high side: significantly amplified
      ratio_lo > 1                               ~ "constitutive_amplified",
      
      # CI contains 1: no statistically supported change
      TRUE                                       ~ "constitutive_stable"
    )
  )

for (i in 1:3) {
  this_cytokine = cytokines[i]
  n_constitutive_attenuated <- nrow(df2[which(df2$cytokine==this_cytokine & df2$gene_class_new=="constitutive_attenuated"),])
  n_constitutive_switched <- nrow(df2[which(df2$cytokine==this_cytokine & df2$gene_class_new=="constitutive_switched"),])
  n_constitutive_amplified <- nrow(df2[which(df2$cytokine==this_cytokine & df2$gene_class_new=="constitutive_amplified"),])
  n_constitutive_stable <- nrow(df2[which(df2$cytokine==this_cytokine & df2$gene_class_new=="constitutive_stable"),])
  n_emergent <- nrow(df2[which(df2$cytokine==this_cytokine & df2$gene_class_new=="emergent"),])
  n_vanishing <- nrow(df2[which(df2$cytokine==this_cytokine & df2$gene_class_new=="vanishing"),])
  
  message(sprintf(
    "ct=%s | cytokine=%s | n_constitutive_attenuated=%i | n_constitutive_switched = %i | n_constitutive_amplified = %i | n_constitutive_stable = %i, n_emergent = %i, n_vanishing = %i",
    ct, cytokines[i], n_constitutive_attenuated, n_constitutive_switched, n_constitutive_amplified, n_constitutive_stable, n_emergent, n_vanishing))
}


df2 %>%
  count(cytokine, gene_class_new) %>%
  ggplot(aes(x = gene_class_new, y = n, fill = gene_class_new)) +
  geom_col() +
  facet_wrap(~ cytokine, scales = "free_y") +
  labs(
    x = "Gene class",
    y = "Number of genes",
    title = "Reclassified QTL gene classes by cytokine"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

top_examples <- data.frame()
for (i in 1:3){
  this_cytokine=cytokines[i]
  
  temp.top <- df2 %>%
    filter(cytokine == this_cytokine) %>%
    group_by(gene_class_new) %>%
    filter(!is.na(anchorSNP_betaComparison_z)) %>%
    arrange(anchorSNP_betaComparison_z, .by_group = TRUE) %>%
    slice_head(n = 3) 
  
  temp.bottom <- df2 %>%
    filter(cytokine == this_cytokine) %>%
    group_by(gene_class_new) %>%
    filter(!is.na(anchorSNP_betaComparison_z)) %>%
    arrange(desc(anchorSNP_betaComparison_z), .by_group = TRUE) %>%
    slice_head(n = 3) 
  
  temp <- rbind(temp.top, temp.bottom) %>% distinct()
  top_examples <- rbind(top_examples, temp)
}
top_examples_trim <- top_examples %>% dplyr::select(gene, anchorSNP) %>% distinct()

fwrite(top_examples_trim, file=paste0(dir,"/reclassified/top_examples.txt"), quote=F, sep="\t")

# ---- Explore ranking metric of assinging gene-prioritizing scores ----
# abs(beta_cyto) * abs_delta * -log10(p_cytokine)



# ---- Function Helpers ----
# Colour palette — one colour per class
class_colors <- c(
  "constitutive_stable"     = "#4393C3",
  "constitutive_amplified"  = "#2166AC",
  "constitutive_attenuated" = "#92C5DE",
  "constitutive_switched"   = "#762A83",
  "emergent_eQTL"           = "#D6604D",
  "vanishing_eQTL"          = "#F4A582"
)

# ---- Q1: How many genes with a valid class have 1 vs 2 independent SNPs? ----
for (i in 1:3) {
  temp <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==cytokines[i]) %>%
    dplyr::select(c(gene,anchorSNP_n_independent))
  
  n_1snp <- length(unique(temp[which(temp$anchorSNP_n_independent==1),]$gene))
  n_2snp <- length(unique(temp[which(temp$anchorSNP_n_independent==2),]$gene))
  
  message(sprintf(
    "ct=%s | cytokine=%s | nGenes_1SNP = %s | nGenes_2SNPs = %s",
    ct, cytokines[i], n_1snp, n_2snp
  ))
}

# ---- Q2: For genes with 2 independent SNPs: 1 unique class vs 2 unique classes? ----
for (i in 1:3) {
  this_cytokine = cytokines[i]
  temp2 <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    dplyr::filter(anchorSNP_n_independent==2) 
  
  gene_lst <- unique(temp2$gene)
  collapsed_table <- data.frame()
  for (i in 1:length(gene_lst)) {
    this_slice <- temp2[temp2$gene==gene_lst[i],]
    if (this_slice[1,]$gene_class==this_slice[2,]$gene_class) { 
      this_res <- this_slice %>% arrange(desc(abs(anchorSNP_betaComparison_z))) %>% dplyr::slice(1)
      this_res$nclass <- 1
    } else {
      this_res <- this_slice
      this_res$nclass <- 2
    }
    collapsed_table <- rbind(collapsed_table, this_res)
  }
  
  n_1class <- length(unique(collapsed_table[which(collapsed_table$nclass==1),]$gene))
  n_2class <- length(unique(collapsed_table[which(collapsed_table$nclass==2),]$gene))
  
  message(sprintf(
    "ct=%s | cytokine=%s | nGenes_with_2SNPs_1class = %s | nGenes_with_2SNPs_2class = %s",
    ct, this_cytokine, n_1class, n_2class
  ))
}

# ---- Q3: Bar chart of unique genes per class ----
# Universe: genes with a valid class AND a single unambiguous class label.
# This includes:
#   (a) genes with 1 independent SNP                    (already 1 class by definition)
#   (b) genes with 2 independent SNPs but same class    (collapsed above to 1 row)
# Genes with 2 independent SNPs AND 2 different classes are EXCLUDED (ambiguous).

for (i in 1:3){
  this_cytokine=cytokines[i]
  
  temp <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    dplyr::select(c(gene,anchorSNP_n_independent))
  
  n_1snp <- length(unique(temp[which(temp$anchorSNP_n_independent==1),]$gene))
  n_2snp <- length(unique(temp[which(temp$anchorSNP_n_independent==2),]$gene))
  
  temp2 <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    dplyr::filter(anchorSNP_n_independent==2) 
  
  gene_lst <- unique(temp2$gene)
  collapsed_table <- data.frame()
  for (i in 1:length(gene_lst)) {
    this_slice <- temp2[temp2$gene==gene_lst[i],]
    if (this_slice[1,]$gene_class==this_slice[2,]$gene_class) { 
      this_res <- this_slice %>% arrange(desc(abs(anchorSNP_betaComparison_z))) %>% dplyr::slice(1)
      this_res$nclass <- 1
    } else {
      this_res <- this_slice
      this_res$nclass <- 2
    }
    collapsed_table <- rbind(collapsed_table, this_res)
  }
  
  n_1class <- length(unique(collapsed_table[which(collapsed_table$nclass==1),]$gene))
  n_2class <- length(unique(collapsed_table[which(collapsed_table$nclass==2),]$gene))

  # For unambiguous genes (genes with one unique genes per class), count n in each class
  unambiguous_genes <- unambiguous_genes_all_cyto %>% dplyr::filter(cytokine==this_cytokine)
  
  n_constitutive_stable <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_stable"),"gene"]))
  n_constitutive_amplified <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_amplified"),"gene"]))
  n_constitutive_attenuated <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_attenuated"),"gene"]))
  n_constitutive_switched <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_switched"),"gene"]))
  n_emergent_eQTL <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="emergent"),"gene"]))
  n_vanishing_eQTL <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="vanishing"),"gene"]))
  
  message(sprintf("Total unambiguous genes for bar chart: %d", nrow(unambiguous_genes)))
  
  # Order bars by count descending
  class_counts <- data.table(gene_class = c("constitutive_stable", "constitutive_amplified",
                                            "constitutive_attenuated", "constitutive_switched",
                                            "emergent_eQTL", "vanishing_eQTL"),
                             n_genes = c(n_constitutive_stable, n_constitutive_amplified,
                                         n_constitutive_attenuated, n_constitutive_switched,
                                         n_emergent_eQTL, n_vanishing_eQTL))
  
  class_counts$gene_class <- factor(class_counts$gene_class,
                                    levels=c("constitutive_stable", "constitutive_amplified",
                                             "constitutive_attenuated", "constitutive_switched",
                                             "emergent_eQTL", "vanishing_eQTL"))
  
  p <- ggplot(class_counts, aes(x = gene_class, y = n_genes, fill = gene_class)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    geom_text(aes(label = n_genes), vjust = -0.4,
              size = 3.5, fontface = "bold", color = "grey25") +
    scale_fill_manual(values = class_colors, na.value = "grey70") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = sprintf("%s %s - QTL gene class distribution", ct, this_cytokine),
      subtitle = sprintf(
        "Unambiguous genes only: 1-SNP genes (n=%d) + 2-SNP same-class genes (n=%d) | excludes 2-SNP 2-class genes (n=%d)",
        n_1snp, n_1class, n_2class),
      x = "Gene class",
      y = "Number of unique genes"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position    = "none",
      axis.text.x        = element_text(angle = 35, hjust = 1, size = 11),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 8.5, color = "grey45"),
      panel.grid.major.y = element_line(color = "grey90")
    )
  
  out_plot <- paste0(dir, "/reclassified/", ct, "_", this_cytokine, "_class_barplot.pdf")
  ggsave(out_plot, p, width = 8, height = 5)
  message(sprintf("Saved: %s", out_plot))
}
plot_dir <- paste0(dir,"/reclassified/")
out_pdf  <- file.path(plot_dir, paste0(ct,"_class_barplots.pdf"))

pdfs <- list.files(plot_dir, pattern = "\\class_barplot.pdf$", full.names = TRUE)
pdfs <- sort(pdfs)   # IMPORTANT: controls page order
qpdf::pdf_combine(pdfs, output = out_pdf)

file.remove(list.files(plot_dir, pattern = "_class_barplot\\.pdf$", full.names = TRUE))

# ---- Q4: Stacked Bar chart of unique genes per class ----
# One grouped bar per cytokine (IFNB, IFNG, TNF), stacked by class, 
# show whether IFNγ drives more emergent eQTLs than TNFα, whether IFNβ is dominated by vanishing, etc.
# Requires running the unambiguous_genes logic for all 3 cytokines first,
# then rbinding with a cytokine column before plotting.
class_counts_all_cyto <- data.frame()
unambiguous_genes_all_cyto <- data.frame()
for (i in 1:3){
  this_cytokine=cytokines[i]
  
  # Genes with 1 SNP — one row each, class is unambiguous
  single_snp_genes <- f %>%
    dplyr::filter(gene_class!="") %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    dplyr::filter(anchorSNP_n_independent == 1) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%           # already 1 row per gene
    dplyr::select(gene, gene_class)
  
  # Genes with 2 SNPs, same class — already collapsed to 1 representative row
  two_snp_1class_genes <- collapsed_table %>%
    dplyr::filter(nclass == 1) %>%
    dplyr::select(gene, gene_class)
  
  # Combine and count unique genes per class
  unambiguous_genes <- bind_rows(single_snp_genes, two_snp_1class_genes) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>% # shouldn't overlap, but defensive
    dplyr::mutate(cytokine=this_cytokine)
  
  unambiguous_genes_all_cyto <- rbind(unambiguous_genes_all_cyto, unambiguous_genes)
  
  n_constitutive_stable <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_stable"),"gene"]))
  n_constitutive_amplified <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_amplified"),"gene"]))
  n_constitutive_attenuated <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_attenuated"),"gene"]))
  n_constitutive_switched <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="constitutive_switched"),"gene"]))
  n_emergent_eQTL <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="emergent"),"gene"]))
  n_vanishing_eQTL <- nrow(unique(unambiguous_genes[which(unambiguous_genes$gene_class=="vanishing"),"gene"]))
  
  # Order bars by count descending
  class_counts <- data.table(gene_class = c("constitutive_stable", "constitutive_amplified",
                                            "constitutive_attenuated", "constitutive_switched",
                                            "emergent_eQTL", "vanishing_eQTL"),
                             n_genes = c(n_constitutive_stable, n_constitutive_amplified,
                                         n_constitutive_attenuated, n_constitutive_switched,
                                         n_emergent_eQTL, n_vanishing_eQTL))
  
  class_counts$gene_class <- factor(class_counts$gene_class,
                                    levels=c("constitutive_stable", "constitutive_amplified",
                                             "constitutive_attenuated", "constitutive_switched",
                                             "emergent_eQTL", "vanishing_eQTL"))
  
  class_counts$cytokine <- this_cytokine
  
  
  # rbind to main table
  class_counts_all_cyto <- rbind(class_counts_all_cyto, class_counts)
}

ggplot(class_counts_all_cyto, aes(x = cytokine, y = n_genes, fill = gene_class)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = class_colors)

# ---- Q5: In Constitutive class, does current grouping make sense? ----
f_unambiguous_genes_all_cyto %>%
  ggplot(
    aes(
      x = anchorSNP_betaComparison_ratio,
      y = anchorSNP_betaComparison_z,
      color = gene_class
    )
  ) +
  geom_point(alpha = 0.7, size = 1.8) +
  facet_wrap(~ gene_class, scales = "free") +
  theme_bw() +
  labs(
    x = "Anchor SNP beta comparison ratio",
    y = "Anchor SNP beta comparison Z",
    title = "Anchor SNP beta comparison by gene class"
  )

# ---- Q6: Canonical reQTL enrichment by class ----
for (i in 1:3){
  this_cytokine=cytokines[i]
  
  f_unambiguous_genes_all_cyto %>%
    dplyr::filter(cytokine==this_cytokine) %>%
    count(gene_class, has_cytokine_reQTL) %>%
    group_by(gene_class) %>%
    mutate(pct = n / sum(n))
}





