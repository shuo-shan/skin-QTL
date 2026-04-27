library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)

ct <- "FRB"
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/", ct)
traits_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc/traits.txt"
traits <- fread(traits_file, header = FALSE)$V1

# function ----
summarize_one_trait <- function(trait, ct, condition, QTLtype, dir, eGene_list) {
  coloc_file <- paste0(dir, "/coloc/summary/coloc_", trait, "_", condition, "_", QTLtype, ".txt")
  
  if (!file.exists(coloc_file)) {
    return(data.frame(
      trait = trait,
      celltype = ct,
      condition = condition,
      QTLtype = QTLtype,
      n_eGenes = NA_integer_,
      n_H4 = NA_integer_,
      n_H4_GWAS = NA_integer_,
      genes_H4 = NA_character_,
      genes_H4_GWAS = NA_character_
    ))
  }
  
  coloc_res <- fread(coloc_file) %>%
    filter(gene %in% eGene_list) %>%
    arrange(desc(PP.H4))
  
  if (nrow(coloc_res) == 0) {
    return(data.frame(
      trait = trait,
      celltype = ct,
      condition = condition,
      QTLtype = QTLtype,
      n_eGenes = 0,
      n_H4 = 0,
      n_H4_GWAS = 0,
      genes_H4 = "",
      genes_H4_GWAS = ""
    ))
  }
  
  coloc_res <- coloc_res %>%
    mutate(
      pass_H4 = PP.H4 > 0.7,
      pass_GWAS = GWAS_best_p < 5e-8,
      pass_both = pass_H4 & pass_GWAS
    )
  
  genes_H4 <- coloc_res %>%
    filter(pass_H4) %>%
    pull(gene) %>%
    unique()
  
  genes_H4_GWAS <- coloc_res %>%
    filter(pass_both) %>%
    pull(gene) %>%
    unique()
  
  gene_leadQTL_H4 <- coloc_res %>%
    filter(pass_H4) %>%
    mutate(temp = paste0(gene,":",QTL_best_rsid)) %>%
    pull(temp) %>%
    unique()
  
  data.frame(
    trait = trait,
    celltype = ct,
    condition = condition,
    QTLtype = QTLtype,
    n_eGenes = nrow(coloc_res),
    n_H4 = sum(coloc_res$pass_H4, na.rm = TRUE),
    n_H4_GWAS = sum(coloc_res$pass_both, na.rm = TRUE),
    genes_H4 = paste(genes_H4, collapse = ","),
    genes_H4_GWAS = paste(genes_H4_GWAS, collapse = ","),
    gene_leadQTL_H4 = paste(gene_leadQTL_H4, collapse = ",")
  )
}

summarize_all_traits_all_QTLs <- data.frame()
# run function across all QTL ----
for (condition in c("PBS","IFNB","IFNG","TNF")) {
  for (QTLtype in c("eQTL", "reQTL")) {
    if (!(condition=="PBS" && QTLtype=="reQTL")) {
      
      message(paste0(condition, QTLtype))
      eGene_file <- paste0(dir, "/eigenMT/results/", ct, "_", condition, "_", QTLtype, "_gene_fdr05_genelist.txt")
      eGene_list <- tryCatch(
        {
          fread(eGene_file, header = FALSE)$V1
        },
        error = function(e) {
          character(0)
        }
      )
      summary_all_traits <- map_dfr(
        traits,
        summarize_one_trait,
        ct = ct,
        condition = condition,
        QTLtype = QTLtype,
        dir = dir,
        eGene_list = eGene_list
      )
      
      summary_all_traits <- summary_all_traits %>%
        dplyr::filter(trait!="height") %>%
        arrange(desc(n_H4_GWAS), desc(n_H4), trait)
      
      summarize_all_traits_all_QTLs <- dplyr::bind_rows(summarize_all_traits_all_QTLs, summary_all_traits)
      fwrite(
        summary_all_traits,
        file = paste0(dir, "/coloc/summary/", ct, "_", condition, "_", QTLtype, "_coloc_trait_summary.tsv"),
        sep = "\t"
      )
      
    }
  }
}

fwrite(
  summarize_all_traits_all_QTLs,
  file = paste0(dir, "/coloc/summary/master_coloc_trait_summary.tsv"),
  sep = "\t"
)


# heatmap ----
df <- summarize_all_traits_all_QTLs %>%
  dplyr::mutate(
    cond_QTL = paste(condition, QTLtype, sep = "_"),
    label = paste0(n_H4_GWAS, ", ", n_H4)
  )

# set desired column order
df$cond_QTL <- factor(
  df$cond_QTL,
  levels = c(
    "PBS_eQTL",
    "IFNB_eQTL",
    "IFNG_eQTL",
    "TNF_eQTL",
    "IFNB_reQTL",
    "IFNG_reQTL",
    "TNF_reQTL"
  )
)

# order traits by total n_H4_GWAS across all columns
trait_order <- df %>%
  group_by(trait) %>%
  summarise(total_H4_GWAS = sum(n_H4_GWAS, na.rm = TRUE)) %>%
  arrange(desc(total_H4_GWAS)) %>%
  pull(trait)

df$trait <- factor(df$trait, levels = trait_order)

# add eGene number to column label
col_labels_df <- summarize_all_traits_all_QTLs %>%
  dplyr::distinct(condition, QTLtype, n_eGenes) %>%
  dplyr::mutate(
    cond_QTL = paste(condition, QTLtype, sep = "_"),
    label = paste0(condition, "_", QTLtype, "\n(", n_eGenes, " eGenes)")
  )
label_map <- setNames(col_labels_df$label, col_labels_df$cond_QTL)

# plot
p <- ggplot(df, aes(x = cond_QTL, y = trait, fill = n_H4_GWAS)) +
  geom_tile(color = "grey85", linewidth = 0.3) +   # ← faint grid lines
  geom_text(aes(label = label), size = 3) +
  scale_x_discrete(labels = label_map) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    x = NULL,
    y = NULL,
    fill = "n_H4_GWAS",
    title = paste0(ct, ": GWAS-colocalized QTL genes across traits"),
    subtitle = "Each cell: # eGenes with PP.H4 > 0.7 & GWAS p < 5e-8\nright = # eGenes with PP.H4 > 0.7"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),  # keep background clean
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    plot.title = element_text(face = "bold")
  )

print(p)

pdf(file=paste0(dir,"/coloc/summary/coloc_summary_heatmap.pdf"), width=10, height=8)
print(p)
dev.off()











# explore annotation table ----
dir="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/coloc/summary/best/"
gene="IRF3"
snp="rs3768018"

df=fread(paste0(dir,gene,"/annotation_big.txt"))

# which SNP has TF motif overlap?
this_df <- df %>% 
  dplyr::filter(`n_TFmotif_1E-4`>0) %>%
  dplyr::select(QTL,
                `n_TFmotif_1E-8`, `TFmotif_1E-8`,
                peak_overlapping_TF,
                AlleleSpecificMark,
                cRE_name,  cRE_dynamic_FRB,
                Func.refGene, Gene.refGene, dist_to_tss,
                `n_TFmotif_1E-6`, `TFmotif_1E-6`,
                `n_TFmotif_1E-4_ENCODETF`, `TFmotif_1E-4_ENCODETF`,
                `n_TFmotif_1E-4`, `TFmotif_1E-4`) %>%
  arrange(desc(`n_TFmotif_1E-8`), desc(`n_TFmotif_1E-6`), desc(`n_TFmotif_1E-4`))
View(this_df)

