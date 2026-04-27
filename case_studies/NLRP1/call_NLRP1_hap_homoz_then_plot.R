#!/usr/bin/env Rscript
# Usage: Rscript call_NLRP1_hap_homoz_then_plot.R

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

dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/case_studies/NLRP1"
tags_f <- "/pi/manuel.garber-umw/human/skin/eQTLs/case_studies/NLRP1/NLRP1.hap.tags.txt"
defs_f <- "/pi/manuel.garber-umw/human/skin/eQTLs/case_studies/NLRP1/NLRP1_haplotype_definitions.tsv"
outp   <- "NLRP1"
CPM_FILE   <- paste0(dir, "/cpm.txt")
META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"


# ---------------------- NLRP1 Haplotype Calls ---------------------- ####
message("[1] Determine NLRP1 Haplotype from donors...")
tags <- read.table(tags_f, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
defs <- read.table(defs_f, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Merge by rsid
m <- merge(tags, defs, by.x="ID", by.y="rsid", all=FALSE)

donors <- setdiff(colnames(m), c("CHROM","POS","ID","REF","ALT","hap1","hap2A"))

# Identify the SNPs that differ between Hap1 and Hap2A
diff_idx <- which(m$hap1 != m$hap2A)
m_diff <- m[diff_idx, ]
K <- nrow(m_diff)
message(sprintf("[Hap dosage] Using %d SNPs that differ between Hap1 and Hap2A", K))

# Convert each donor’s 0/1/2 into “# Hap1 alleles” per SNP
# vector of TRUE/FALSE: at each diff SNP, is Hap1 allele the ALT?
hap1_is_alt <- (m_diff$hap1 == m_diff$ALT)

hap1_dosage <- sapply(donors, function(d){
  gt <- as.integer(m_diff[[d]])  # 0/1/2
  # if missing, make it NA (or treat as 0; your choice)
  if (any(is.na(gt))) return(NA_integer_)
  
  hap1_alleles_per_snp <- ifelse(hap1_is_alt, gt, 2L - gt)
  sum(hap1_alleles_per_snp)
})

# compute the “SNP-count dosage” (0..K)
# hap1_dosage == 0 → consistent with Hap2A/Hap2A alleles across diff SNPs
# hap1_dosage == 2K → consistent with Hap1/Hap1 across diff SNPs
hap1_hom_count <- sapply(donors, function(d){
  gt <- as.integer(m_diff[[d]])
  if (any(is.na(gt))) return(NA_integer_)
  
  hap1_alleles_per_snp <- ifelse(hap1_is_alt, gt, 2L - gt)
  sum(hap1_alleles_per_snp == 2L)
})

# Determine required homozygote code for hap1 and hap2A at each SNP
req_code <- function(hap_allele, ref, alt){
  if (hap_allele == ref) return(0L)
  if (hap_allele == alt) return(2L)
  return(NA_integer_)
}

# compute req codes on full m first
m$req_hap1  <- mapply(req_code, m$hap1,  m$REF, m$ALT)
m$req_hap2A <- mapply(req_code, m$hap2A, m$REF, m$ALT)

# then define m_diff
diff_idx <- which(m$hap1 != m$hap2A)
m_diff <- m[diff_idx, ]

hap1_homs <- donors[sapply(donors, function(d) is_hom_all(as.integer(m_diff[[d]]), m_diff$req_hap1))]
hap2A_homs <- donors[sapply(donors, function(d) is_hom_all(as.integer(m_diff[[d]]), m_diff$req_hap2A))]

# Function to test “matches all SNPs” for a donor
is_hom_all <- function(gt_vec, req_vec){
  # gt_vec are numeric 0/1/2; allow NA to fail
  all(!is.na(gt_vec) & gt_vec == req_vec)
}

writeLines(hap1_homs, paste0(dir,"/",outp, ".Hap1_Hap1.donors.txt"))
writeLines(hap2A_homs, paste0(dir,"/",outp, ".Hap2A_Hap2A.donors.txt"))

# Also output a per-donor call table
call <- data.frame(
  donor = donors,
  Hap1_Hap1 = donors %in% hap1_homs,
  Hap2A_Hap2A = donors %in% hap2A_homs,
  hap1_dosage_alleles = as.integer(hap1_dosage),     # 0..2K
  hap1_dosage_homSNPs = as.integer(hap1_hom_count),  # 0..K (optional)
  stringsAsFactors = FALSE
)

call$diplotype <- ifelse(call$Hap1_Hap1, "Hap1/Hap1",
                         ifelse(call$Hap2A_Hap2A, "Hap2A/Hap2A", "Other"))

call$gradientH1toH2A <- (10-call$hap1_dosage_alleles)/10

write.table(call, paste0(dir,"/",outp, ".diplotype_calls.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("Done.\n")
cat("Hap1/Hap1:", length(hap1_homs), "\n")
cat("Hap2A/Hap2A:", length(hap2A_homs), "\n")

# ---------------------- NLRP1 Gene Expression (CPM) ---------------------- ####
message("[2] Loading CPM expression...")
cpm <- fread(CPM_FILE) %>%
  select(-c("gene","name")) %>%
  distinct(final_gene, .keep_all = TRUE) %>%
  column_to_rownames("final_gene") %>%
  rownames_to_column("gene_id") 

expr_long <- as.data.table(cpm) %>%
  melt(id.vars = "gene_id", variable.name = "sample", value.name = "CPM")
rm(cpm); invisible(gc())

# ---------------------- Sample Metadata ---------------------- ####
message("[3] Loading metadata...")
meta <- fread(META_FILE)
expr_long <- expr_long %>%
  inner_join(meta, by = "sample") %>%
  mutate(condition = as.character(condition),
         celltype  = as.character(celltype))

# ---------------------- Merge Haplotype Dosage with Expression & Metadata ---------------------- ####
expr_long <- expr_long %>%
  inner_join(call, by = "donor")

# ---------------------- Plot Helper (PBS vs condition panel) ----------------------
make_subpanel <- function(df, ct, stim){
  
  sub <- df %>% filter(celltype == ct, condition %in% c("PBS", stim))
  if(nrow(sub) == 0) return(NULL)
  
  # Count donors per genotype + condition
  counts <- sub %>%
    group_by(diplotype, condition) %>%
    summarize(n = n_distinct(donor), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
    mutate(strip_label = paste0(diplotype, "\nPBS=", PBS, " ", stim, "=", !!sym(stim)))
  
  # join counts into sub
  sub <- sub %>%
    left_join(counts, by = c("diplotype"))
  
  # ---- enforce facet order by numeric dosage ----
  sub <- sub %>%
    arrange(gradientH1toH2A) %>%   # ensures 0 → 1 → 2 before converting to factor
    mutate(
      genotype_label = factor(strip_label, levels = unique(strip_label))
    )
  
  # ensure PBS appears left on x-axis
  sub$condition <- factor(sub$condition, levels = c("PBS", stim))
  
  # ---- plotting ----
  subtitle <- paste0(ct, " — PBS vs ", stim)
  
  ggplot(sub, aes(x = condition, y = CPM, group = donor)) +
    geom_line(alpha = 0.25, color = "grey70") +
    geom_point(aes(color = condition), size = 2.2, alpha = 0.9) +
    scale_color_manual(values = c(
      "PBS" = "gray35",
      "IFNG" = "#0072B2",
      "IFNB" = "#E69F00",
      "TNF" = "#D55E00"
    )) +
    facet_wrap(~ genotype_label, scales = "fixed", nrow = 1) +
    labs(title = subtitle, y = "CPM", x = "") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 3),
      strip.text = element_text(size = 8, lineheight = 0.8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.major.x = element_blank(),   # already removing vertical major grid
      panel.grid.major.y = element_blank(),   # ← removes horizontal major grid
      panel.grid.minor.y = element_blank()    # ← removes horizontal minor grid
    )
  
}

plot_cpm_by_haplotype_gradient <- function(df, ct){
  sub <- df %>% 
    filter(celltype == ct) 
  
  sub$condition <- ordered(sub$condition, levels = c("PBS","IFNG","IFNB","TNF"))
  
  subtitle <- paste0(ct)
  
  ggplot(sub, aes(x = gradientH1toH2A, y = CPM, color = diplotype)) +
    geom_point(aes(color = condition), size = 2.2, alpha = 0.9) +
    scale_color_manual(values = c(
      "PBS" = "gray35",
      "IFNG" = "#0072B2",
      "IFNB" = "#E69F00",
      "TNF" = "#D55E00"
    )) +
    facet_wrap(~ condition, scales = "fixed", nrow = 1) +
    labs(title = subtitle, y = "CPM", x = "") +
    lims(x=c(-0.1, 1.1)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 3),
      strip.text = element_text(size = 8, lineheight = 0.8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.major.x = element_blank(),   # already removing vertical major grid
      panel.grid.major.y = element_blank(),   # ← removes horizontal major grid
      panel.grid.minor.y = element_blank()    # ← removes horizontal minor grid
    )
}
# ---------------------- Main Plot for a SNP : gene ----------------------
make_plot <- function(){
  
  df <- expr_long 
  
  ct_order <- c("FRB","KRT","MEL")
  stim_order <- c("IFNG","IFNB","TNF")
  
  panels <- list()
  for(ct in ct_order){
    for(stim in stim_order){
      p <- make_subpanel(df, ct, stim)
      if(!is.null(p)) panels[[paste(ct, stim)]] <- p
    }
  }
  
  combined <- wrap_plots(panels, ncol=3)
  
  combined <- combined + plot_annotation(
    title = "NLRP1 mRNA level by haplotype"
  )
  combined & theme(
    plot.title = element_text(size=7, face="bold", hjust=0.5)  # affects ONLY main title now
  )
}

make_cpm_by_haplotype_gradient_plot <- function(){
  
  df <- expr_long 
  
  ct_order <- c("FRB","KRT","MEL")
  
  panels <- list()
  for(ct in ct_order){
      p <- plot_cpm_by_haplotype_gradient(df, ct)
      if(!is.null(p)) panels[[paste(ct)]] <- p
  }
  
  combined <- wrap_plots(panels, nrow=3)
  
  combined <- combined + plot_annotation(
    title = "NLRP1 mRNA level by haplotype gradient (1 to 2A)"
  )
  combined & theme(
    plot.title = element_text(size=7, face="bold", hjust=0.5)  # affects ONLY main title now
  )
}
# ---------------------- Output ----------------------
message("[5] Writing PDF: ", OUT_PDF)
OUT_PDF <- paste0(dir,"/response_plot.pdf")
pdf(OUT_PDF, width=16, height=12)
print(make_plot())
dev.off()


OUT_PDF <- paste0(dir,"/NLRP1_vs_gradient_plot.pdf")
pdf(OUT_PDF, width=16, height=12)
print(make_cpm_by_haplotype_gradient_plot())
dev.off()
message("Done")

