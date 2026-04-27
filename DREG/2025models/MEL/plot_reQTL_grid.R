#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

# ---------------------- Argument Parsing ----------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[1])
}

pairs_file <- get_arg("--pairs")
GENO_FILE  <- get_arg("--geno")
STATS_FILE <- get_arg("--stats")
CPM_FILE   <- get_arg("--cpm")
META_FILE  <- get_arg("--meta")
OUT_PDF    <- get_arg("--out")

# ---------------------- Optional Hardcoded Paths ----------------------
# Uncomment to bypass command line input:
STATS_FILE <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/results/reQTL_IFNG_pvalE-09.txt"
GENO_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots/genotype_reQTL_IFNG_pvalE-09.tsv"
pairs_file <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots/pairs_reQTL_IFNG_pvalE-09.txt"
CPM_FILE   <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
META_FILE  <- "/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
OUT_PDF    <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/plots/reQTL_by_condition.pdf"

stopifnot(!is.null(pairs_file), !is.null(GENO_FILE), !is.null(STATS_FILE),
          !is.null(CPM_FILE), !is.null(META_FILE))

# ---------------------- Load SNP-gene Pairs ----------------------
message("[1] Loading SNP-gene pairs...")
pairs <- fread(pairs_file, header = FALSE)
colnames(pairs) <- c("snp","gene")

# ---------------------- Load Expression ----------------------
message("[2] Loading CPM expression...")
cpm <- fread(CPM_FILE) %>%
  select(-c("gene","name")) %>%
  distinct(final_gene, .keep_all = TRUE) %>%
  column_to_rownames("final_gene") %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id %in% pairs$gene)

expr_long <- as.data.table(cpm) %>%
  melt(id.vars = "gene_id", variable.name = "sample", value.name = "CPM")
rm(cpm); invisible(gc())

# ---------------------- Metadata ----------------------
message("[3] Loading metadata...")
meta <- fread(META_FILE)
expr_long <- expr_long %>%
  inner_join(meta, by = "sample") %>%
  mutate(condition = as.character(condition),
         celltype  = as.character(celltype))

# ---------------------- Genotype Dosage ----------------------
message("[4] Loading genotype dosage...")
geno_dt <- fread(GENO_FILE)
fixed_cols <- c("CHROM","POS","ID","REF","ALT")

old_cols <- setdiff(colnames(geno_dt), fixed_cols)
new_cols <- gsub("skineQTL-", "", old_cols)
new_cols <- gsub("^F0", "F", new_cols)
setnames(geno_dt, old_cols, new_cols)

geno_samples <- setdiff(colnames(geno_dt), fixed_cols)
expr_donors <- unique(expr_long$donor)
keep <- intersect(geno_samples, expr_donors)
geno_dt <- geno_dt[, c(fixed_cols, keep), with=FALSE]
setkey(geno_dt, ID)

dosage_to_label <- function(d, ref, alt){
  factor(ifelse(d==0, paste0(ref,ref),
                ifelse(d==1, paste0(ref,alt),
                       ifelse(d==2, paste0(alt,alt), NA))),
         levels=c(paste0(ref,ref),paste0(ref,alt),paste0(alt,alt)))
}

get_genotype_df <- function(snp){
  row <- geno_dt[.(snp)]
  if (nrow(row) == 0) return(NULL)
  
  ref <- row$REF[1]
  alt <- row$ALT[1]
  
  # Extract numeric dosage vector correctly
  dos <- unlist(row[, ..keep], use.names = TRUE)
  dos <- as.numeric(dos)   # convert "0","1","2" to numeric
  names(dos) <- keep       # make sure donor names are preserved
  
  tibble(
    donor = names(dos),
    genotype = dosage_to_label(dos, ref, alt)
  )
}


# ---------------------- Stats ----------------------
message("[5] Loading reQTL stats (no header)...")
stats <- fread(STATS_FILE, header=FALSE)
colnames(stats) <- c("celltype","snp","gene",
                     "eQTL_beta_PBS","eQTL_beta_IFNG","eQTL_beta_IFNB","eQTL_beta_TNF",
                     "reQTL_dbeta_IFNG","reQTL_p_IFNG",
                     "reQTL_dbeta_IFNB","reQTL_p_IFNB",
                     "reQTL_dbeta_TNF","reQTL_p_TNF")

get_stats <- function(s, g, ct, stim){
  row <- stats %>% filter(snp == s, gene == g, celltype == ct)
  if(nrow(row) == 0) return(list(beta = NA, p = NA))
  list(
    beta = row[[paste0("reQTL_dbeta_", stim)]][1],
    p    = row[[paste0("reQTL_p_", stim)]][1]
  )
}


fmt_stat <- function(beta,p){
  if(is.na(beta)|is.na(p)) return("p = NA, β = NA")
  stars <- case_when(p<1e-4~"***", p<1e-3~"**", p<1e-2~"*", TRUE~"")
  paste0("p=", format(p, digits=1, scientific=TRUE),
         ", beta=", sprintf("%.1f",beta),", ",stars)
}

# ---------------------- Plot Function ----------------------
# ---------------------- Plot Helper (PBS vs condition panel) ----------------------
make_subpanel <- function(df, ct, stim, snp, gene){
  
  sub <- df %>% filter(celltype == ct, condition %in% c("PBS", stim))
  if(nrow(sub) == 0) return(NULL)
  
  # Count donors per genotype+condition
  counts <- sub %>%
    group_by(genotype, condition) %>%
    summarize(n = n_distinct(donor), .groups="drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
    mutate(strip_label = paste0(genotype, "\nPBS=", PBS, " ", stim, "=", !!sym(stim))) %>%
    select(genotype, strip_label)
  
  # Attach counts to genotype label
  sub <- sub %>%
    left_join(counts, by="genotype") %>%
    mutate(genotype = factor(strip_label, levels=unique(strip_label))) %>%
    select(-strip_label)
  
  # Ensure PBS is left
  sub$condition <- factor(sub$condition, levels=c("PBS", stim))
  
  # Stats
  st <- get_stats(snp, gene, ct, stim)
  subtitle <- paste0(ct, " — PBS vs ", stim, "  (", fmt_stat(st$beta, st$p), ")")
  
  ggplot(sub, aes(x=condition, y=CPM, group=donor)) +
    geom_line(alpha=0.25, color="grey70") +
    geom_point(aes(color=condition), size=2.2, alpha=0.9) +
    scale_color_manual(values=c(
      "PBS"="gray35",
      "IFNG"="#0072B2",
      "IFNB"="#E69F00",
      "TNF"="#D55E00"
    )) +
    facet_wrap(~genotype, scales="fixed", nrow=1) +
    labs(title = subtitle, y="CPM", x="") +
    theme_bw() +
    theme(
      plot.title = element_text(size=3),   # ↓ smaller subpanel titles
      strip.text = element_text(size=8, lineheight=0.8), # ↓ smaller genotype + count display
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=8),
      axis.text.y = element_text(size=8),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
}


# ---------------------- Main Plot for a SNP : gene ----------------------
make_plot <- function(snp, gene){
  
  df <- expr_long %>% filter(gene_id == gene)
  geno <- get_genotype_df(snp)
  df <- df %>% inner_join(geno, by="donor")
  
  ct_order <- c("FRB","KRT","MEL")
  stim_order <- c("IFNG","IFNB","TNF")
  
  panels <- list()
  for(ct in ct_order){
    for(stim in stim_order){
      p <- make_subpanel(df, ct, stim, snp, gene)
      if(!is.null(p)) panels[[paste(ct, stim)]] <- p
    }
  }
  
  combined <- wrap_plots(panels, ncol=3)
  
  # Apply global title WITHOUT overriding sub-panel title sizes
  combined <- combined + plot_annotation(
    title = paste0(snp, " : ", gene)
  )
  
  combined & theme(
    plot.title = element_text(size=14, face="bold", hjust=0.5)  # affects ONLY main title now
  )
}

# ---------------------- Output ----------------------
# message("[6] Writing PDF: ", OUT_PDF)
# pdf(OUT_PDF, width=16, height=8)
# #for(i in seq_len(nrow(pairs))){
# for(i in seq_len(10)){
#   print(make_plot(pairs$snp[i], pairs$gene[i]))
# }
# dev.off()
# message("Done")

# ----------------------- Output best SNP per gene ----
lead <- stats %>%
  tidyr::pivot_longer(
    cols = c(starts_with("reQTL_p_"), starts_with("reQTL_dbeta_")),
    names_to = c(".value","stim"),
    names_pattern = "(reQTL_.*?)_(IFNG|IFNB|TNF)"
  ) %>%
  rename(p = reQTL_p, dbeta = reQTL_dbeta) %>%
  group_by(celltype, gene, stim) %>%
  slice_min(p, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(c("snp","gene")) %>%
  distinct

message("[6] Writing PDF: ", OUT_PDF)
pdf(OUT_PDF, width=16, height=8)
#for(i in seq_len(nrow(lead))){
for(i in seq_len(20)){
  print(make_plot(lead$snp[i], lead$gene[i]))
  message(paste0(lead$snp[i],":",lead$gene[i]))
}
dev.off()
message("Done")
