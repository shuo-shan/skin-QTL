library(tidyverse)
library(coloc)
library(ggplot2)
library(patchwork)
library(rtracklayer)

# ============================================================
# SETTINGS
# ============================================================
TICAM1_chr   <- 19
TICAM1_start <- 4315932   # gene center - 500kb
TICAM1_end   <- 5331712   # gene center + 500kb

gwas_dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc"
out_dir <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/data/vitiligo_vs_otherGWAS_TICAM1"

# all traits to compare against vitiligo
traits <- list(
  psoriasis                    = "psoriasis/standardized_psoriasis.chr19.tsv.gz",
  systemic_lupus               = "systemic_lupus_erythematosus/standardized_systemic_lupus_erythematosus.chr19.tsv.gz",
  cutaneous_lupus              = "cutaneous_lupus_erythematosus/standardized_cutaneous_lupus_erythematosus.chr19.tsv.gz",
  atopic_dermatitis            = "atopic_dermatitis/standardized_atopic_dermatitis.chr19.tsv.gz",
  alopecia_areata              = "alopecia_areata/standardized_alopecia_areata.chr19.tsv.gz",
  rheumatoid_arthritis         = "rheumatoid_arthritis/standardized_rheumatoid_arthritis.chr19.tsv.gz",
  sunburn                      = "sunburn/standardized_sunburn.chr19.tsv.gz",
  basal_cell_carcinoma = "basal_cell_carcinoma/standardized_basal_cell_carcinoma.chr19.tsv.gz",
  crohns_disease = "crohns_disease/standardized_crohns_disease.chr19.tsv.gz",
  Melanomas_of_skin_dx_or_hx = "Melanomas_of_skin_dx_or_hx/standardized_Melanomas_of_skin_dx_or_hx.chr19.tsv.gz",
  skin_pigmentation = "skin_pigmentation/standardized_skin_pigmentation.chr19.tsv.gz",
  squamous_cell_carcinoma = "squamous_cell_carcinoma/standardized_squamous_cell_carcinoma.chr19.tsv.gz"
)

# ============================================================
# HELPER: read + subset to TICAM1 locus
# ============================================================
read_locus <- function(f, chr, start, end) {
  read_tsv(f, show_col_types = FALSE) %>%
    filter(chr == !!chr, pos >= start, pos <= end) %>%
    filter(!is.na(beta), !is.na(varbeta), !is.na(MAF), MAF > 0) %>%
    mutate(snp_id = if_else(!is.na(rsid) & rsid != ".", rsid, variant_id)) %>%
    filter(!duplicated(snp_id))
}

read_locus_jin <- function(chr_num, start, end) {
  f <- file.path(out_dir, sprintf("jin2016_vitiligo_chr%s_hg38.tsv.gz", chr_num))
  read_tsv(f, show_col_types = FALSE) %>%
    filter(pos >= start, pos <= end) %>%
    filter(!is.na(beta), !is.na(varbeta), !is.na(maf), maf > 0) %>%
    mutate(snp_id = snp) %>%
    filter(!duplicated(snp_id))
}



# ============================================================
# HELPER: build coloc dataset list
# ============================================================
make_coloc_dataset <- function(df) {
  trait_type <- df$type[1]
  
  if (trait_type == "quant") {
    cat("    type = quant, using N only\n")
    return(list(
      snp     = df$snp_id,
      beta    = df$beta,
      varbeta = df$varbeta,
      MAF     = df$MAF,
      type    = "quant",
      N       = df$N[1]
    ))
  }
  
  # cc: compute s if missing
  s_val <- df$s[1]
  if (is.na(s_val)) {
    s_val <- df$num_cases[1] / (df$num_cases[1] + df$num_controls[1])
    cat(sprintf("    s computed from case/control counts: %.4f\n", s_val))
  }
  
  list(
    snp     = df$snp_id,
    beta    = df$beta,
    varbeta = df$varbeta,
    MAF     = df$MAF,
    type    = "cc",
    s       = s_val,
    N       = df$N[1]
  )
}

make_coloc_dataset_jin <- function(df) {
  list(
    snp     = df$snp_id,
    beta    = df$beta,
    varbeta = df$varbeta,
    MAF     = df$maf,
    type    = "cc",
    s       = jin_s,
    N       = jin_N
  )
}
# ============================================================
# HELPER: align two dfs on shared SNPs
# ============================================================
align_pair <- function(df1, df2) {
  shared <- intersect(df1$snp_id, df2$snp_id)
  cat(sprintf("    Shared SNPs: %d\n", length(shared)))
  list(
    d1 = df1 %>% filter(snp_id %in% shared) %>% arrange(snp_id),
    d2 = df2 %>% filter(snp_id %in% shared) %>% arrange(snp_id)
  )
}

# ============================================================
# LOAD vitiligo once
# ============================================================
# load jin2016 locus
cat("Loading jin2016 TICAM1 locus...\n")
# s for jin2016
jin_s <- 2853 / (2853 + 37405)
jin_N <- 2853 + 37405

jin_locus <- read_locus_jin(TICAM1_chr, TICAM1_start, TICAM1_end)
vit <- jin_locus
cat(sprintf("  vitiligo jin2016 SNPs in locus: %d\n", nrow(jin_locus)))

# ============================================================
# RUN COLOC for each trait
# ============================================================
results <- list()

for (trait_name in names(traits)) {
  cat(sprintf("\n=== Vitiligo × %s at TICAM1 ===\n", trait_name))
  
  f <- file.path(gwas_dir, traits[[trait_name]])
  
  # gracefully skip if file missing
  if (!file.exists(f)) {
    cat(sprintf("  WARNING: file not found, skipping: %s\n", f))
    next
  }
  
  trait_df <- read_locus(f, TICAM1_chr, TICAM1_start, TICAM1_end)
  cat(sprintf("  %s SNPs in locus: %d\n", trait_name, nrow(trait_df)))
  
  if (nrow(trait_df) < 10) {
    cat("  WARNING: too few SNPs, skipping\n")
    next
  }
  
  pair <- align_pair(vit, trait_df)
  
  if (length(pair$d1$snp) < 10) {
    cat("  WARNING: too few shared SNPs, skipping\n")
    next
  }
  
  res <- coloc.abf(
    dataset1 = make_coloc_dataset(pair$d1),
    dataset2 = make_coloc_dataset(pair$d2)
  )
  
  print(res$summary)
  results[[trait_name]] <- res
}

# ============================================================
# SUMMARY TABLE
# ============================================================
summary_tbl <- map_dfr(names(results), function(trait_name) {
  s <- results[[trait_name]]$summary
  tibble(
    comparison = paste0("vitiligo_vs_", trait_name),
    locus      = "TICAM1",
    nsnps      = s["nsnps"],
    PP.H0      = s["PP.H0.abf"],
    PP.H1      = s["PP.H1.abf"],
    PP.H2      = s["PP.H2.abf"],
    PP.H3      = s["PP.H3.abf"],
    PP.H4      = s["PP.H4.abf"]
  )
}) %>%
  arrange(desc(PP.H4))

cat("\n\n========== FINAL SUMMARY (sorted by PP.H4) ==========\n")
print(summary_tbl, n = Inf)

# save
out_file <- file.path(out_dir, "TICAM1_coloc_vitiligo_vs_all_traits.tsv")
write_tsv(summary_tbl, out_file)
cat(sprintf("\nSaved to: %s\n", out_file))

# ============================================================
# GWAS × GWAS locus plot
# ============================================================

plot_gwas_gwas_locus <- function(df1, df2,
                                 trait1_name, trait2_name,
                                 coloc_res,
                                 gene_name = "TICAM1",
                                 gene_pos  = 4823822) {
  shared <- intersect(df1$snp_id, df2$snp_id)
  d1 <- df1 %>% filter(snp_id %in% shared) %>% arrange(snp_id)
  d2 <- df2 %>% filter(snp_id %in% shared) %>% arrange(snp_id)
  
  wide <- d1 %>%
    transmute(snp = snp_id, pos = pos, DIST = pos - gene_pos,
              P_trait1 = p, MAF = MAF) %>%
    left_join(d2 %>% transmute(snp = snp_id, P_trait2 = p), by = "snp") %>%
    mutate(
      P_trait1      = pmax(P_trait1, 1e-300),
      P_trait2      = pmax(P_trait2, 1e-300),
      log10p_trait1 = -log10(P_trait1),
      log10p_trait2 = -log10(P_trait2)
    )
  
  coloc_snps <- as_tibble(coloc_res$results) %>%
    mutate(snp = str_to_lower(snp)) %>%
    select(snp, SNP.PP.H4)
  
  wide <- wide %>%
    mutate(snp_lower = str_to_lower(snp)) %>%
    left_join(coloc_snps, by = c("snp_lower" = "snp")) %>%
    select(-snp_lower) %>%
    arrange(SNP.PP.H4)  # low PP.H4 first, high on top
  
  lead_vitiligo <- wide$snp[which.max(wide$log10p_trait1)]
  s    <- coloc_res$summary
  PPH0 <- formatC(s["PP.H0.abf"], digits = 3, format = "f")
  PPH1 <- formatC(s["PP.H1.abf"], digits = 3, format = "f")
  PPH2 <- formatC(s["PP.H2.abf"], digits = 3, format = "f")
  PPH3 <- formatC(s["PP.H3.abf"], digits = 3, format = "f")
  PPH4 <- formatC(s["PP.H4.abf"], digits = 3, format = "f")
  
  subtitle_txt <- paste0(
    "lead vitiligo SNP: ", lead_vitiligo, "\n",
    "PP.H0=", PPH0, " | PP.H1=", PPH1, " | PP.H2=", PPH2,
    " | PP.H3=", PPH3, " | PP.H4=", PPH4
  )
  title_txt     <- paste0(gene_name, " locus  |  ", trait1_name, " × ", trait2_name)
  
  col_scale <- scale_color_viridis_c(
    option = "plasma", limits = c(0, 1),
    oob = scales::squish, na.value = "grey80", name = "SNP PP.H4"
  )
  
  p1 <- ggplot(wide, aes(x = DIST, y = log10p_trait1, color = SNP.PP.H4)) +
    geom_point(alpha = 0.7, size = 1.5) + col_scale +
    labs(title = title_txt, subtitle = subtitle_txt,
         x = NULL, y = paste0("-log10(p) ", trait1_name)) +
    theme_bw(base_size = 12) + theme(legend.position = "top")
  
  p2 <- ggplot(wide, aes(x = DIST, y = log10p_trait2, color = SNP.PP.H4)) +
    geom_point(alpha = 0.7, size = 1.5) + col_scale +
    labs(x = NULL, y = paste0("-log10(p) ", trait2_name)) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  
  p3 <- ggplot(wide, aes(x = DIST, y = SNP.PP.H4, color = SNP.PP.H4)) +
    geom_point(alpha = 0.7, size = 1.5) + col_scale +
    labs(x = "Distance to TICAM1 midpoint (bp)", y = "SNP PP.H4") +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  
  p1 / p2 / p3 + plot_layout(heights = c(1, 1, 0.8))
}

# all traits
for (trait_name in names(results)) {
  cat(sprintf("Plotting vitiligo × %s...\n", trait_name))
  
  if (is.null(trait_dfs[[trait_name]])) {
    f <- file.path(gwas_dir, traits[[trait_name]])
    trait_dfs[[trait_name]] <- read_locus(f, TICAM1_chr, TICAM1_start, TICAM1_end)
  }
  
  p <- plot_gwas_gwas_locus(
    df1         = vit,
    df2         = trait_dfs[[trait_name]],
    trait1_name = "vitiligo",
    trait2_name = trait_name,
    coloc_res   = results[[trait_name]]
  )
  
  out_file <- file.path(out_dir, paste0("TICAM1_coloc_vitiligo_vs_", trait_name, ".pdf"))
  ggsave(out_file, p, width = 9, height = 10)
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\nAll plots done! 🎉\n")

# ============================================================
# plot all comparisons and save PDF
# ============================================================
trait_dfs <- list()  

for (trait_name in names(results)) {
  cat(sprintf("Plotting vitiligo × %s...\n", trait_name))
  
  if (is.null(trait_dfs[[trait_name]])) {
    f <- file.path(gwas_dir, traits[[trait_name]])
    trait_dfs[[trait_name]] <- read_locus(f, TICAM1_chr, TICAM1_start, TICAM1_end)
  }
  
  p <- plot_gwas_gwas_locus(
    df1          = vit,
    df2          = trait_dfs[[trait_name]],
    trait1_name  = "vitiligo",
    trait2_name  = trait_name,
    coloc_res    = results[[trait_name]]
  )
  
  out_file <- file.path(out_dir, paste0("TICAM1_coloc_vitiligo_vs_", trait_name, ".pdf"))
  ggsave(out_file, p, width = 9, height = 10)
  cat(sprintf("  Saved: %s\n", out_file))
}

cat("\nAll plots done! 🎉\n")