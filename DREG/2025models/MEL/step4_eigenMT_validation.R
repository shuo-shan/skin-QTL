# validating eigenMT using 200 IFNG reQTL genes
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(data.table)
  library(magrittr)
  library(ggplot2)
})

ct="MEL"
dir <- paste0("/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/",ct,"/eigenMT/results")

eigenMT <- fread(paste0(dir,"/MEL_IFNG_reQTL.eigenMT.txt"), header=T)

# Choose a floor that's smaller than any observed non-zero p
p_nonzero_min <- min(eigenMT$p_gene_eigenMT[eigenMT$p_gene_eigenMT > 0], na.rm = TRUE)
p_floor <- p_nonzero_min / 10
eigenMT <- eigenMT %>%
  mutate(
    p_gene_eigenMT_safe = pmax(p_gene_eigenMT, p_floor),
    log10p = -log10(p_gene_eigenMT_safe)
  )
summary(eigenMT$log10p)

summary(eigenMT$Meff)
eigenMT <- eigenMT %>%
  mutate(
    sig_tier = case_when(
      p_gene_eigenMT_safe <= 1e-6 ~ "very_strong",
      p_gene_eigenMT_safe <= 1e-4 ~ "strong",
      p_gene_eigenMT_safe <= 1e-2 ~ "moderate",
      TRUE                        ~ "null"
    ),
    Meff_tier = case_when(
      Meff < 30            ~ "low_Meff",        # sparse LD
      Meff >= 30 & Meff < 45 ~ "mid_Meff",       # increasing LD
      Meff >= 45           ~ "high_Meff"   # LD ceiling
    )
  )

eigenMT$sig_tier <- ordered(eigenMT$sig_tier, levels=c("null","moderate","strong","very_strong"))
eigenMT$Meff_tier <- ordered(eigenMT$Meff_tier, levels = c("low_Meff","mid_Meff","high_Meff"))

eigenMT %>%
  dplyr::group_by(Meff_tier) %>%
  dplyr::summarise(
    n = n(),
    Meff_min = min(Meff, na.rm = TRUE),
    Meff_q25 = quantile(Meff, 0.25, na.rm = TRUE),
    Meff_median = median(Meff, na.rm = TRUE),
    Meff_q75 = quantile(Meff, 0.75, na.rm = TRUE),
    Meff_max = max(Meff, na.rm = TRUE)
  )

table(eigenMT$sig_tier, eigenMT$Meff_tier)


# sample 160 genes evenly across 12 strata
set.seed(42)
genes_per_stratum <- ceiling(200 / 12)

sampled_genes <- eigenMT %>%
  group_by(sig_tier, Meff_tier) %>%
  group_modify(~ {
    n_take <- min(nrow(.x), genes_per_stratum)
    dplyr::slice_sample(.x, n = n_take)
  }) %>%
  ungroup()

nrow(sampled_genes)
table(sampled_genes$sig_tier, sampled_genes$Meff_tier)
summary(sampled_genes$Meff)
summary(sampled_genes$p_gene_eigenMT_safe)

# ------- fetch permutation result table ------ ####
permutation.file <- "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/permutation/result/IFNG/reQTL/MEL_IFNG_reQTL.permutation_summary.txt"
permutation <- fread(permutation.file) %>%
  dplyr::select(gene, n_snps, perms, exceed, p_emp, stage) %>%
  set_colnames(c("gene","n_snps", "perms", "exceed","p_gene_permut","stage"))
strong_genes <- unique(permutation[which(permutation$stage=="stopped@B3"),]$gene)
permutation <- permutation[!(gene %in% strong_genes & stage == "needToRunB3")]

# -------- compare eigenMT vs. permutation result ------- ####
df <- left_join( sampled_genes, permutation, by="gene")
# ---- make safe p's for permutation too ----
p_perm_nonzero_min <- min(df$p_gene_permut[df$p_gene_permut > 0], na.rm = TRUE)
p_perm_floor <- p_perm_nonzero_min / 10

df <- df %>%
  mutate(
    p_gene_permut_safe = pmax(p_gene_permut, p_perm_floor),
    x = -log10(p_gene_eigenMT_safe),
    y = -log10(p_gene_permut_safe),
    diff = y - x   # >0 means perm more significant; <0 means eigenMT more significant
  )

ggplot(df, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(color = sig_tier), alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "-log10(p_gene_eigenMT)", y = "-log10(p_gene_permutation)") +
  xlim(0,15) + ylim(0,15) +
  theme_classic()


# calculate correlation coefficient
with(df, cor(x, y, method = "pearson", use = "complete.obs"))
with(df, cor(x, y, method = "spearman", use = "complete.obs"))
fit <- lm(y ~ x, data = df)
summary(fit)$coefficients
#If slope ≈ 1 and intercept ≈ 0, you’re golden.
#If slope < 1, eigenMT is often too aggressive at small p.

# plot by Meff bins
ggplot(df, aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(color = sig_tier), alpha = 0.8, size = 2) +
  facet_wrap(~Meff_tier) +
  labs(x = "-log10(p_gene_eigenMT)", y = "-log10(p_gene_permutation)") +
  theme_classic()

ggplot(df, aes(x = x, y = diff)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(color = Meff_tier), alpha = 0.8, size = 2) +
  labs(x = "-log10(p_gene_eigenMT)", y = "(-log10 perm) - (-log10 eigenMT)") +
  theme_classic()

# diff near 0 → agreement
# diff negative → eigenMT is more significant than permutation (anti-conservative)
# diff positive → eigenMT is less significant (conservative)


df_mid <- df %>%
  filter(y <= 4.2)  # or p_gene_permut > ~1e-4
with(df_mid, cor(x, y, method = "pearson"))
with(df_mid, cor(x, y, method = "spearman"))
df_mid %>%
  group_by(Meff_tier) %>%
  summarise(
    spearman = cor(x, y, method = "spearman"),
    n = n()
  )


# plotting rank instead
# eigenMT preserves ordering, which is what FDR depends on
df %>%
  mutate(
    rank_eigenMT = rank(x),
    rank_perm    = rank(y)
  ) %>%
  ggplot(aes(rank_eigenMT, rank_perm)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_classic()

# compare genes that pass FDR 0.05 in permutation vs. eigenMT
permutation$q_gene = p.adjust(permutation$p_gene_permut, method="BH")
eigenMT$q_gene = p.adjust(eigenMT$p_gene_eigenMT, method="BH")
genes_sig_by_permutation = permutation[which(permutation$q_gene<0.05),]$gene
genes_sig_by_eigenMT = eigenMT[which(eigenMT$q_gene<0.05),]$gene
genes_sig_by_permutation
genes_sig_by_eigenMT
