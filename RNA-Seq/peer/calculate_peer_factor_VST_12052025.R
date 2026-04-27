library(tidyverse)
library(DESeq2)
library(magrittr)
library(stringr)
library(gplots)
library(peer)

Dir="~/Downloads/nl/human/skin/eQTLs/RNA-Seq/peer"
setwd(Dir)
#save.image(file=paste0(Dir,"/myEnvironment_main.RData"))
load(file=paste0(Dir,"/myEnvironment_main.RData"))

#### ----- exclude outliers ----- ####
# F78 and F138 doesn't have genotype data (no more sample left)
# F86 or F104 might also need to be excluded since they failed heterozygosity check. TBD.
outliers_list <- c("F26_KRT_PBS","F39_KRT_PBS","F42_KRT_PBS","F43_KRT_PBS","F60_KRT_PBS","F116_KRT_PBS","F120_KRT_PBS",
                   "F26_KRT_IFNG","F39_KFT_IFNG","F42_KFT_IFNG","F43_KRT_IFNG","F60_KRT_IFNG","F112_KRT_IFNG",
                   "F103_KRT_IFNB","F118_KRT_TNF")
pattern <- str_c(outliers_list, collapse = "|")
outlier_samples_metadata.krt <- metadata %>% dplyr::filter(str_detect(sample, pattern))

outliers_list <- c("F108_MEL","F109_MEL","F110_MEL","F111_MEL","F51_MEL",
                   "F52_MEL","F61_MEL","F118_MEL","F133_MEL","F116_MEL_IFNG","F116_MEL_IFNG",
                   "F142_MEL_TNF","F42_MEL_PBS","F34_MEL_PBS","F30_MEL_IFNG")
lowcounts_list <- c("F116_MEL_IFNG","F142_MEL_TNF","F42_MEL_PBS","F34_MEL_PBS","F30_MEL_IFNG")
pattern <- str_c(outliers_list, collapse = "|")
outlier_samples_metadata.mel <- metadata %>% dplyr::filter(str_detect(sample, pattern))

genotype_outliers_list <- c("F86","F104","F78","F138")
pattern <- str_c(genotype_outliers_list, collapse = "|")
genotype_outlier_samples_metadata <- metadata %>% dplyr::filter(str_detect(sample, pattern))


outlier_samples_to_exclude <- unique(c(
  outlier_samples_metadata.krt$sample,
  outlier_samples_metadata.mel$sample,
  genotype_outlier_samples_metadata$sample
))

VST <- left_join(vst.mel, vst.krt, by="gene") %>%
  left_join( . , vst.frb, by = "gene") %>%
  column_to_rownames("gene") 

# 3. Identify outliers present in the current VST table
outliers_present <- intersect(outlier_samples_to_exclude, colnames(VST))

# 4. Filter the VST table to exclude only the detected outliers
VST <- VST %>%
  # The negative sign (-) on the list drops those columns
  dplyr::select(-all_of(outliers_present))

CPM <- CPM[, colnames(VST)]

metadata = metadata %>% dplyr::filter(sample %in% colnames(VST))

#rm(list = setdiff(ls(envir = .GlobalEnv), c("CPM","log2CPM","metadata","expressed_genes_list","Dir","scenarios","peer_factor_list")), envir = .GlobalEnv)

#### ----- testing if excluding outliers would cause modeling to fail ----- ####
snpinfo <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/2025models/test1/snp.bed") %>% t() %>% 
  as.data.frame() %>%
  rownames_to_column("donor") %>%
  set_colnames(c("donor","GT"))
snpinfo$donor <- gsub("skineQTL-","",snpinfo$donor)
REF <- "T"
ALT <- "C"
genotype <- snpinfo[7:nrow(snpinfo), ]
genotype <- genotype %>%
  mutate(
    genotype.num = case_when(
      GT == "0|0" ~ 0,
      GT %in% c("0|1", "1|0") ~ 1,
      GT == "1|1" ~ 2,
      TRUE ~ NA_real_
    ),
    genotype.nt = case_when(
      GT == "0|0" ~ "TT",
      GT == "0|1" ~ "TC|CT",
      GT == "1|0" ~ "TC|CT",
      GT == "1|1" ~ "CC",
      TRUE ~ NA_character_
    )
  )
genotype$donor <- sub("^F0", "F", genotype$donor)

erap2 <- data.frame(CPM = as.numeric(CPM["ERAP2",]), log2CPM = as.numeric(log2CPM["ERAP2",]), sample = colnames(CPM))

this_matrix <- left_join(erap2, metadata, by="sample") %>%
  left_join( . , genotype, by="donor") %>%
  dplyr::filter(celltype=="MEL") %>%
  dplyr::filter(condition%in%c("PBS","IFNG"))
this_matrix$condition <- factor(this_matrix$condition, levels = c("PBS", "IFNG"), ordered = FALSE)

model13 <- lmerTest::lmer(CPM ~ condition * genotype.num  +
                            (1 | donor), #random intercept for donor
                          data = this_matrix)
broom.mixed::tidy(model13)


model14 <- lmerTest::lmer(log2CPM ~ condition * genotype.num  +
                            (1 | donor), #random intercept for donor
                          data = this_matrix)
broom.mixed::tidy(model14)

this_matrix$fitted <- fitted(model13)
this_matrix$resid <- resid(model13)
this_matrix$genotype_condition <- interaction(this_matrix$genotype.num, this_matrix$condition, sep = "_")

ggplot(this_matrix, aes(x = fitted, y = resid, color = genotype_condition)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40") +
  labs(x = "Fitted values", y = "Residuals") +
  theme_minimal()

ggplot(this_matrix) +
  geom_point(aes(x=genotype.num, y=CPM)) +
  facet_wrap(~condition)

#### ----- PEER factor calculation ----- ####
# install peer
#devtools::install_github("belowlab/r-peer")

#### compile expressed genes (CPM>=1 in >=5 donors in any condition) ####
# Create all scenarios
scenarios <- c(
  "MEL_PBS-IFNG-IFNB-TNF",
  "KRT_PBS-IFNG-IFNB-TNF",
  "FRB_PBS-IFNG-IFNB-TNF"
)

# Initialize list for expressed genes
expressed_genes_list <- list()

# Function to get expressed genes for a set of samples
get_expressed_genes <- function(samples, CPM, min_cpm = 1, min_donors = 5) {
  if (length(samples) == 0) return(character(0))
  # Matrix of CPM >= 1
  expr_mat <- CPM[, samples, drop = FALSE] >= min_cpm
  # Keep genes with >=5 donors with CPM >= 1
  rownames(CPM)[rowSums(expr_mat) >= min_donors]
}

# Loop through each scenario
for (sc in scenarios) {
  parts <- strsplit(sc, "_")[[1]]
  ct <- parts[1]
  
  if (grepl("-", sc)) {
    # Two-condition scenario (PBS and cytokine)
    conds <- strsplit(parts[2], "-")[[1]]
    cond_samples <- c()
    
    for (cond in conds) {
      cond_samples <- c(cond_samples, colnames(CPM)[grepl(paste0("_", ct, "_", cond), colnames(CPM))])
    }
    
    expressed_genes <- get_expressed_genes(cond_samples, CPM)
    
  } else {
    # Single-condition scenario
    cond <- parts[2]
    cond_samples <- colnames(CPM)[grepl(paste0("_", ct, "_", cond), colnames(CPM))]
    expressed_genes <- get_expressed_genes(cond_samples, CPM)
  }
  
  expressed_genes_list[[sc]] <- expressed_genes
}

# Check results
names(expressed_genes_list)
length(expressed_genes_list[["MEL_PBS-IFNG-IFNB-TNF"]])  # number of expressed genes in MEL_PBS

#### calculate peer factor ####
# Directory for saving results
outdir <- "peer_factors"
dir.create(outdir, showWarnings = FALSE)

# Number of factors
n_factors <- 20

peer_factor_list <- list()

for (sc in names(expressed_genes_list)) {
  message("Processing scenario: ", sc)
  
  # Get samples and genes for this scenario
  # Determine cell type and condition(s)
  if (grepl("-", sc)) {
    # Two-or-more-condition scenario, e.g., "FRB_PBS-TNF"
    ct <- strsplit(sc, "_")[[1]][1]    # e.g., "FRB"
    conds <- strsplit(strsplit(sc, "_")[[1]][2], "-")[[1]]  # c("PBS", "TNF")
    
    # Match samples from all condition
    # Create regex:  "^MEL_(PBS|IFNG|TNF|IFNB)"
    pattern <- paste0(ct, "_(IFN|", paste(conds, collapse="|"), ")")
    
    # Get matching samples
    samples <- grep(pattern, colnames(VST), value = TRUE)
  } else {
    # Single-condition scenario
    samples <- colnames(VST)[grepl(sc, colnames(VST))]
  }
  genes <- expressed_genes_list[[sc]]
  genes <- intersect(genes, rownames(VST))
  
  if (length(samples) < 5 || length(genes) < 50) {
    warning("Skipping ", sc, ": too few samples or genes.")
    next
  }
  
  # take top 3000 most variable genes
  expr_full <- VST[genes, samples, drop = FALSE]
  
  # Keep only variable genes (top 3000 by variance)
  vars <- apply(expr_full, 1, var)
  genes_use <- names(sort(vars, decreasing = TRUE))[1:min(3000, length(vars))]
  expr_mat <- expr_full[genes_use, , drop = FALSE]
  
  # Remove zero-variance genes
  keep_genes <- apply(expr_mat, 1, function(x) var(x, na.rm = TRUE) > 0)
  expr_mat <- expr_mat[keep_genes, , drop = FALSE]
  
  # **Standardize per gene (critical)**
  # prevents PEER from absorbing global IFNG/TNF shifts
  expr_mat <- as.matrix(expr_mat)
  expr_mat <- t(scale(t(expr_mat), center = TRUE, scale = TRUE))
  
  # Prepare covariates: one column per condition (no intercept)
  metadata_df <- as.data.frame(metadata)
  match_indices <- match(samples, metadata_df$sample)
  cond <- factor(metadata_df[match_indices, "condition"])
  
  cov <- model.matrix(~ 0 + cond)
  
  # Run PEER
  model <- PEER()
  PEER_setPhenoMean(model, t(expr_mat))   # samples x genes
  PEER_setCovariates(model, cov)          # condition fixed effects matrix
  PEER_setNk(model, n_factors)
  PEER_update(model)
  
  # Extract factors
  factors <- PEER_getX(model)
  n_factors_actual <- ncol(factors)
  colnames(factors) <- paste0("PEER", seq_len(n_factors_actual))
  rownames(factors) <- samples
  
  write.table(factors,
              file = file.path("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/peer/peer_factors", paste0("/peer_factors_", sc, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
  
  peer_factor_list[[sc]] <- factors
}

saveRDS(peer_factor_list, file.path("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/peer/peer_factor_list.rds"))

### plot scree plot similar to PCA
peer_factors <- peer_factor_list[["FRB_PBS-IFNG-IFNB-TNF"]]  # matrix: samples x PEER factors
sdev <- apply(peer_factors, 2, sd)
var_explained <- sdev^2 / sum(sdev^2)

plot(var_explained, type = "b", pch = 19,
     xlab = "PEER factor", ylab = "Fraction of variance explained",
     main = "Scree plot of PEER factors (FRB_PBS-IFNG-IFNB-TNF)")

