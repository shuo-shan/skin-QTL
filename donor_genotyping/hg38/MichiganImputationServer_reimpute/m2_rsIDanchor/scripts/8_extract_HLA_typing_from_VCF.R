# script to extract HLA typing information from postMIS genotype file of chr6
# author: carol lopez

library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

vcf_fp <- "~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/HLA/chr6.dose.vcf.gz"     # tabix-index with: bcftools index -t my.vcf.gz

# 1) Read VCF (Michigan output is bgzipped + tabix indexed)
vcf <- readVcf(vcf_fp, genome="hg19")

# 2) Keep only HLA allele variants (exclude amino-acid markers that start with AA_)
ids <- rownames(vcf)                      # VariantAnnotation stores ID in rownames
is_hla_allele <- grepl("^HLA", ids) & !grepl("^AA_", ids)
table(is_hla_allele)
vcf_hla <- vcf[is_hla_allele, ]

# 3) Pull per-sample FORMAT: use DS if present; else build DS from GP (dosage = p(het)+2*p(homALT))
fmt <- geno(vcf_hla)
have_ds <- "DS" %in% names(fmt)
if (have_ds) {
  DS <- fmt$DS  # matrix: variants x samples
} else if ("GP" %in% names(fmt)) {
  # GP is 3 probs per sample: P(0/0), P(0/1), P(1/1)
  GP <- fmt$GP  # array: variants x samples x 3
  DS <- GP[,,2] + 2*GP[,,3]
} else if ("GT" %in% names(fmt)) {
  # Fallback from hard calls
  GT <- fmt$GT
  DS <- apply(GT, 2, function(col) {
    # Convert 0/0, 0|1 etc -> 0,1,2
    gsub_res <- gsub("[|/]", " ", col)
    sapply(strsplit(gsub_res, " "), function(a) sum(as.integer(a), na.rm=TRUE))
  })
  if (!is.matrix(DS)) DS <- matrix(DS, nrow=nrow(GT), ncol=ncol(GT))
  rownames(DS) <- rownames(GT); colnames(DS) <- colnames(GT)
} else {
  stop("VCF lacks DS, GP, and GT in FORMAT; cannot derive dosages.")
}

# # 4) Parse gene + allele from the variant ID
# # Michigan IDs commonly look like "HLA_A*01:01" or "HLA_A_01_01"
# parse_hla_id <- function(x) {
#   # Normalize underscores to star/colon if needed
#   # HLA_A_01_01 -> HLA_A*01:01
#   x2 <- gsub("^HLA_([A-Z0-9]+)_(\\d{2})_(\\d{2}).*$", "HLA_\\1*\\2:\\3", x)
#   x2 <- gsub("^HLA-([A-Z0-9]+)_(\\d{2})_(\\d{2}).*$", "HLA_\\1*\\2:\\3", x2)
#   # Already star/colon form stays as-is
#   x2
# }
# 
# ids_norm <- parse_hla_id(ids[is_hla_allele])

# 4) Normalize IDs to IMGT two-field format with HLA- prefix
normalize_hla <- function(x) {
  # Case 1: Michigan format with underscores, e.g. HLA_A_01_01 -> HLA-A*01:01
  x <- gsub("^HLA[_-]([A-Z0-9]+)_(\\d{2})_(\\d{2})(.*)$", "HLA-\\1*\\2:\\3\\4", x)
  
  # Case 2: Already HLA- prefixed, but may have extra fields
  x <- gsub("^HLA[_-]?([A-Z0-9]+)\\*", "HLA-\\1*", x)
  
  # Case 3: Missing prefix, e.g. A*02:06 -> HLA-A*02:06
  x <- ifelse(grepl("^[A-Z]+\\*", x), paste0("HLA-", x), x)
  
  # Drop fields beyond two-field resolution, e.g. HLA-A*02:06:01 -> HLA-A*02:06
  x <- sub("^(HLA-[A-Z0-9]+\\*\\d+:\\d+).*", "\\1", x)
  
  return(x)
}

ids_norm <- normalize_hla(ids[is_hla_allele])

# Extract gene (A,B,C,DRB1,...) and allele (e.g., A*01:01)
gene   <- sub("^HLA[_-]?([^*_]+).*", "\\1", ids_norm)
allele <- ids_norm
#allele <- sub("^HLA[_-]?([^*_]+)\\*(\\d+:\\d+).*", "\\1*\\2", ids_norm)

# 5) Build a long table of dosages per (sample, gene, allele)
ds_long <- as.data.frame(DS, check.names = FALSE) |>
  tibble::rownames_to_column("var_id") |>
  mutate(gene = gene, allele = allele) |>
  pivot_longer(-c(var_id, gene, allele),
               names_to = "sample", values_to = "dosage")

# 6) For each (sample, gene), pick the top 2 alleles by dosage
# Require at least 0.5 dosage to keep an allele
calls <- ds_long |>
  group_by(sample, gene) |>
  arrange(desc(dosage), .by_group = TRUE) |>
  filter(dosage >= 0.5) |>         # <-- threshold
  slice_head(n = 2) |>
  mutate(rank = row_number()) |>
  ungroup() |>
  select(sample, gene, rank, allele, dosage)

# # 7) Cast to wide: H1/H2 per gene
# calls_wide <- calls |>
#   mutate(field2 = sub("^([^*]+\\*\\d+:\\d+).*", "\\1", allele)) |>   # force two-field
#   mutate(label = ifelse(rank == 1, "H1", "H2")) |>
#   select(sample, gene, label, field2, dosage) |>
#   pivot_wider(names_from = label,
#               values_from = c(field2, dosage),
#               names_sep = ".")
# 7) Cast to wide: H1/H2 per gene
calls_wide <- calls |>
  mutate(label = ifelse(rank == 1, "H1", "H2")) |>
  select(sample, gene, label, allele, dosage) |>
  pivot_wider(names_from = label,
              values_from = c(allele, dosage),
              names_sep = ".")

# 8) Finally, one row per sample with all genes
hla_wide <- calls_wide |>
  arrange(sample, gene) |>
  pivot_wider(id_cols = sample,
              names_from = gene,
              values_from = c(allele.H1, allele.H2, dosage.H1, dosage.H2),
              names_sep = ".")

# Save if you like:
write_tsv(hla_wide, "~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/HLA/HLA_typing_from_Michigan.tsv")

hla_wide[1:5, ]  # peek
class(hla_wide)
library(dplyr)
library(purrr)
library(readr)

# 1) Turn any list-column into a single string (comma-separated). Empty -> NA.
hla_wide_flat <- hla_wide %>%
  mutate(across(where(is.list),
                ~ purrr::map_chr(.x, ~ if (is.null(.x) || length(.x) == 0) NA_character_
                                 else paste0(.x, collapse = ","))))

# 2) (Optional) ensure all are atomic vectors
hla_wide_flat <- data.frame(hla_wide_flat, check.names = FALSE)
dim(hla_wide_flat)
# 3) Write to disk
write_tsv(hla_wide_flat,
          "~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/HLA/HLA_typing_persample_Hg19_output.tsv"
)

#write.table(unlist(hla_wide),"~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/HLA/HLA_typing_persample_Hg19_output.tsv",sep="\t",quote=F, col.names=F)
