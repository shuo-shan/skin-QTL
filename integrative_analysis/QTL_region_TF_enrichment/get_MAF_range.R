library(vcfR)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
maf=args[1]

dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment"
vcf <- read.vcfR(paste0(dir,"/annotated_region.vcf.gz"))

lower_bound <- maf - maf * 0.2
upper_bound <- maf + maf * 0.2

# Assuming 'af' contains the allele frequencies calculated or extracted earlier
snps_within_range <- which(vcf@fix[,"AF"] >= lower_bound & vcf@fix[,"AF"] <= upper_bound)
df <- vcf@fix %>% tidyr::separate_rows(INFO,sep=";")
df <- as.data.frame(vcf@fix)
df$AF <- lapply(df$INFO, function(x) unlist(strsplit(x, ";"))[3]) %>% gsub("AF=","",.) %>% as.numeric()
snps_within_range <- which(df$AF >= lower_bound & df$AF <= upper_bound)

# Extract these SNPs
filtered_snps <- vcf[snps_within_range, ]
write.vcf(filtered_snps, file=paste0(dir,"/SNP_in_region_and_MAF_range.vcf"))