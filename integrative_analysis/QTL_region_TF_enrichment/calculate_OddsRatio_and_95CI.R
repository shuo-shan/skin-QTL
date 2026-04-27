#!/usr/bin/env Rscript
# written by shuo.shan@umassmed.edu 03/2024

# Input arguments for counts from the contingency table
args = commandArgs(trailingOnly=TRUE)
qtl_in_region <- as.numeric(args[1])
qtl_not_in_region <- as.numeric(args[2])
nonqtl_in_region <- as.numeric(args[3])
nonqtl_not_in_region <- as.numeric(args[4])
QTL_type <- as.character(args[5])
region_type <- as.character(args[6])
output_file <- as.character(args[7])

#qtl_in_region <- 10
#qtl_not_in_region <- 476
#nonqtl_in_region <- 9709
#nonqtl_not_in_region <- 1426632
# QTL_type <- "MEL_reQTL"
# region_type <- "enhancer"
# output_file <- "~/Downloads/nl/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment/new_ODDS_RATIO_95CI.txt"
# 


# Calculating Odds Ratio (OR)
OR <- (qtl_in_region * nonqtl_not_in_region) / (qtl_not_in_region * nonqtl_in_region)

# Calculating the 95% Confidence Interval for the OR
# Standard error
SE_log_OR <- sqrt((1/qtl_in_region) + (1/qtl_not_in_region) + (1/nonqtl_in_region) + (1/nonqtl_not_in_region))

# Z-score for 95% confidence
Z <- 1.96

# Calculating the log of OR and Confidence Interval
log_OR <- log(OR)

# Calculating Confidence interval
lower_CI <- exp(log_OR - Z*SE_log_OR)
upper_CI <- exp(log_OR + Z*SE_log_OR)

# perform fisher's exact test
counts <- matrix(
  c(qtl_in_region, qtl_not_in_region,       # QTL counts (in region, not in region)
    nonqtl_in_region, nonqtl_not_in_region # nonQTL counts (in region, not in region)
  ),
  nrow = 2,         # Number of rows in the matrix
  byrow = TRUE,     # Fill matrix by rows
  dimnames = list(
    c("QTL", "nonQTL"),               # Row names
    c("In_Region", "Not_in_Region")   # Column names
  )
)
fisher_result <- fisher.test(counts)

# write to file
results <- data.frame(
  QTL_Type = QTL_type,
  Region_Type = region_type,
  OR = round(OR, 2),
  Lower_CI = round(lower_CI, 2),
  Upper_CI = round(upper_CI, 2),
  QTL_in_region = qtl_in_region,
  QTL_notIn_region = qtl_not_in_region,
  nonQTL_in_region = nonqtl_in_region,
  nonQTL_notIn_region = nonqtl_not_in_region,
  fisher_test_pval = fisher_result$p.value
)
data.table::fwrite(results,output_file,append=TRUE,quote=F,sep="\t",col.names =F)
