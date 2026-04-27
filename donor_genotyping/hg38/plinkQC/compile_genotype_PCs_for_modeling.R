library(dplyr)

indir="/Users/crystal/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata"
name <- "data.only" # Because your files are test.bed, test.bim, test.fam
path2plink <- "/Users/crystal/Downloads/plink_mac_20250615/plink"

# list of PC files of my pruned data: (all donors)
# /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata/data.only.eigenval,
# /pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata/data.only.eigenvec
# Format .eigenvec into usable dataframe in R
pca_df <- read.table(file.path(indir, paste0(name, ".eigenvec")), header = FALSE)
colnames(pca_df) <- c("FID", "donor", paste0("PC", 1:10))

# If FID/IID are the same as donor IDs, you can drop one or rename as needed. 
# here subset to the list of donors for celltype of interest.
pca_df$donor_num <- pca_df$donor %>%
  gsub("skineQTL-", "", .) %>%
  gsub("^F0","", .) %>%
  gsub("^F", "", .)
pca_df$donor <- paste0("F",pca_df$donor_num)
pca_df <- pca_df[, c(c("FID", "donor", "donor_num", paste0("PC", 1:10)))]

# write table out
data.table::fwrite(pca_df, file="~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/genotype_PCs_for_modeling_07242025.txt", 
                   quote=F, sep="\t")