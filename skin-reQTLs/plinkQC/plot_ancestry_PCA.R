library(plinkQC)
indir="/Users/crystal/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata"
name <- "data" # Because your files are test.bed, test.bim, test.fam
path2plink <- "/Users/crystal/Downloads/plink_mac_20250615/plink"
refname <- 'HapMapIII_CGRCh38'
prefixMergedDataset <- paste(name, ".", refname, sep="")

exclude_ancestry <-
  evaluate_check_ancestry(indir=indir, name=name,
                          prefixMergedDataset=prefixMergedDataset,
                          refSamplesFile=paste(indir, "/HapMap_ID2Pop.txt",
                                               sep=""), 
                          refColorsFile=paste(indir, "/HapMap_PopColors.txt",
                                              sep=""),
                          verbose=TRUE,
                          interactive=TRUE)
