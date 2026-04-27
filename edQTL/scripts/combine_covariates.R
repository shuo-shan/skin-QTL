library(tidyverse)

# parsing input arguments
args = commandArgs(trailingOnly=TRUE)
peerFile=args[1]
genotypePCFile=args[2]
outFile=args[3] # specify full file path

#peerFile="~/Downloads/nl/human/skin/eQTLs/edQTL/output/foreskin.edMat.10cov.20samps.noXYM.qqnorm.PEER_covariates.txt"
#genotypePCFile="~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/analysis/pca/pca_eigenvec_table.txt"
#outFile="~/Downloads/nl/human/skin/eQTLs/edQTL/output/combined_covariates.txt"

peer = read.table(peerFile, header=TRUE, sep="\t")

genotypePC = read.table(genotypePCFile, header=TRUE, sep="\t")
genotypePC.new = genotypePC[,4:ncol(genotypePC)]
rownames(genotypePC.new) = genotypePC[,1]
colnames(genotypePC.new) = gsub("X","gPC_",colnames(genotypePC.new))
genotypePC.new = t(genotypePC.new)
genotypePC.new = cbind(data.frame(ID=rownames(genotypePC.new)),genotypePC.new)
genotypePC.new = genotypePC.new[,colnames(peer)]

combined_table = rbind(peer, genotypePC.new)

write.table(combined_table, file=outFile,quote=FALSE, sep="\t", row.names=FALSE, col.names = TRUE)
