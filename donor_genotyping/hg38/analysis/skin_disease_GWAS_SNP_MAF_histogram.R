# crystal shan 06/12/2022
# create histogram for minor allele frequency of skin disease SNPs

library(dplyr)
library(magrittr)
library(tidyverse)
library(ggplot)
library(ggpubr)

f="/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/skin_disease_GWAS_SNPs_all_donors_merged.MAF.txt"
df=read.table(f,header=FALSE)
colnames(df)=c("position","variant","MAF")

gghistogram(df, x = "MAF", fill = "lightgray",add = "mean", rug = TRUE)

hist(df$MAF,breaks=20,main="MAF of 2,592 skin disease GWAS SNPs",xlab="MAF")
abline(v=0.135,col="royalblue4",lwd=2) # 5 donors out of 37 have homozygous minor allele: 5/37=0.135. 1,674 SNPs satisfy that.
abline(v=0.167,col="royalblue2",lwd=2) # 5 donors out of 30 have homozygous minor allele: 5/30=0.167. 1,538 SNPs satisfy this.
abline(v=0.25,col="royalblue",lwd=2) # 5 donors out of 20 have homozygous minor allele: 5/25=0.25. 1,151 SNPs.


f="/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/skin_disease_GWAS_SNPs_all_donors_merged.MAF.txt"
f=
  df=read.table(f,header=FALSE)
colnames(df)=c("position","variant","MAF")

gghistogram(df, x = "MAF", fill = "lightgray",add = "mean", rug = TRUE)

hist(df$MAF,breaks=20,main="MAF of 2,592 skin disease GWAS SNPs",xlab="MAF")
abline(v=0.135,col="royalblue4",lwd=2) # 5 donors out of 37 have homozygous minor allele: 5/37=0.135. 1,674 SNPs satisfy that.
abline(v=0.167,col="royalblue2",lwd=2) # 5 donors out of 30 have homozygous minor allele: 5/30=0.167. 1,538 SNPs satisfy this.
abline(v=0.25,col="royalblue",lwd=2) # 5 donors out of 20 have homozygous minor allele: 5/25=0.25. 1,151 SNPs.