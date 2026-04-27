library("haploR")
library(data.table)
library(dplyr)
library(tidyr)
setwd("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping")
gwas=fread("/project/umw_garberlab/vangala/bin/gwas_catalog_v1.0.2-associations_e96_r2019-11-21.tsv",header = T,sep="\t")
gwas=gwas[grep("rs",gwas$SNPS,invert = F),]
snps=fread("tmp",sep="\t",header = T)
length(unique(snps$`DISEASE/TRAIT`))

gwas_snps=gwas[gwas$`DISEASE/TRAIT`%in%snps$`DISEASE/TRAIT`,c(8,22)]
#gwas_snps=inner_join(gwas_snps,snps[,c(1,2)],by=c("DISEASE/TRAIT"="V1"))
length(unique(gwas_snps$`DISEASE/TRAIT`))

gwas_snps=gwas_snps %>% separate_rows(SNPS)
gwas_snps=gwas_snps[grep("rs",gwas_snps$SNPS),]
gwas_snps=gwas_snps[!duplicated(gwas_snps),]
#write.table(gwas_snps,"gwas_skinDiseases",sep="\t",row.names = F,col.names = F,quote = F)

s=unique(gwas_snps$SNPS)
ld=c()
skip_to_next=F
i=0
k=c(seq(1,length(s),by = 100),length(s))
for(i in 2:length(k))
{
  skip_to_next=F
  tryCatch(queryHaploreg(query=gwas_snps$SNPS[c(k[i-1]:k[i])]), error =function(e)  { skip_to_next <<- TRUE})
  if(skip_to_next) { next }else{
    l=queryHaploreg(query=gwas_snps$SNPS[c(k[i-1]:k[i])])
    ld=rbind(l,ld)
    print(i)
  }
}
ld_disease=inner_join(ld,gwas_snps,by=c("query_snp_rsid"="SNPS"))

rm(gwas,snps,k)

dbSNP=fread("gunzip -c /project/umw_garberlab/vangala/bin/common_all_20180423.vcf.gz")
ld_hg19=dbSNP[dbSNP$ID%in%ld$rsID,]
rm(dbSNP)
ld_hg19_dis=inner_join(ld_hg19[,c(1:3)],ld_disease,by=c("ID"="rsID"))
ld_hg19_dis$`#CHROM`=gsub("^","chr",ld_hg19_dis$`#CHROM`)
ld_hg19_dis=ld_hg19_dis[!duplicated(ld_hg19_dis),]
write.table(ld_hg19_dis[,c(1:2,2:38)],"ld_snp_hg19_skin_diseases.tsv",row.names = F,sep="\t",quote=F)

######### GWAS skin 
skin=scan("/project/umw_garberlab/vangala/bin/GWAS_skin.txt",what="character",sep="\t")
skin=tolower(skin)
gwas$MAPPED_TRAIT=tolower(gwas$MAPPED_TRAIT)
skin_gwas=gwas[gwas$MAPPED_TRAIT%in%skin,]
skin_gwas=skin_gwas[grep("NR",skin_gwas$`RISK ALLELE FREQUENCY`,invert = T),]
skin_gwas$`RISK ALLELE FREQUENCY`=as.numeric(skin_gwas$`RISK ALLELE FREQUENCY`)
ggplot(skin_gwas[skin_gwas$MAPPED_TRAIT=="vitiligo",],aes(`RISK ALLELE FREQUENCY`))+geom_histogram()+geom_vline(xintercept =0.13)
