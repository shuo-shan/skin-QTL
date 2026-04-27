library(dplyr)
Dir="/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/analysis/pca"
setwd(Dir)

## Load data
dist_populations<-read.table("dataForPCA.mdist",header=F)
ancestry<-read.table("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/ancestry.txt",header=F)
colnames(ancestry)=c("donor","ancestry")
### Extract breed names
fam <- data.frame(famids=read.table("dataForPCA.mdist.id")[,1])
### Extract individual names 
famInd <- data.frame(IID=read.table("dataForPCA.mdist.id")[,2])

## Perform PCA using the cmdscale function
mds_populations <- cmdscale(dist_populations,eig=T,10)

## Extract the eigen vectors
eigenvec_populations <- cbind(fam,famInd,mds_populations$points)
eigenvec_populations <- left_join(ancestry,eigenvec_populations,by=c("donor"="IID")) %>% filter(donor!="F47")
data.table::fwrite(eigenvec_populations,"pca_eigenvec_table.txt",sep="\t",row.names=FALSE, col.names = TRUE)

## Proportion of variation captured by each eigen vector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)

# plot PCA
pdf("hg38_genotype_PCA_PC1n2.pdf")
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`,color = ancestry), show.legend = TRUE ) + 
  labs(title = "PCA of 37 donor genotype in skin-eQTL project",
       x = paste0("Principal component 1 (",eigen_percent[1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2]," %)")) + 
  geom_text(aes(x = `1`, y = `2`,label=donor), size=1) +
  theme_minimal() + theme(legend.position="bottom")
dev.off()

pdf("hg38_genotype_PCA_PC3n4.pdf")
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `3`, y = `4`,color = ancestry), show.legend = TRUE ) + 
  labs(title = "PCA of 37 donor genotype in skin-eQTL project",
       x = paste0("Principal component 3 (",eigen_percent[3]," %)"),
       y = paste0("Principal component 4 (",eigen_percent[4]," %)")) + 
  geom_text(aes(x = `3`, y = `4`,label=donor), size=1) +
  theme_minimal() + theme(legend.position="bottom")
dev.off()

pdf("hg38_genotype_PCA_PC5n6.pdf")
ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `5`, y = `6`,color = ancestry), show.legend = TRUE ) + 
  labs(title = "PCA of 37 donor genotype in skin-eQTL project",
       x = paste0("Principal component 5 (",eigen_percent[5]," %)"),
       y = paste0("Principal component 6 (",eigen_percent[6]," %)")) + 
  geom_text(aes(x = `5`, y = `6`,label=donor), size=1) +
  theme_minimal() + theme(legend.position="bottom")
dev.off()