setwd("/nl/umw_manuel_garber/human/skin/eQTLs/ATAC-seq/fastqs/report3336/atac")

dat=fread("/nl/umw_manuel_garber/human/skin/eQTLs/ATAC-seq/fastqs/report3336/atac/filt_atac_annotated_snp.tsv",header=T)

unique(dat$F22)
dat=dat[,c(1:6,55:61,7:54)]
ceil(nrow(dat)/10000)
k=seq(from=1,to=nrow(dat),by = 10000)
k=c(k,nrow(dat))
collapsed_snps=c()
for(i in 1:(length(k)-1))
{
  x=melt(dat[c(k[i]:(k[i+1]-1)),],id=1:17)
  y=x %>% group_by(CHROM,POS,ID,REF,ALT,chr,start,end,FIB,KRT,MEL,region,value) %>% summarise(freq=n(),donors = paste(variable, collapse = ','))
  x=y %>% pivot_wider(id_cols = c("CHROM","POS","ID","REF","ALT","chr","start","end","FIB","KRT","MEL","region"),names_from = value,values_from = c(freq,donors))
  collapsed_snps=rbind(collapsed_snps,x)
  print(i)
}
rm(x)
rm(y)

