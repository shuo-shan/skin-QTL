Dir="/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/FRB"
setwd("/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/FRB")

p1=read.table(file=paste0(Dir,"/GSE63525_IMR90_HiCCUPS_liftover_hg19tohg38_p1.bed"),header=FALSE) %>% set_colnames(c("chr","start","end","id","extra"))
p2=read.table(file=paste0(Dir,"/GSE63525_IMR90_HiCCUPS_liftover_hg19tohg38_p2.bed"),header=FALSE) %>% set_colnames(c("chr","start","end","id","extra"))

pairs=inner_join(p1,p2,by="id") %>% 
  dplyr::select(-c("extra.x","extra.y","id"))

write.table(pairs,file=paste0(Dir,"/GSE63525_IMR90_HiCCUPS_liftover_hg19tohg38.bedpe"),
            quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
