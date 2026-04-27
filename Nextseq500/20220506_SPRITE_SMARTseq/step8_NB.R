# Written by Pranitha Vangala
# Adapted by Crystal Shan July 2022

# >>> ABOUT: 
#     This script takes pairwise interaction table of 5kb tiles annotated with sequencing coverage and GC content,
#     fit a negative binomial model on coverage, GC content, and distance between a pair of interacting tiles, to 
#     model the raw interaction counts.

# >>> FORMAT: 
#     region called _1, _p1, or _p -> it's the anchor point, i.e. the 5kb tile containing a promoter
#     region called _2, _p2        -> it's the region interacting with anchor. 
#                                     either a promoter (with gene name), enhancer (overlap k27ac peak), or none

# >>> MODEL: 
#     Negative Binomial model. two parameters: mu (mean) and theta (shape param). two rounds of model fitting.
#     Round #1: take 20% of data to train model. then use model to fit to all data.
#     Round #2: 

# >>> INPUT FILES:
# hg19.5kb.coverage_2-100.......read coverage and GC content for each 5kb tile, defined in /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq/step8_NBmodel.sh
# all_cis.......................pairwise interaction table, defined in /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq/step7_tile-contacts.sh
# k27ac_peaks_annotations_5k....5kb tiles with their annotations (promoter, enhancer, or none; defined by /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq/step7_tile-contacts.sh)
# 
# >>> OUTPUT FILES: 
# all_cis_2MB.rds...............rds file for all cis pair interactions within 2M bp distance
# all_cis_2MB...................text file of above
# NB_all........................text file, all pairs fitted by NegBinom model
# NB_sig_all....................text file, significant pairs fitted by NegBinom model with p.adj <= 0.05
# NB_sig_active.bedpe...........BEDPE file of promoter-active enhancer significant pairs for IGV viewing
# NB_sig_all.bedpe..............BEDPE file of all significant pairs for IGV viewing

### ******************************************************
indir="/nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams"
#indir="/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq"
outdir="/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq"
options(stringsAsFactors = F)
sort_paste=function(x,y){
  key= paste(sort(c(x, y)), collapse=":")
}
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
setwd(outdir)

### ******************************************************
# step 1. Get sprite read coverage in each 5kb tile. Keep tiles with sprite reads.
#bg=fread(paste0(indir,"/hg19.5kb.coverage_2-100"),header=T) # pranitha coverage for /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams
#bg=fread(paste0(indir,"/hg38.5kb.coverage_2-100"),header=T) 
colnames(bg)=c("chrom_2","start_2","end_2","bg","GC")
bg <- bg %>% dplyr::filter(bg>100 & bg<10000) # require coverage to be between 100 and 10,000
#summary(bg) # uncomment this during interactive session
#head(bg) # uncomment this during interactive session
bg$id=paste(bg$chrom_2,bg$start_2,bg$end_2,sep="_")

### ******************************************************
# step 2. Get pairwise interaction information between two 5kb tiles. 
# Then refine pairwise interaction table. Only keep tiles that have enough coverage.
cis=fread(paste0(outdir,"/all_cis"))
colnames(cis)=c("chrom_2","start_2","end_2","gene_2","rawCounts","m","id_p1")
cis$id_p2=paste(cis$chrom_2,cis$start_2,cis$end_2,sep="_")
cis=cis[,c("chrom_2","start_2","end_2","gene_2","rawCounts","id_p1","id_p2","m")]

cis_p1=inner_join(cis,bg[,c("bg","GC","id")],by=c("id_p1"="id"))
colnames(cis_p1)=c(colnames(cis),"bg_p1","GC_p1")

cis_p2=inner_join(cis_p1,bg[,c("bg","GC","id")],by=c("id_p2"="id"))
colnames(cis_p2)=c(colnames(cis_p1),"bg_p2","GC_p2")

cis_final=cis_p2[,c("chrom_2","start_2","end_2","gene_2","rawCounts",
                    "id_p1","id_p2","bg_p1","GC_p1","bg_p2","GC_p2","m")]
rm(cis,cis_p1,cis_p2)

cis_final=cis_final %>% mutate(k=cis_final$id_p1) %>% tidyr::separate(k,c("chr_p","start_p","end_p"),sep="_")
#cis$jointD=as.numeric(cis$depth_p1) * as.numeric(cis$bg) # commented out by PV in original script. kept for bookkeeping.
cis_final$start_p=as.numeric(as.character(cis_final$start_p)) # start of promoter
cis_final$end_p=as.numeric(as.character(cis_final$end_p)) # end of promoter
cis_final$center=cis_final$start_p+round((cis_final$end_p-cis_final$start_p)/2,0) # center of promoter
cis_final$d=abs(cis_final$center-cis_final$start_2)-round((cis_final$end_p-cis_final$start_p)/2,digits = 0) # distance between inner edges of promoter and enhancer.
#summary(cis_final$d) # uncomment this during interactive session
print("pre-processing done; Fitting negative binomial model")

# limiting interactions to 2Mbp window, also remove the 5kb bins that are adjacent to the promoter bin.
cis_s=cis_final[which(cis_final$d<=2000000 & cis_final$d >5000),]
rm(bg,cis_final)

colnames(cis_s)=c("chrom_2","start_2","end_2","gene_2","rawCounts",
                  "id_p1","id_p2","bg_p1","GC_p1","bg_p2","GC_p2",
                  "chr_p","start_p","end_p","center","d","m")
cis_s=cis_s[,c(1:11,13:17,12)]
cis_s=cis_s %>%
  dplyr::rowwise() %>%
  mutate(key= paste(sort(c(id_p1, id_p2)), collapse=":")) %>% distinct(key,.keep_all=T)

saveRDS(cis_s,paste0(outdir,"/all_cis_2MB.rds"))
write.table(cis_s,paste0(outdir,"/all_cis_2MB"),sep="\t",quote = F,row.names = F,col.names = T)

### ******************************************************
### Step 3.  fit Negative Binomial Model
cis_s=readRDS(paste0(indir,"/all_cis_2MB.rds"))
cis_s=cis_s[sample(1:nrow(cis_s),100000),
            c("#chrom_2","start_2","end_2","gene_2","rawCounts",
              "id_p1","id_p2","bg_p1","GC_p1","bg_p2","GC_p2","chr_p",
              "start_p","end_p","center","d","m","key")] # I added this to run the script faster.

# create training set (20% of data)
s=sample(1:nrow(cis_s),nrow(cis_s)*0.2)
# uncomment this during interactive session
#k=sample(1:nrow(cis_s),10000) 
#ggplot(cis_s[k,],aes(GC_p1,bg_p1))+geom_point()+scale_x_log10()+scale_y_log10() # low coverage in GC poor regions.

# train model using 20% of data
mod_negB=MASS::glm.nb(rawCounts ~ log(bg_p2)+log(bg_p1)+log(d)+log(GC_p1)+log(GC_p2),data=cis_s[s,])

#summary(mod_negB) # uncomment this during interactive session

# predict raw interaction counts using the model for 100% of data
cis_s$pred_negB=predict(mod_negB,dispersion=1/mod_negB$theta,cis_s,type="response")

## taking the top scoring bins as true positives and refitting the model after removing
cis_s$outliers=(qnbinom(0.70, mu=cis_s$pred_negB,size=mod_negB$theta)) # get the predicted raw contacts at 0.70 probability
z_refit=cis_s[cis_s$rawCounts<cis_s$outliers,] # keep low-mid scoring bins to refit model. (all observed raw contacts lower than the Pr=0.70 probability highly fited data)
# ^ My understanding is refitting without true positives could prevent model from over-fitting.
s=sample(1:nrow(z_refit),nrow(z_refit)*0.2) # take 20% of the low-mid scoring samples to train new model.
mod_negB_re=MASS::glm.nb(rawCounts ~  log(bg_p2)+log(bg_p1)+log(d)+log(GC_p1)+log(GC_p2),data=z_refit[s,])
#summary(mod_negB_re) # uncomment this during interactive session

# refit all data using new model: pred_negB_reFit2 (refitted model predictions with low-mid scoring bins).
cis_s$pred_negB_reFit2=predict(mod_negB_re,cis_s,type="response",dispersion=1/mod_negB_re$theta)
# pvalue is calculated from the refitted model.
cis_s$pval_NB=1-pnbinom((cis_s$rawCounts),size=mod_negB_re$theta,mu=cis_s$pred_negB_reFit2,lower.tail = T)
cis_s$padj_NB=p.adjust(cis_s$pval_NB,method = "fdr")
cis_s$OE_NB=cis_s$rawCounts/cis_s$pred_negB_reFit2
cis_s$res=cis_s$rawCounts-cis_s$pred_negB_reFit2
cis_s$logpadj=-log10(round(cis_s$padj_NB,digits=15))
colnames(cis_s)[1]="#chrom_2"

cis_s$logpadj[is.infinite(cis_s$logpadj)]=max(cis_s$logpadj[!is.infinite(cis_s$logpadj)])+3 # add pseudocount to padj
cis_sig=cis_s[cis_s$padj_NB<=0.05,]
#cis_sig=cis_sig[cis_sig$OE_NB>2,] # commented out by PV in original script. kept for bookkeeping.
#cis_sig=cis_sig[cis_sig$rawCounts>4,] # commented out by PV in original script. kept for bookkeeping.

#table(cis_sig$gene_2) # uncomment this during interactive session
#table(cis_s$gene_2) # uncomment this during interactive session

k27ac=fread(paste0(indir,"/k27ac_peaks_annotations_5k")) # PV code
#k27ac=fread(paste0(indir,"/k27ac_peaks_annotations_5k.txt")) # Crystal code
k27ac= k27ac %>% group_by(V4) %>%
  summarise(V5=paste(V5, collapse=';'))
colnames(k27ac)=c("id_p1","element_p1")
cis_s=left_join(cis_s,k27ac,by="id_p1")
colnames(k27ac)=c("id_p2","element_p2")
cis_s=left_join(cis_s,k27ac,by=("id_p2"))
cis_s$element_p2[is.na(cis_s$element_p2)]="none"

#cis_s=cis_s[,c(1:4,28,6,27,5,7:16,18:26,17)] # v reorder columns to below
cis_s=cis_s[,c("#chrom_2","start_2","end_2","gene_2","element_p2",
               "id_p1","element_p1","rawCounts","id_p2","bg_p1","GC_p1","bg_p2","GC_p2",
               "chr_p","start_p","end_p","center","d",
               "key","pred_negB","outliers","pred_negB_reFit2","pval_NB","padj_NB",
               "OE_NB","res","logpadj","m")]
cis_sig=cis_s[cis_s$padj_NB<=0.05,]

write.table(cis_s,paste0(outdir,"/NB_all"),sep="\t",quote = F,row.names = F)
write.table(cis_sig,paste0(outdir,"/NB_sig_all"),sep="\t",quote = F,row.names = F)

### ******************************************************
# step 4. BEDPE file creation
# create BEDPE file for IGV viewing loops for active enhancers and promoters. (ignore promoter-interacting-regions that don't have k27ac peak)
# taken from script /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/washU.sh
bedpe.active=cis_sig %>% 
  dplyr::filter(element_p2!="none") %>%
  dplyr::select(c("#chrom_2","start_2","end_2","chr_p","start_p","end_p","logpadj")) %>%
  dplyr::mutate(key=paste0(`#chrom_2`,":",start_2,"-",end_2,"::",chr_p,":",start_p,"-",end_p)) %>%
  dplyr::mutate(extracol1="*",
                extracol2="*")
write.table(bedpe.active,paste0(outdir,"/NB_sig_active.bedpe"),sep="\t",quote = F,col.names = F ,row.names = F)


bedpe.active=cis_sig %>% 
  dplyr::select(c("#chrom_2","start_2","end_2","chr_p","start_p","end_p","logpadj")) %>%
  dplyr::mutate(key=paste0(`#chrom_2`,":",start_2,"-",end_2,"::",chr_p,":",start_p,"-",end_p)) %>%
  dplyr::mutate(extracol1="*",
                extracol2="*")
write.table(bedpe.active,paste0(outdir,"/NB_sig_all.bedpe"),sep="\t",quote = F,col.names = F ,row.names = F)
