#!/usr/bin/env Rscript
# written by Crystal Shan 02/2022
cat("heya!!!! =D, loading libraries now... \n")
# import packages
#suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))

cat("starting to run the code... Good luck!!! \n")

# parsing arguments
args = commandArgs(trailingOnly=TRUE)

Dir=args[1] # make sure path doesn't have the "/" in the end
genotypeF=args[2] # genotype bed file (genotype and dosage)
headerF=args[3]
outF=paste0("new.",genotypeF) # outputs a new.xxx file.

dfName=paste0(Dir,"/",genotypeF)
df=read.table(dfName,header=FALSE,sep="\t")
colnames(df)=read.table(paste0(Dir,"/",headerF),header=TRUE,sep="\t") %>% colnames(.)
df=df %>% distinct(.,ID,.keep_all = TRUE)

snp.genotype=df %>% dplyr::select(c("ID",ends_with(".GT"))) %>% column_to_rownames("ID")
snp.anno=df %>% dplyr::select(!c(ends_with(".GT")))
  
output=data.frame()
for (i in 1:nrow(snp.genotype)) {
  this.row=snp.genotype[i,]
  df=this.row %>% t() %>% table %>% as.data.frame()
  colnames(df)=c("genotype","Freq")
  df=df[which(df$genotype!="./."),] # remove unknown genotypes
  df=df %>% arrange(desc(Freq)) # rearrange by top Frequency genotype
  has_3_genotypes=(nrow(df)==3) 
  has_at_least_2_genotypes=(nrow(df)>=2)
  is_well_spread=(min(df[,2])>=3) # min frequency to be 3
  is_balanced=(abs((df[1,2]-df[2,2])/max(df[1:2,2]))<0.5) # 1st highest and 2nd highest genotype freq don't differ by 50%
  if (has_at_least_2_genotypes && is_well_spread && is_balanced) {
    output=rbind(output,this.row)
  }
}

output=output %>% rownames_to_column("ID")

output2=right_join(snp.anno,output,by="ID")
colnames(output2)=gsub(".GT","",colnames(output2))

write.table(output2,file=paste0(Dir,"/",outF),
            sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

cat("finished processing file ",genotypeF," , have a nice day! \n")
