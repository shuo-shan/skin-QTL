
counts_table=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_08172022/analysis/batch-corrected-counts.txt")
cpm_table=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_08172022/analysis/TMM-normalized-batch-corrected-CPM.txt")
metadata2=read.csv("/nl/umw_manuel_garber/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_08172022/analysis/metadata.txt",sep="\t")

# step 1. load batch-corrected-counts of selected genes, metadata, and SNP genotype tables.
# 1a. import gene count matrix used for modeling.
counts2=counts_table["IRF4",]
samples.pbs=colnames(counts2)[grepl("PBS",colnames(counts2))]
samples.ifn=colnames(counts2)[grepl("IFN",colnames(counts2))]
counts2=counts2 %>% as.data.frame() %>% 
  dplyr::filter(rowSds(as.matrix(.[,samples.pbs])) > 0) %>% 
  dplyr::filter(rowSds(as.matrix(.[,samples.ifn])) > 0) 
adjusted_counts.pbs <- counts2 %>% as.data.frame() %>% dplyr::select(contains("PBS"))
adjusted_counts.ifn <- counts2 %>% as.data.frame() %>% dplyr::select(contains("IFN")) 

# 1b. import CPM matrix used for plotting.
CPM2=cpm_table["IRF4",]
CPM2=CPM2[rownames(counts2),]

# 1c. import genotype for snp of interest
snp.table=read.table("/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/rs12203592_genotype_hg19.bed",header=TRUE) 
snp.table$ALT=gsub("TRUE","T",snp.table$ALT)
View(t(snp.table[,c("CHROM","START","END","ID","REF","ALT","F25","F49","F55")]))
snp=snp.table[1,"ID"]
snp.ref=snp.table[1,"REF"]
snp.alt=snp.table[1,"ALT"] %>% gsub("TRUE","T",.)

genotype=snp.table %>% 
  dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
  rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))
genotype$donor=genotype$donor %>% gsub(".GT","", . )
genotype$genotype = genotype$genotype %>%
  gsub("0/0",0, . ) %>%
  gsub("0/1",1, . ) %>%
  gsub("1/1",2, . ) %>%
  gsub("./.",NA, . )

# 1d. import donor metadata
metadata2=metadata2 %>% right_join(genotype, . ,by="donor")

mod_test_table=data.frame()
g="IRF4"

df.pbs <- t(adjusted_counts.pbs[g,]) %>% as.data.frame() %>% set_colnames("expression") %>%
  dplyr::mutate(sample=rownames(.)) %>% left_join( . , metadata2, by="sample")
df.ifn <- t(adjusted_counts.ifn[g,]) %>% as.data.frame() %>% set_colnames("expression") %>% 
  dplyr::mutate(sample=rownames(.)) %>% left_join( . , metadata2, by="sample")

# build model for pbs
model.pbs <- glm.nb(expression ~ genotype, data=df.pbs)
beta1=summary(model.pbs)$coefficients[2,1]
se1=summary(model.pbs)$coefficients[2,2]
p1=summary(model.pbs)$coefficients[2,4]
# build model for ifn
model.ifn <- glm.nb(expression ~ genotype, data=df.ifn)
beta2=summary(model.ifn)$coefficients[2,1]
se2=summary(model.ifn)$coefficients[2,2]
p2=summary(model.ifn)$coefficients[2,4]
# then z-test on genotype coefficient beta
z=(beta1-beta2)/sqrt(se1^2 + se2^2)
pval=2*pnorm(-abs(z))

# prepare paired cpm value table. each donor is a row. very important table, used in next step.
this.cpm=t(CPM2[g,]) %>% as.data.frame() %>% rownames_to_column("sample")  %>%
  set_colnames(c("sample","batch_corrected_cpm")) %>% 
  left_join( . ,metadata2,by="sample") %>%
  mutate(condition=ordered(.$condition, levels=c("PBS","IFN"))) %>%
  filter(!is.na(genotype))
this.cpm$genotype=this.cpm$genotype %>%
  gsub("0",paste0(snp.ref,snp.ref),.) %>%
  gsub("1",paste0(snp.ref,snp.alt),.) %>%
  gsub("2",paste0(snp.alt,snp.alt),.) 
this.cpm$genotype=ordered(this.cpm$genotype,levels=c(paste0(snp.ref,snp.ref),paste0(snp.ref,snp.alt),paste0(snp.alt,snp.alt)))
this.cpm.pbs=this.cpm %>% dplyr::filter(condition=="PBS")
this.cpm.ifn=this.cpm %>% dplyr::filter(condition=="IFN")
this.cpm.paired= pivot_wider( this.cpm[,c("donor","genotype","batch_corrected_cpm","condition")], 
                              names_from="condition", values_from="batch_corrected_cpm")

this.metrics=this.cpm.paired
this.metrics$diff=this.cpm.paired$IFN - this.cpm.paired$PBS
this.metrics$fc=this.cpm.paired$IFN / this.cpm.paired$PBS
this.metrics$pctChange=(abs(this.metrics$diff) / this.metrics$PBS) * 100

# classify the IFNg responses by repsonse direction: up or down. Value of 0 will be considered "mild" later on.
this.metrics$change=this.metrics$diff %>% as.character()
this.metrics$change[this.metrics$diff < 0]<-"down"
this.metrics$change[this.metrics$diff > 0]<-"up"
# responses less than 25% of the PBS value are considered to be "mild".
this.metrics$change[this.metrics$pctChange <= 25] <- "mild" 
this.metrics$genotype=gsub("TRUE","T",this.metrics$genotype)

pdf("/nl/umw_manuel_garber/human/skin/eQTLs/DREG/MEL/rs12203592_IRF4_MEL.pdf")
p<-ggpaired(this.metrics,cond1="PBS",cond2="IFN",facet.by="genotype",
            color="condition",palette="aaas",line.color = "gray", line.size = 0.4,
            title=paste0(snp,":",g),width=0,
            xlab="condition",ylab="batch-corrected TMM-normalized CPM, melanocytes")
p$layers<-p$layers[-1] # remove geom_boxplot layer 
print(p)
dev.off()

