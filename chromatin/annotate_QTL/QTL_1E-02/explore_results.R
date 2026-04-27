library(tidyverse)
library(dplyr)
library(magrittr)

####### load data
# also load DREG/scripts/make_tiny_QTL_paired_plots.R plotting data and function.
dir="~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02"
bigtable.krt <- readRDS("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_krt.rds")
bigtable.mel <- readRDS("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_mel.rds")
bigtable.frb <- readRDS("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/bigtable_frb.rds")
betacomp.krt <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/KRT/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE)) %>%
  dplyr::filter(p.betaComp10KPermut<0.01)
betacomp.mel <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/MEL/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE)) %>%
  dplyr::filter(p.betaComp10KPermut<0.01)
betacomp.frb <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/FRB/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE)) %>%
  dplyr::filter(p.betaComp10KPermut<0.01)
ldtrait <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results/QTLs_LDtrait_association_ALLpopulation_r20.1.txt")
deg <- data.table::fread(paste0(dir,"/yuqing_disease_DEG/four_diseases_DEG_padj0.05_log2FC1.txt")) %>%
  set_colnames(c("gene","n_disease","disease"))
load.image("~/Downloads/nl/human/skin/eQTLs/website/data/plotting_data.rds")
#save.image("~/Downloads/skin_eQTL_website_data.rds")

#### pick top QTL for ERAP2 #### 
df <- bigtable.krt %>% dplyr::filter(gene=="ERAP2")

#### pick disease DE genes and top reQTL #### 
genes_exclude = c("ERAP2",deg$gene[grep("HLA-",deg$gene)])
df <- bigtable.krt %>% 
  dplyr::filter(gene %in% deg$gene) %>%
  dplyr::filter(!is.na(reQTL_pval)) %>%
  dplyr::filter(!gene %in% genes_exclude) %>%
  dplyr::filter(reQTL_pval < 0.0000001)

#### GBP3 fine-mapping #### 
df <- bigtable.krt %>% dplyr::filter(gene=="GBP3") %>% dplyr::filter(!is.na(reQTL_pval))

#### ITGA1 fine-mapping #### 
df <- bigtable.krt %>% dplyr::filter(gene=="ERAP2") %>% dplyr::filter(!is.na(reQTL_pval))

#### SEPHS2 fine-mapping #### 
df <- bigtable.mel %>% dplyr::filter(gene=="SEPHS2") %>% dplyr::filter(!is.na(reQTL_pval))


#### focus on SNPs that have high polyPhen score (possibly damaging to probably damaging) #### 
snp_list <- c("rs2235794","rs600377","rs12283300","rs7165988","rs4889244","rs9912644",
              "rs6728493","rs430665","rs923828","rs9912852")
df <- bigtable.mel %>% dplyr::filter(QTL %in% snp_list) %>% dplyr::filter(gene %in% deg$gene)


#### focus on SNPs that overlap the most TF motifs ####
df <- bigtable.frb %>% 
  dplyr::filter(gene %in% deg$gene) %>%
  dplyr::filter(!gene %in% genes_exclude) %>%
  dplyr::filter(`n_TFmotif_1E-4` > 0) %>%
  dplyr::filter(reQTL_pval < 0.000001)


#### focus on disease gene SNPs that have allele specificity ####
df <- bigtable.krt %>% 
  dplyr::filter(gene %in% deg$gene) %>%
  dplyr::filter(!gene %in% genes_exclude) %>%
  dplyr::filter(AlleleSpecificMarkCount > 0) %>%
  dplyr::filter(reQTL_pval < 0.00001 | PBSeQTL_pval < 0.00001 | IFNeQTL_pval < 0.00001 )

### any SNPs that might disrupt gene body?


#### any SNPs that might disrupt TF motif binding? ####
df <- bigtable.krt %>%
  dplyr::filter(`n_TFmotif_1E-4_ENCODETF` > 0) %>%
  dplyr::filter(gene %in% deg$gene) %>%
  dplyr::filter(reQTL_pval < 0.0001 | PBSeQTL_pval < 0.0001 | IFNeQTL_pval < 0.0001 )

#### SNP:gene pairs filtered by beta-comparison ####
betacomp.krt <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/KRT/masteroutput_with_colnames.txt") %>%
  mutate(tag=paste0(SNP,"_",GENE)) %>%
  dplyr::filter(p.betaComp10KPermut<0.01)

temp <- bigtable.krt %>% mutate(tag=paste0(QTL,"_",gene))

genes_exclude = c("ERAP2","HLA-DRB5","ITGA1","GBP3")

df <- left_join(betacomp.krt, temp , by="tag") %>%
  dplyr::filter(!gene %in% genes_exclude) %>%
  dplyr::filter(p.betaComp10KPermut<0.005) %>%
  dplyr::filter(!is.na(gene))

gene_list <- unique(df$gene)
write(df$gene,file="~/Downloads/temp_genes.txt")

df2 <- data.frame(snp=character(), gene=character(), p=numeric())
for (i in 1:length(gene_list)) {
  this_df <- df[which(df$gene==gene_list[i]),] %>% arrange(p.betaCompPnorm,p.betaComp10KPermut)
  # Only if this_df is not empty, proceed to rbind
  if(nrow(this_df) > 0) {
    # Create a new dataframe with the correct column names and the first row of sorted this_df
    new_row <- data.frame(snp=this_df[1,]$SNP, 
                          gene=this_df[1,]$GENE, 
                          p=this_df[1,]$p.betaComp10KPermut,
                          stringsAsFactors = FALSE)
    
    # Bind the new_row dataframe to df2
    df2 <- rbind(df2, new_row)
  }
}
df = df2


#### SNPs (lm betaComp handpicking) affecting TFs and downstream genes ####
f <- list.files("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/plot_keep/")
gene_lst <- unique( gsub(".png","",gsub("pairedPlot_3donorColored_CPM_rs.*_","",f)) )
write(gene_lst,"~/Downloads/nl/human/skin/eQTLs/chromatin/ANNOVAR/QTLs_1E-02/temp_genes.txt")
tf <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/literature/Lambert_2018_human_TFs.txt",header=F) %>%
  dplyr::filter(V1 %in% gene_lst) %>% pull(V1) %>% unique()
df <- bigtable.krt %>% 
  dplyr::filter(gene %in% tf)

#### TFmotif same as eGene ####
df <-  data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/temp_filtered_KRT.txt",header=F)
colnames(df) <- colnames(bigtable.krt)
# checked the plots. don't look good.

#### SNPs that have allele specificity, disrupt TF, interesting gene
f <- list.files("~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/plot_keep/")
gene_lst <- unique( gsub(".png","",gsub("pairedPlot_3donorColored_CPM_rs.*_","",f)) )
df <- bigtable.krt %>% 
  dplyr::filter(gene %in% gene_lst) %>%
  dplyr::filter(AlleleSpecificMarkCount>0) %>%
  dplyr::filter(num_of_peak_overlapping_TF>0) %>%
  dplyr::filter(`n_TFmotif_1E-4_ENCODETF`>0)
#### SNPs that are in high LD to vitiligo GWAS SNPs ####
snps <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results/vitiligo_highLD_SNPs.txt",header=F)$V1
df <- bigtable.krt %>% 
  dplyr::filter(QTL %in% snps) %>%
  dplyr::filter(p.betaComp10KPermut<0.01) %>%
  dplyr::filter(gene=="DHX58")
#### make plots ####
snp="rs2910686" # rs1874417, rs10864012, rs370174977,
gene="DEFB1"
make_reQTL_plot_CPM_3cts(snp,gene)


df <- bigtable.krt %>% 
  dplyr::filter(gene=="ZNF226") %>% 
  dplyr::filter(AlleleSpecificMarkCount>0 | num_of_peak_overlapping_TF>0 | `n_TFmotif_1E-4_ENCODETF`>0)

# check which genes are correlated to the TF
genes <- rownames(zlogCPM.heatmap)
this_tf <- "ETV7"
cor_result <- c()
for (i in 1:length(genes)) {
  this_gene <- genes[i]
  this_gene_CPM <- CPM_rna[this_gene,] %>% dplyr::select(contains("KRT")) %>% as.numeric()
  TF_CPM <- CPM_rna[this_tf,] %>% dplyr::select(contains("KRT")) %>% as.numeric()
  cor_result<- c(cor_result, cor(this_gene_CPM, TF_CPM, method="spearman"))
}
cor_result_df <- data.frame(gene=genes, spearmanCor=cor_result) %>% dplyr::filter(abs(spearmanCor)>0.85)
df <- data.frame(gene=cor_result_df$gene, QTL="rs704960")
#### make plots ####
g <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/literature/KEGG/GSEA_hallmark-genes-upregulated-in-response-to-IFNG.txt",header=F)$V1
df <- data.frame(QTL="rs7208907",gene=c("DHX58"))
this_outDir <- "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/plot_temp"
for (i in 1:nrow(df)) {
  this_snp <- as.character(df[i,"QTL"])
  this_gene <- as.character(df[i,"gene"])
  tryCatch({
    make_all_the_plots(this_snp, this_gene, this_outDir)
  }, error = function(e) {
    cat("An error occurred with SNP:", this_snp, "Gene:", this_gene, "\nError message:", e$message, "\n")
  })
  print(paste("[",as.character(i),"]",this_snp,this_gene,Sys.time()))
}


#### investigate SNPs ####
df1 <- bigtable.krt[which(bigtable.krt$gene=="DHX58"),] 
df2 <- df1[,c("QTL","reQTL_pval","p.betaComp10KPermut","p.betaCompPnorm","start")]
df2$distance <- df2$start - 42101404
df3 <- bigtable.krt[which(bigtable.krt$gene=="DHX58"),] %>%
  dplyr::filter(reQTL_pval<0.001 | p.betaComp10KPermut<0.001)
df2$threshold <- ifelse(df2$QTL %in% df3$QTL,"yes","no")

ggplot(df2, aes(x=distance, y=-log10(p.)))

#### investigate rs7503405:SOCS3 ####
this.snp="rs7503405"
this.gene="SOCS3"
CPM.pbs=CPM.pbs.krt
CPM.ifn=CPM.ifn.krt
# pick the snp
this.snp.genotype <- genotype %>% dplyr::filter(ID==this.snp)
snp.ref=this.snp.genotype$REF
snp.alt=this.snp.genotype$ALT

this.genotype=this.snp.genotype %>% 
  dplyr::select(-c("CHROM","START","END","ID","REF","ALT")) %>% t() %>% as.data.frame() %>%
  rownames_to_column("donor") %>% set_colnames(c("donor","genotype"))

this.genotype$donor=this.genotype$donor %>% gsub(".GT","", . )
this.genotype$genotype = this.genotype$genotype %>%
  gsub("0/0",paste0(snp.ref,snp.ref), . ) %>%
  gsub("0/1",paste0(snp.ref,snp.alt), . ) %>%
  gsub("1/1",paste0(snp.alt,snp.alt), . ) %>%
  gsub("./.",NA, . )
this.genotype$genotype <- ordered(this.genotype$genotype, 
                                  levels <- c(paste0(snp.ref,snp.ref),
                                              paste0(snp.ref,snp.alt),
                                              paste0(snp.alt,snp.alt)))
this.genotype = this.genotype %>% na.omit()
rm(this.snp.genotype)

# pick the gene
this.CPM.PBS <- CPM.pbs[this.gene,] %>% t() %>% as.data.frame() %>%
  rownames_to_column("donor") %>% set_colnames(c("donor","PBS"))
this.CPM.IFN <- CPM.ifn[this.gene,] %>% t() %>% as.data.frame() %>%
  rownames_to_column("donor") %>% set_colnames(c("donor","IFN"))
this.CPM <- inner_join(this.CPM.PBS, this.CPM.IFN, by="donor") %>% 
  inner_join( . , this.genotype, by="donor")

# make heatmap with genotype in columns
load("~/Downloads/nl/human/skin/eQTLs/integrative_analysis/heatmap_RNA_round2/myEnvironment_heatmap_RNAseq_DEgenes_14kmm-hc_04182023_padj0.05_log2FC1.5_avgCPM10.RData")

# 4. load genes of interest
g <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/literature/KEGG/GSEA_hallmark-genes-upregulated-in-response-to-IFNG.txt",header=F)$V1
zlogCPM.heatmap.new.v1 <- km_res.rna.new[,1:ncol(CPM_rna_DEG)]
temp <- colnames(zlogCPM.heatmap.new.v1) %>% dplyr::mutate()
zlogCPM.heatmap_with_gene <- zlogCPM.heatmap.new %>% dplyr::filter(rownames(.) %in% g)
row_anno_index <- which(rownames(zlogCPM.heatmap.new) %in% rownames(zlogCPM.heatmap_with_gene))
row_anno_label <- rownames(zlogCPM.heatmap.new)[row_anno_index]

# 5. plot heatmap
column_ha = HeatmapAnnotation(
  condition = colnames(zlogCPM.heatmap.new) %>% 
    lapply( . , function(x) unlist(strsplit(x , "_", ))[2]) %>% unlist(),
  celltype = colnames(zlogCPM.heatmap.new) %>% 
    lapply( . , function(x) unlist(strsplit(x , "_", ))[3]) %>% unlist(),
  col = list(condition = c("PBS" = "black", "IFN" = "#FA8072"),
             celltype = c("KRT" = "#DAF7A6", "MEL" = "#964B00", "FRB" = "#CCCCFF")),
  show_legend = c(TRUE,TRUE))

row_ha = rowAnnotation(
  foo = anno_mark(at = row_anno_index, 
                  labels = row_anno_label,
                  which="row",
                  labels_gp = gpar(fontsize=10),
                  padding=unit(0.5,"mm")))

pdf(paste0(Dir,"/heatmap_RNA_round2/heatmap_reordered_RNAseq_DEgenes_14kmm-hc_04182023_padj0.05_log2FC1.5_avgCPM10.pdf"),width=24,height=30)
heatmap.rna <- draw(Heatmap(zlogCPM.heatmap.new, 
                            name = "z-score log2CPM",
                            top_annotation = column_ha,
                            right_annotation = row_ha,
                            column_title = "DE genes in 3 celltypes, control vs IGNg stimulated, known IFNg inducible genes labeled",
                            cluster_columns = FALSE,
                            cluster_rows = FALSE,
                            column_split = c(rep(1,length(colnames(CPM_rna_DEG)[grepl("PBS_KRT",colnames(CPM_rna_DEG))])), 
                                             rep(2,length(colnames(CPM_rna_DEG)[grepl("PBS_MEL",colnames(CPM_rna_DEG))])), 
                                             rep(3,length(colnames(CPM_rna_DEG)[grepl("PBS_FRB",colnames(CPM_rna_DEG))])), 
                                             rep(4,length(colnames(CPM_rna_DEG)[grepl("IFN_KRT",colnames(CPM_rna_DEG))])), 
                                             rep(5,length(colnames(CPM_rna_DEG)[grepl("IFN_MEL",colnames(CPM_rna_DEG))])),
                                             rep(6,length(colnames(CPM_rna_DEG)[grepl("IFN_FRB",colnames(CPM_rna_DEG))]))),
                            row_split = split_instructions.new, 
                            show_row_names = FALSE,
                            row_title_rot = 0,
                            column_names_rot = 45,
                            row_gap = unit(0, "mm"),
                            border = TRUE,
                            column_names_gp = gpar(fontsize = 10),
                            column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                            row_title_gp = gpar(fontsize = 10, fontface = "bold")))
dev.off()


#### count number of candidate reQTLs per cell type ####
df <- bigtable.krt %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut)))
write(unique(df$QTL), "~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/KRT_candidate_reQTLs_052024.txt")

df <- bigtable.frb %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut)))
write(unique(df$QTL), "~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/FRB_candidate_reQTLs_052024.txt")

df <- bigtable.mel %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut)))
write(unique(df$QTL), "~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/MEL_candidate_reQTLs_052024.txt")

#### investigate rs2304206 ####
snps <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/all_candidate_reQTLs_052024_also_skinDis_melanoma_autoImmuneDis_snps.txt",header=F)$V1
df1 <- bigtable.mel %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut))) %>%
  dplyr::filter(QTL %in% snps)%>%
  dplyr::select(c(QTL,gene))
df1$celltype="MEL"

df2 <- bigtable.krt %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut))) %>%
  dplyr::filter(QTL %in% snps)%>%
  dplyr::select(c(QTL,gene))
df2$celltype="KRT"

df3 <- bigtable.frb %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut))) %>%
  dplyr::filter(QTL %in% snps) %>%
  dplyr::select(c(QTL,gene))
df3$celltype="FRB"

df <- rbind(df1,rbind(df2,df3)) %>%
  group_by(QTL,gene) %>%
  summarize(celltype=paste(celltype, collapse=","))
data.table::fwrite(df, file="~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/all_candidate_reQTL_gene_pairs_052024_also_skinDis_melanoma_autoImmuneDis_snps.txt", append=T, sep="\t",quote=F)

gwas_snp <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/lite-GWAS-catalog-all-associations-autosomal-snp.tsv")
gwas_snp <- gwas_snp %>% 
  tidyr::separate( . , col=`STRONGEST_SNP-RISK_ALLELE`, into=c("SNP","risk_allele"),sep="-") %>%
  group_by(SNP,risk_allele) %>%
  summarize(`DISEASE/TRAIT`=paste(`DISEASE/TRAIT`, collapse=","))


df4 <- left_join(df, gwas_snp, by=c("QTL"="SNP"))


df5 <- df4[,c("QTL","gene")] %>% distinct()

this_outDir <- "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/plot_temp"
for (i in 38:nrow(df5)) {
  this_snp <- as.character(df5[i,"QTL"])
  this_gene <- as.character(df5[i,"gene"])
  tryCatch({
    make_all_the_plots(this_snp, this_gene, this_outDir)
  }, error = function(e) {
    cat("An error occurred with SNP:", this_snp, "Gene:", this_gene, "\nError message:", e$message, "\n")
  })
  print(paste("[",as.character(i),"]",this_snp,this_gene,Sys.time()))
}

df6 <- df4 %>% dplyr::filter(QTL %in% c("rs11742570","rs11150589","rs7428430","rs11574938","rs1830610"))


#### investigate high LD tags of GWAS SNPs ####
df.krt <- bigtable.krt %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter(!(is.na(reQTL_pval) & is.na(p.betaComp10KPermut))) %>%
  #dplyr::filter(reQTL_pval<0.00001 & !is.na(reQTL_pval)) %>%
  dplyr::filter(!is.na(GWAS_Trait)) %>%
  mutate(celltype="KRT") %>%
  dplyr::select(-c(cRE_overlap_KRT,cRE_dynamic_KRT,model_zscore,model_zscore_pval))
df.krt <- df.krt[grep("vitiligo|lupus|psoriasis|derma|autoimmune|melanoma",df.krt$GWAS_Trait,ignore.case = TRUE),]

df.frb <- bigtable.frb %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter(!(is.na(reQTL_pval) & is.na(p.betaComp10KPermut))) %>%
  #dplyr::filter(reQTL_pval<0.00001 & !is.na(reQTL_pval)) %>%
  dplyr::filter(!is.na(GWAS_Trait)) %>%
  mutate(celltype="FRB") %>%
  dplyr::select(-c(cRE_overlap_FRB,cRE_dynamic_FRB,model_zscore,model_zscore_pval))
df.frb <- df.frb[grep("vitiligo|lupus|psoriasis|derma|autoimmune|melanoma",df.frb$GWAS_Trait,ignore.case = TRUE),]

df.mel <- bigtable.mel %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter(!(is.na(reQTL_pval) & is.na(p.betaComp10KPermut))) %>%
  #dplyr::filter(reQTL_pval<0.00001 & !is.na(reQTL_pval)) %>%
  dplyr::filter(!is.na(GWAS_Trait)) %>%
  mutate(celltype="MEL") %>%
  dplyr::select(-c(cRE_overlap_MEL,cRE_dynamic_MEL))
df.mel <- df.mel[grep("vitiligo|lupus|psoriasis|derma|autoimmune|melanoma",df.mel$GWAS_Trait,ignore.case = TRUE),]

df <- rbind(df.krt, rbind(df.frb, df.mel)) %>%
  dplyr::filter(!gene %in% c("ERAP2","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5")) %>%
  group_by(gene) %>%        # Group the data by the gene column
  filter(p.betaComp10KPermut == min(p.betaComp10KPermut)) 
this_outDir <- "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/plot_temp"
for (i in 1:nrow(df)) {
  this_snp <- as.character(df[i,"QTL"])
  this_gene <- as.character(df[i,"gene"])
  tryCatch({
    make_all_the_plots(this_snp, this_gene, this_outDir)
  }, error = function(e) {
    cat("An error occurred with SNP:", this_snp, "Gene:", this_gene, "\nError message:", e$message, "\n")
  })
  print(paste("[",as.character(i),"]",this_snp,this_gene,Sys.time()))
}

df1 <- rbind(df.krt, rbind(df.frb, df.mel))
df2 <- df1[grep("melanoma",df1$GWAS_Trait,ignore.case = TRUE),] %>%
  pull(QTL) %>% unique()
write(unique(df1$QTL),"~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs/old_reQTL_linked_to_GWAS.txt")

#### MEL-specific reQTLs ####
df.krt <- bigtable.krt %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut))) %>%
  mutate(celltype="KRT") %>%
  dplyr::select(-c(cRE_overlap_KRT,cRE_dynamic_KRT,model_zscore,model_zscore_pval))

df.frb <- bigtable.frb %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut))) %>%
  mutate(celltype="FRB") %>%
  dplyr::select(-c(cRE_overlap_FRB,cRE_dynamic_FRB,model_zscore,model_zscore_pval))

df.mel <- bigtable.mel %>%
  dplyr::filter(reQTL_pval<0.00001 | p.betaComp10KPermut<0.001) %>%
  dplyr::filter((!is.na(reQTL_pval)) & (!is.na(p.betaComp10KPermut))) %>%
  mutate(celltype="MEL") %>%
  dplyr::select(-c(cRE_overlap_MEL,cRE_dynamic_MEL))

df <- rbind(df.krt, rbind(df.frb, df.mel)) %>%
  group_by(QTL,gene,REF,ALT) %>%
  summarize(celltype=paste(celltype, collapse=","))

df1 <- df %>% dplyr::filter(celltype=="MEL")

df2 <- rbind(df.krt, rbind(df.frb, df.mel)) %>%
  dplyr::filter(QTL %in% df1$QTL) %>%
  arrange(reQTL_pval,p.betaComp10KPermut)

this_outDir <- "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/plot_temp"
for (i in 5988:5988) {
  this_snp <- as.character(df[i,"QTL"])
  this_gene <- as.character(df[i,"gene"])
  tryCatch({
    make_all_the_plots(this_snp, this_gene, this_outDir)
  }, error = function(e) {
    cat("An error occurred with SNP:", this_snp, "Gene:", this_gene, "\nError message:", e$message, "\n")
  })
  print(paste("[",as.character(i),"]",this_snp,this_gene,Sys.time()))
}

make_all_the_plots("rs8081327", "CXCL11", this_outDir)
