#!/usr/bin/env Rscript
# written by Crystal Shan 08/2024
# This script joins the modeling result table for DHX58 KRT reQTL candidate rs8081327 vs all KRT DEG
# and explores the possibility of it being a trans-QTL for IFN responding genes (esp. type I IFN)

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)

dir="~/Downloads/nl/human/skin/eQTLs/case_studies/DHX58/rs8081327/"
##### 1. load modeling tables and compile ##### 
## beta comparison modeling
tab.betacomp <- data.table::fread(paste0(dir,"/model_betacomparison/masteroutput_with_colnames.txt"), header=T)
## log2FC linear modeling
tab.linear <- data.table::fread(paste0(dir,"/model_log2FClinear/masteroutput_with_colnames.txt")) 
## join tables
tab.joined <- merge(tab.betacomp, tab.linear, by=c("SNP","GENE"), all = TRUE)
tab.joined[tab.joined == "."] <- NA
tab.joined <- tab.joined %>% dplyr::filter(!is.na(z.betaComp) | !is.na(tab.joined$log2FC_featureSelected_reQTL_genotype_beta))

##### 2. annotate genes by type I, type II, or type I or II from Carol's data ##### 
IFNB_genes <- data.table::fread(paste0(dir,"/Carol_clusters/genes_cluster_4.txt"), header=F) %>%
  separate_rows(V2, sep = "/") %>% 
  dplyr::filter(V2 %in% tab.joined$GENE) %>%
  pull(V2)

##### 3. Fisher's Exact test for whether significant hits are
View(tab.joined[which(tab.joined$p.betaComp10KPermut < 0.01),])
tab.joined$CPM_minimalModel_IFNeQTL_genotype_pval <- as.numeric(tab.joined$CPM_minimalModel_IFNeQTL_genotype_pval)
View(tab.joined[tab.joined$CPM_minimalModel_IFNeQTL_genotype_pval < 0.05])
# nevermind. the pvalue is so bad for the top antiviral genes.

##### 4. GSEA of IFNB responding genes in pannus keratinocytes based on the beta comparison model p.permut value
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)

# convert gene symbols to Entrez IDs
gene_symbols <- tab.joined$GENE
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# prepare a named vector of p-values
df <- tab.joined %>% dplyr::filter(GENE %in% gene_entrez$SYMBOL) %>% dplyr::select(c("GENE","p.betaComp10KPermut")) %>%
  left_join( . , gene_entrez, by=c("GENE"="SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

p_values <- df$p.betaComp10KPermut
names(p_values) <- df$ENTREZID
p_values <- -log10(p_values)

# Create a gene list ranked by -log10(p-value)
gene_list <- sort(p_values, decreasing = TRUE)

# Create a custom gene set list
gene_set_list <- list(
  IFNB_genes = gene_entrez$ENTREZID[gene_entrez$SYMBOL %in% IFNB_genes]
)
# Perform GSEA
gsea_result <- fgsea(pathways = gene_set_list,
                     stats = p_values,
                     scoreType = "pos",
                     nperm = 1000)
gsea_result

##### 5. GSEA of all KEGG pathways, based on p.permut of IFNG reQTL in KRT in DHX58 reQTL across all genes.
# Use clusterProfiler to get KEGG pathways
library(msigdbr)
library(enrichplot)
library(ggupset)

# convert gene symbols to Entrez IDs
gene_symbols <- tab.joined$GENE
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# prepare a named vector of p-values
df <- tab.joined %>% dplyr::filter(GENE %in% gene_entrez$SYMBOL) %>% dplyr::select(c("GENE","p.betaComp10KPermut")) %>%
  left_join( . , gene_entrez, by=c("GENE"="SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

p_values <- df$p.betaComp10KPermut
names(p_values) <- df$ENTREZID
p_values <- -log10(p_values)


selected_genes <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = df$GENE,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL") %>% na.omit()

gse <- enrichKEGG(gene=selected_genes$ENTREZID,
                  pvalueCutoff = 0.05)
gse.sig <- gse@result[which(gse@result$p.adjust<0.05),]
selected_pathways <- gse$Description[1:nrow(gse.sig)]

# bubble plot
p.dot <- dotplot(gse, showCategory=selected_pathways, font.size=10) + 
  ggtitle("KEGG pathway enrichment for gene set:")
pdf(file=paste0(dir,"/GSEA_KEGGpathways_KRT_rs8081327_vs_allGenes_p.permut_bubbleplot.pdf"), width=8,height=14)
p.dot
dev.off()

# bar plot
p.bar <- gse %>% dplyr::mutate( . , qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore", showCategory=selected_pathways, font.size=10) +
  xlab("-log10(p.adj)")
pdf(file=paste0(dir,"/GSEA_KEGGpathways_KRT_rs8081327_vs_allGenes_p.permut_barplot.pdf"), width=8,height=12)
p.bar
dev.off()
# upset plot
p.upset <- upsetplot(gse, n=nrow(gse.sig))
pdf(file=paste0(dir,"/GSEA_KEGGpathways_KRT_rs8081327_vs_allGenes_p.permut_upsetplot.pdf"), width=8,height=12)
p.upset
dev.off()

### make plot for each pathway
# Select a pathway of interest, say "hsa04621" (e.g., NOD-like receptor signaling pathway)
pathway_of_interest <- "hsa04621"

# Extract the genes for that pathway
pathway_genes <- setReadable(gse, 'org.Hs.eg.db', keyType = 'ENTREZID')
selected_pathway <- subset(pathway_genes, pathway_genes@result$ID == pathway_of_interest)

# Create a plot with the genes on x-axis and pvalue on the y-axis
# First, extract the p-values and genes
genes_in_pathway <- unlist(strsplit(selected_pathway$geneID, "/"))
pvalue <- selected_pathway$p.adjust
pathway_name <- selected_pathway$Description

# Plot using ggplot2 for visualization
ggplot(selected_pathway, aes(x = reorder(geneID, p.adjust), y = -log10(p.adjust))) +
  geom_point(aes(size = Count), color = "blue") +
  coord_flip() +
  labs(title = pathway_name,
       x = "Genes",
       y = "-log10(p.adjust)",
       size = "Gene Count") +
  theme_minimal()

# Plot enrichment for the pathway of interest
plotEnrichment(pathway = pathway_of_interest,
               stats = p_values,
               title = paste("Enrichment Plot for", pathway_of_interest))

gseaplot(gse, geneSetID = selected_pathway$ID)




gse <- gseKEGG(gene=selected_genes$ENTREZID,
                  pvalueCutoff = 0.05)
gse.sig <- gse@result[which(gse@result$p.adjust<0.05),]
selected_pathways <- gse$Description[1:nrow(gse.sig)]



##### new code
# load the MSigDB gene sets into a named list
pathways.kegg <- gmtPathways("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/GSEA/MSigDB/c2.cp.kegg.v7.5.1.symbols.gmt")
pathways.GO_BP <- gmtPathways("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/GSEA/MSigDB/c5.go.bp.v7.5.1.symbols.gmt")
pathways <- pathways.kegg
MSigDBtype="KEGG" # type of MSigDB gene set

# convert gene symbols to Entrez IDs
gene_symbols <- tab.joined$GENE
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# prepare a named vector of p-values
df <- tab.joined %>% dplyr::filter(GENE %in% gene_entrez$SYMBOL) %>% dplyr::select(c("GENE","p.betaComp10KPermut")) %>%
  left_join( . , gene_entrez, by=c("GENE"="SYMBOL")) %>%
  distinct(ENTREZID, .keep_all = TRUE)

p_values <- df$p.betaComp10KPermut
names(p_values) <- df$GENE
p_values <- -log10(p_values)

# Create a gene list ranked by -log10(p-value)
gene_list <- sort(p_values, decreasing = TRUE)


# Perform GSEA
gsea_result <- fgseaMultilevel(pathways = pathways,
                     stats = gene_list,
                     scoreType = "pos",
                     minSize = 10,
                     nPermSimple = 10000)
View(gsea_result)
hist(sapply(pathways, length),breaks=100)

fgseaResTidy <- gsea_result %>% as_tibble %>%
  arrange(desc(NES))


pathway.name="KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY"
this.padj=fgseaResTidy[which(fgseaResTidy$pathway==pathway.name),"padj"]
this.NES=fgseaResTidy[which(fgseaResTidy$pathway==pathway.name),"NES"]
leadingEdge=fgseaResTidy[which(fgseaResTidy$pathway==pathway.name),"leadingEdge"] %>% unlist() %>% as.data.frame()
data.table::fwrite(leadingEdge,paste0(dir,"/GSEA_KEGGpathway_KRT_rs8081327_vs_allGenes_p.permut_in_",pathway.name,".txt"),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

# Enrichment plot for a given pathway/gene-set
p <- plotEnrichment(pathway=pathways[[pathway.name]], p_values) +
  labs(title=pathway.name,
       subtitle=paste0("padj = ",round(this.padj,3),", NES = ",round(this.NES,3))) +
  theme(plot.title=element_text(size=8, face="bold", color="black")) +
  theme(plot.subtitle=element_text(size=8, face="italic", color="black"))
p
#pdf(paste0(Dir,"/gsea_enrichment_plot_",pathway.name,"_TGFB1_DEG.pdf"))
#print(p)
#dev.off()

# create a named vector [ranked genes]
ranks <- f$log2FoldChange
names(ranks) <- f$SYMBOL

# run fgsea algorithm:
# NES: enrichment score normalized to mean enrichment of random samples of the same size
fgseaRes <- fgseaMultilevel(pathways=pathways, stats=ranks)
#View(fgseaRes[which(fgseaRes$padj<0.05),])
res.sig <- fgseaRes %>% arrange(padj) %>% as.matrix()

#write.table(res.sig,file="~/Downloads/nl/human/skin/melanocyte_chemokine/RNAseq_analysis/TGFB/result_TGFB2_DEGenes_GSEA_pathwayKEGG.txt",,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
#View(head(fgseaRes[order(padj), ], n=10))
#View(head(fgseaRes[order(padj, abs(NES)), ], n=10))

# tidy the results:
fgseaResTidy <- fgseaRes %>% as_tibble %>%
  arrange(desc(NES)) # order by normalized enrichement score

# to see what genes are in each of these pathways:
gene.in.pathway <- pathways %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols=c(SYMBOL)) %>%
  inner_join(f, by="SYMBOL")



###### archive code #####

## gene annotation 
gene_clusters <- data.table::fread(paste0(dir,"/Carol_clusters/gene_cluster_annotation.txt"), header=F) %>%
  set_colnames(c("GENE","cluster"))
gene_clusters$cluster <- gsub("cluster6","IFNB",gene_clusters$cluster)
gene_clusters$cluster <- gsub("cluster7","IFNG",gene_clusters$cluster)
gene_clusters$cluster <- gsub("cluster8","IFNG",gene_clusters$cluster)
gene_clusters$cluster <- gsub("cluster5","IFNBandIFNG",gene_clusters$cluster)
gene_clusters$cluster[grep("cluster*",gene_clusters$cluster)] <- "TNF"

##### 3. t-test? #####
tab.joined.anno <- left_join(tab.joined, gene_clusters, by="GENE")

df <- tab.joined.anno[,c("GENE","z.betaComp","p.betaCompPnorm",
                         "p.betaComp10KPermut","log2FC_featureSelected_reQTL_genotype_beta",
                         "rankNormCPM_minimalModel_IFNeQTL_genotype_beta",
                         "rankNormCPM_minimalModel_IFNeQTL_genotype_pval",
                         "cluster")]
df$cluster[is.na(df$cluster)] <- "otherDEG"
df$neglog10Ppermut <- -log10(df$p.betaComp10KPermut)

ggplot(df, aes(x=cluster, y=neglog10Ppermut, color=cluster)) +
  geom_violin() +
  geom_jitter(width=0.1)

ggplot(df, aes(x=cluster, y=z.betaComp, color=cluster)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.4)

df$log2FC_featureSelected_reQTL_genotype_beta <- as.numeric(df$log2FC_featureSelected_reQTL_genotype_beta)
df %>%
  na.omit() %>%
  ggplot(., aes(x=cluster, y=log2FC_featureSelected_reQTL_genotype_beta, color=cluster)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.4)

df$rankNormCPM_minimalModel_IFNeQTL_genotype_beta <- as.numeric(df$rankNormCPM_minimalModel_IFNeQTL_genotype_beta)
df %>%
  na.omit() %>%
  ggplot(., aes(x=cluster, y=rankNormCPM_minimalModel_IFNeQTL_genotype_beta, color=cluster)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.4)

##### 3. plot log2FC of genes based on the rs8081327 genotype #####
CPM.PBS <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_PBS.txt") %>% column_to_rownames("V1")
CPM.IFN <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt") %>% column_to_rownames("V1")

logCPM.pbs <- log2(CPM.PBS+1)
logCPM.ifn <- log2(CPM.IFN+1)
log2FC <- logCPM.ifn - logCPM.pbs

genotype <- data.table::fread(paste0(dir,"/model_log2FClinear/genotype.txt")) %>% as.data.frame()
DEG <- data.table::fread(paste0(dir,"/KRT_DEgenes.txt"),header=F)$V1
log2FC.AA <- log2FC[DEG,] %>% na.omit %>% dplyr::select(genotype[genotype$genotype==2,]$donor)
log2FC.notAA <- log2FC[DEG,] %>% na.omit %>% dplyr::select(genotype[genotype$genotype!=2,]$donor)

res <- data.frame(gene=as.character(), t.test.pval=as.numeric())
for (i in 1:length(rownames(log2FC.AA))) {
  this_gene <- rownames(log2FC.AA)[i]
  this_res <- t.test(as.numeric(log2FC.AA[this_gene,]),
                     as.numeric(log2FC.notAA[this_gene,]))
  res <- rbind(res, c(this_gene, this_res$p.value))
}
colnames(res) <- c("gene","t.test.pval")
res$t.test.pval <- as.numeric(res$t.test.pval)
res.annotated <- left_join(res, gene_clusters, by=c("gene"="GENE"))
res.annotated$cluster[is.na(res.annotated$cluster)] <- "otherDEG"
res.annotated$t.test.negLog10Pval <- -log10(res.annotated$t.test.pval)

ggplot(res.annotated, aes(x=cluster, y=t.test.negLog10Pval, color=cluster)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.4)


##### 4. plot CPM of IFNG condition of genes based on the rs8081327 genotype #####
CPM.IFN <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/CPM_expressedGenes_IFN.txt") %>% column_to_rownames("V1")

genotype <- data.table::fread(paste0(dir,"/model_log2FClinear/genotype.txt")) %>% as.data.frame()
DEG <- data.table::fread(paste0(dir,"/KRT_DEgenes.txt"),header=F)$V1
CPM.IFN.AA <- CPM.IFN[DEG,] %>% na.omit %>% dplyr::select(genotype[genotype$genotype==2,]$donor)
CPM.IFN.notAA <- CPM.IFN[DEG,] %>% na.omit %>% dplyr::select(genotype[genotype$genotype!=2,]$donor)

res.ifn <- data.frame(gene=as.character(), t.test.pval=as.numeric())
for (i in 1:length(rownames(CPM.IFN.AA))) {
  this_gene <- rownames(CPM.IFN.AA)[i]
  this_res <- t.test(as.numeric(CPM.IFN.AA[this_gene,]),
                     as.numeric(CPM.IFN.notAA[this_gene,]))
  res.ifn <- rbind(res.ifn, c(this_gene, this_res$p.value))
}
colnames(res.ifn) <- c("gene","t.test.pval")
res.ifn$t.test.pval <- as.numeric(res.ifn$t.test.pval)
res.ifn.annotated <- left_join(res.ifn, gene_clusters, by=c("gene"="GENE"))
res.ifn.annotated$cluster[is.na(res.ifn.annotated$cluster)] <- "otherDEG"
res.ifn.annotated$t.test.negLog10Pval <- -log10(res.ifn.annotated$t.test.pval)
res.ifn.annotated$t.test.padj <- p.adjust(res.ifn.annotated$t.test.pval)

ggplot(res.ifn.annotated, aes(x=cluster, y=t.test.negLog10Pval, color=cluster)) +
  geom_violin() +
  geom_jitter(width=0.1, alpha=0.4)
