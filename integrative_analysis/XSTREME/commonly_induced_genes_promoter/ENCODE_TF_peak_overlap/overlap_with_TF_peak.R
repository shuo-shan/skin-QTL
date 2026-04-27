library(tidyverse)
library(dplyr)
library(magrittr)

Dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis"
link_gene_promoter = readRDS(paste0(Dir,"/heatmap_RNA_round2/link_gene_promoter_within_300kbp_of_TSS_allcts.Rdata"))
link_gene_enhancer = readRDS(paste0(Dir,"/heatmap_RNA_round2/link_gene_enhancer_within_300kbp_of_TSS_allcts.Rdata"))
bed <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/chromatin/RegulatoryElements/method3/regulatory_regions.bed")
colnames(bed) <- c("chr","start","end","region")
load(paste0(Dir,"/heatmap_RNA_round2/myEnvironment_heatmap_RNAseq_DEgenes_14kmm-hc_04182023_padj0.05_log2FC1.5_avgCPM10.RData"))
outDir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/XSTREME/commonly_induced_genes_promoter/ENCODE_TF_peak_overlap"

# this chunk of code is taken from this:
# ~/Downloads/nl/human/skin/eQTLs/integrative_analysis/XSTREME/commonly_induced_genes_promoter/fetch_promoters.Rmd
# get gene list
gene_list <- c()
for (i in c(1,2,3)){
  this_cluster_genes = rownames(km_res.rna.new)[which(km_res.rna.new$newclass==i)]
  gene_list <- c(gene_list, this_cluster_genes)
}

# get promoter list
promoter_df <- link_gene_promoter %>% dplyr::filter(gene %in% gene_list)
# only pick promoters that exist in all three cell types
promoter_list <- table( promoter_df[,c("promoter")]) %>% as.data.frame() %>% dplyr::filter(Freq==3) %>% pull(promoter) %>% as.character()
# 187 out of 263 genes have at least 1 promoter that's active in all three celltypes
# update promoter list based on new info
pruned_promoter_df <- promoter_df[which(promoter_df$promoter %in% promoter_list),]

#### load TF ChIPseq overlap result
TF_overlap <- data.table::fread(paste0(outDir,"/promoter_overlapping_TF_peaks.bed"))
colnames(TF_overlap) <- c("chr","start","end","promoter","gene","entry")
tab <- left_join(pruned_promoter_df, TF_overlap, by="promoter") %>% na.omit()

#### check which TF are linked to which gene
gene_TF <- tab %>% dplyr::select(c(gene,TF)) %>% distinct()
gene_TF_n <- gene_TF %>% group_by(TF) %>% mutate(count = n()) %>% dplyr::select(c(TF,count)) %>% distinct()
gene_TF_summary <- gene_TF %>% group_by(TF) %>% summarize(gene=paste(gene,collapse=",")) %>%
  left_join(gene_TF_n, . , by="TF")

data.table::fwrite(gene_TF_summary, file=paste0(outDir,"/ENCODE_TF_peak_intersecting_promoters_of_commonly_induced_genes_summary.txt"),
                   quote=F, sep="\t", row.names=F, col.names=T)


#### regulatory network construction
library(igraph)
library(plotly)
edges <- data.frame(TF=gene_TF$TF, Genes=gene_TF$gene)
g <- igraph::graph_from_data_frame(edges, directed = TRUE)

ig <- plot_ly() %>% 
  add_trace(data = as_data_frame(g, what = "edges"), 
            x = ~from_id, y = ~to_id, type = 'scatter', mode = 'lines') %>% 
  add_trace(data = as_data_frame(g, what = "vertices"), 
            x = ~name, y = ~name, type = 'scatter', mode = 'markers', text = ~name, marker = list(size = 10))


library(networkD3)
# Sample data
src <- c("A", "A", "A", "A", "B", "C", "D", "D")
target <- c("B", "C", "D", "J", "E", "F", "G", "H")
df <- gene_TF[grep("IRF",gene_TF$gene),]
networkData <- data.frame(df$TF, df$gene)

# Create a simple network graph
sn <- simpleNetwork(networkData)
# Custom CSS to change font size
sn <- htmlwidgets::onRender(
  sn,
  '
  function(el, x) {
    // Reduce font size
    d3.select(el).selectAll(".node text").style("font-size", "8px");

    // Center text inside nodes by adjusting x and y attributes
    d3.select(el).selectAll(".node text")
      .attr("x", 0)    // Centers horizontally; adjust as needed
      .attr("y", 3)    // Vertically center; small offset helps align better visually
      .style("text-anchor", "middle"); // Ensures text is centered on coordinates
  }
  '
)
print(sn)





