dir="~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/yuqing_disease_DEG"
file_list <- list.files(dir)

dds.dm <- data.table::fread(paste0(dir,"/",file_list[1]))
dds.lupus <- data.table::fread(paste0(dir,"/",file_list[2]))
dds.psoriasis <- data.table::fread(paste0(dir,"/",file_list[3]))
dds.vit <- data.table::fread(paste0(dir,"/",file_list[4]))

DEG.dm <- dds.dm %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(V1) %>% mutate(disease="DM") %>% set_colnames(c("gene","disease"))
DEG.lup <- dds.lupus %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(V1) %>% mutate(disease="lupus") %>% set_colnames(c("gene","disease"))
DEG.pso <- dds.psoriasis %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(V1) %>% mutate(disease="pso") %>% set_colnames(c("gene","disease"))
DEG.vit <- dds.vit %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(V1) %>% mutate(disease="vit") %>% set_colnames(c("gene","disease"))

DEG.1 <- rbind(DEG.dm, DEG.lup, DEG.pso, DEG.vit) %>% 
  group_by(gene) %>%
  mutate(count = n()) %>% 
  select(-disease) %>%
  distinct()

DEG.2 <- rbind(DEG.dm, DEG.lup, DEG.pso, DEG.vit) %>% 
  group_by(gene) %>%
  summarize(disease = paste(disease, collapse=",")) %>%
  distinct()

DEG <- left_join(DEG.1, DEG.2, by="gene")

data.table::fwrite(DEG, file=paste0(dir,"/four_diseases_DEG_padj0.05_log2FC1.txt"), quote=F, col.names = F, sep="\t")
