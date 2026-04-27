# written by shuo.shan@umassmed.edu 03/2024
library(tidyverse)
library(magrittr)

# load QTL and TF CHIPseq intersect result
dir="~/Downloads/nl/human/skin/eQTLs/chromatin/TF_motif_overlapping/QTL_1E-02"
# lenient q-value 1E-04 but also the TF from ENCODE TF ChIPseq dataset
f1 = data.table::fread(paste0(dir,"/fimo_output_SNP-containing_qval1E-4_ENCODETF.bed"))
colnames(f1) <- c("chr","start","end","ID","REF","ALT","chr.TF","start.TF","end.TF","ID2","strand","TF_motif","score","pval","qval","sequence")
f1$TF <- gsub(".H12CORE.*", "", f1$TF_motif)
data = f1[,c("ID","TF")] %>% distinct()
data.collapsed.n = data %>% 
  group_by(ID) %>% 
  mutate(count = n()) %>% 
  distinct( . , ID, .keep_all=TRUE) %>%
  select(-TF)
data.collapsed = data %>% 
  arrange(TF) %>%
  group_by(ID) %>% 
  summarize(TF=paste(TF,collapse=','))
data.joined <- left_join(data.collapsed.n , data.collapsed, by="ID")
colnames(data.joined) <- c("ID","n_TFmotif_1E-4_ENCODETF","TFmotif_1E-4_ENCODETF")
f1.result <- data.joined
rm(f1,data.collapsed,data.collapsed.n,data.joined,data)

# lenient q-value 1E-04 but any TF from the HOCOMOCOv12 core database
f2 = data.table::fread(paste0(dir,"/fimo_output_SNP-containing_qval1E-4.bed"))
colnames(f2) <- c("chr","start","end","ID","REF","ALT","chr.TF","start.TF","end.TF","ID2","strand","TF_motif","score","pval","qval","sequence")
f2$TF <- gsub(".H12CORE.*", "", f2$TF_motif)
data = f2[,c("ID","TF")] %>% distinct()
data.collapsed.n = data %>% 
  group_by(ID) %>% 
  mutate(count = n()) %>% 
  distinct( . , ID, .keep_all=TRUE) %>%
  select(-TF)
data.collapsed = data %>% 
  arrange(TF) %>%
  group_by(ID) %>% 
  summarize(TF=paste(TF,collapse=','))
data.joined <- left_join(data.collapsed.n , data.collapsed, by="ID")
colnames(data.joined) <- c("ID","n_TFmotif_1E-4","TFmotif_1E-4")
f2.result <- data.joined
rm(f2,data.collapsed,data.collapsed.n,data.joined,data)

# moderate-strictness q-value 1E-06 but any TF from the HOCOMOCOv12 core database
f3 = data.table::fread(paste0(dir,"/fimo_output_SNP-containing_qval1E-6.bed"))
colnames(f3) <- c("chr","start","end","ID","REF","ALT","chr.TF","start.TF","end.TF","ID2","strand","TF_motif","score","pval","qval","sequence")
f3$TF <- gsub(".H12CORE.*", "", f3$TF_motif)
data = f3[,c("ID","TF")] %>% distinct()
data.collapsed.n = data %>% 
  group_by(ID) %>% 
  mutate(count = n()) %>% 
  distinct( . , ID, .keep_all=TRUE) %>%
  select(-TF)
data.collapsed = data %>% 
  arrange(TF) %>%
  group_by(ID) %>% 
  summarize(TF=paste(TF,collapse=','))
data.joined <- left_join(data.collapsed.n , data.collapsed, by="ID")
colnames(data.joined) <- c("ID","n_TFmotif_1E-6","TFmotif_1E-6")
f3.result <- data.joined
rm(f3,data.collapsed,data.collapsed.n,data.joined,data)

# strict q-value 1E-08 but any TF from the HOCOMOCOv12 core database
f4 = data.table::fread(paste0(dir,"/fimo_output_SNP-containing_qval1E-8.bed"))
colnames(f4) <- c("chr","start","end","ID","REF","ALT","chr.TF","start.TF","end.TF","ID2","strand","TF_motif","score","pval","qval","sequence")
f4$TF <- gsub(".H12CORE.*", "", f4$TF_motif)
data = f4[,c("ID","TF")] %>% distinct()
data.collapsed.n = data %>% 
  group_by(ID) %>% 
  mutate(count = n()) %>% 
  distinct( . , ID, .keep_all=TRUE) %>%
  select(-TF)
data.collapsed = data %>% 
  arrange(TF) %>%
  group_by(ID) %>% 
  summarize(TF=paste(TF,collapse=','))
data.joined <- left_join(data.collapsed.n , data.collapsed, by="ID")
colnames(data.joined) <- c("ID","n_TFmotif_1E-8","TFmotif_1E-8")
f4.result <- data.joined
rm(f4,data.collapsed,data.collapsed.n,data.joined,data)

# compile output
bedf=data.table::fread(paste0(dir,"/QTL.bed"))
colnames(bedf) <- c("chr","start","end","ID","REF","ALT")
res1 <- left_join(bedf, f1.result, by="ID")
res1$`n_TFmotif_1E-4_ENCODETF`[is.na(res1$`n_TFmotif_1E-4_ENCODETF`)] <- 0

res2 <- left_join(res1, f2.result, by="ID")
res2$`n_TFmotif_1E-4`[is.na(res2$`n_TFmotif_1E-4`)] <- 0

res3 <- left_join(res2, f3.result, by="ID")
res3$`n_TFmotif_1E-6`[is.na(res3$`n_TFmotif_1E-6`)] <- 0

res4 <- left_join(res3, f4.result, by="ID")
res4$`n_TFmotif_1E-8`[is.na(res4$`n_TFmotif_1E-8`)] <- 0

final_result <- res4

# write to file
data.table::fwrite(final_result, file=paste0(dir,"/QTL_1E-02_overlapping_TFmotif.bed"),
                   quote = F, sep="\t", row.names = F, col.names = T)
