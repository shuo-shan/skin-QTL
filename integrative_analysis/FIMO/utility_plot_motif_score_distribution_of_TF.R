library(dplyr)
library(magrittr)
library(tidyverse)

Dir="~/Downloads/nl/human/skin/eQTLs/integrative_analysis/FIMO"
setwd(Dir)

TF <- "STAT1"
dir.pos <- paste0(TF,"_",TF,"HighCorrelationGenesPromoters")
dir.neg <- paste0(TF,"_",TF,"LowCorrelationGenesPromoters")

fileName <- paste0("fimo_",TF,".bed")


file.pos <- read.table(paste0(Dir,"/",dir.pos,"/",fileName)) %>%
  set_colnames(c("chr","start","end","motif_start","motif_end","strand","score","pval")) %>%
  mutate(tag="HighCorrelationGenesPromoters") %>%
  dplyr::filter(pval <= 1)

file.neg <- read.table(paste0(Dir,"/",dir.neg,"/",fileName)) %>%
  set_colnames(c("chr","start","end","motif_start","motif_end","strand","score","pval")) %>%
  mutate(tag="LowCorrelationGenesPromoters") %>%
  dplyr::filter(pval <= 1)

df <- rbind(file.pos[,c("score","tag")],
            file.neg[,c("score","tag")])

# histogram plot
ggplot(df, aes(x=score, fill=tag)) +
  geom_histogram(color="white",alpha=0.6, position="identity", binwidth = 1) +
  xlab("FIMO motif matching score") +
  ggtitle(paste0(TF," motif matching score distribution, all")) +
  theme_classic() +
  theme(legend.position="top")


# density plot
ggplot(df, aes(x=score, fill=tag)) +
  geom_density(alpha=0.6) +
  xlab("FIMO motif matching score") +
  ggtitle(paste0(TF," motif matching score distribution density, all")) +
  theme_classic() +
  theme(legend.position="top")









