library(dplyr)
library(ggplot2)
library(ggrepel)

het <- data.table::fread("~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/plinkQC/genodata/data.het")
het$donor <- het$FID %>%
  gsub("skineQTL-","", .)


muF  <- mean(het$F, na.rm = TRUE)
sdF  <- sd(het$F,  na.rm = TRUE)
low  <- muF - 3 * sdF
high <- muF + 3 * sdF
het$fail_het <- het$F < low | het$F > high

ggplot(het, aes(x = donor, y = F, color = fail_het, label = donor)) +
  geom_point() +
  geom_hline(yintercept = c(low, muF, high),
             linetype = c("dashed","solid","dashed"),
             color = c("grey60","grey30","grey60")) +
  geom_text_repel(size = 2.5, max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(x = "Sample (donor ID)", y = "F (heterozygosity)",
       color = "Fail heterozygosity?",
       title = "Per-individual heterozygosity (F)") +
  theme_bw()
