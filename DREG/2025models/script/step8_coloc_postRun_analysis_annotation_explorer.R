# explore annotation table
library(tidyverse)

dir="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/coloc/summary/best/"
gene="CTSS"
snp="rs3768018"

df=fread(paste0(dir,gene,"/annotation_big.txt"))

# which SNP has TF motif overlap?
this_df <- df %>% 
  dplyr::filter(`n_TFmotif_1E-4`>0) %>%
  dplyr::select(QTL,
                `n_TFmotif_1E-8`, `TFmotif_1E-8`,
                peak_overlapping_TF,
                AlleleSpecificMark,
                cRE_name,  cRE_dynamic_FRB,
                Func.refGene, Gene.refGene, dist_to_tss,
                `n_TFmotif_1E-6`, `TFmotif_1E-6`,
                `n_TFmotif_1E-4_ENCODETF`, `TFmotif_1E-4_ENCODETF`,
                `n_TFmotif_1E-4`, `TFmotif_1E-4`) %>%
  arrange(desc(`n_TFmotif_1E-8`), desc(`n_TFmotif_1E-6`), desc(`n_TFmotif_1E-4`))
View(this_df)


