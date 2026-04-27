library(tidyverse)
library(magrittr)
library(ggpubr)

dir="~/Downloads/nl/human/skin/eQTLs/GWAS_SNPs"

# intergenic file
inF="gwas_catalog_v1.0_all_associations_intergenic.txt"
tab=data.table::fread(paste0(dir,"/",inF),header=T)
tab.expand = tab %>% separate_rows( . , MAPPED_GENE, sep=" - ")
tab.expand[tab.expand == ''] <- NA
tab.expand[tab.expand == 'NA'] <- NA
tab.expand = na.omit(tab.expand)
tab.expand$SNP_ID_CURRENT = paste0("rs",tab.expand$SNP_ID_CURRENT)
tab.expand = tab.expand %>% separate(., col=STRONGEST_SNP_RISK_ALLELE, into=c("SNP","RISK_ALLELE"),sep="-")
tab.expand = na.omit(tab.expand)
tab.expand = tab.expand %>% dplyr::filter(RISK_ALLELE %in% c("A","C","G","T"))
tab.expand = tab.expand %>% dplyr::select(-c("SNP","SNPS"))
tab.expand$GENE_POS = "intergenic"
tab.intergenic = tab.expand
rm(tab,tab.expand)

# intragenic file
inF="gwas_catalog_v1.0_all_associations_intragenic.txt"
tab=data.table::fread(paste0(dir,"/",inF),header=T)
tab=tab[,1:8]
tab.expand = tab %>% separate_rows( . , MAPPED_GENE, sep=", ") %>%
  mutate(temp = paste0(MAPPED_GENE,":",SNPS)) %>%
  distinct(temp, .keep_all = TRUE) %>%
  dplyr::select(-temp)
tab.expand[tab.expand == ''] <- NA
tab.expand[tab.expand == 'NA'] <- NA
tab.expand = na.omit(tab.expand)
tab.expand$SNP_ID_CURRENT = paste0("rs",tab.expand$SNP_ID_CURRENT)
tab.expand = tab.expand %>% separate(., col=STRONGEST_SNP_RISK_ALLELE, into=c("SNP","RISK_ALLELE"),sep="-")
tab.expand = na.omit(tab.expand)
tab.expand = tab.expand %>% dplyr::filter(RISK_ALLELE %in% c("A","C","G","T"))
tab.expand = tab.expand %>% dplyr::select(-c("SNP","SNPS"))
tab.expand$GENE_POS = "intragenic"
tab.intragenic = tab.expand
rm(tab,tab.expand)

# joining intergenic and intragenic files
outtab = rbind(tab.intergenic, tab.intragenic)

# write to file
data.table::fwrite(outtab, file = paste0(dir,"/gwas_catalog_v1.0_all_associations_snps.txt"),
                   quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)



