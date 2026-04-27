library(dplyr)
library(magrittr)
library(ggplot2)
# load data
dir="~/Downloads/nl/human/skin/eQTLs/donor_genotyping/hg38/072024_donors"
dict <- data.table::fread(paste0(dir,"/file_donor_lookup.txt"), header = F) %>%
  set_colnames(c("file","donorID"))
dict$fileID <- gsub(".vcf.gz","",dict$file)
dict$donorID <- gsub("UMMSGarberSkineQTL","",dict$donorID)


# raw coverage overview 
coverage_table <- data.table::fread(paste0(dir,"/merged_qc_data.tsv")) %>%
  dplyr::filter(quality_control_type.key %in% c("raw_coverage","effective_coverage_min"))
coverage_table$fileID <- gsub("_qc.json","",coverage_table$file)
coverage_table <- left_join(coverage_table[,c("fileID","quality_control_type.key","quality_control.value_measured")], dict,
                            by="fileID")
coverage_table <- coverage_table %>%
  arrange(if_else(quality_control_type.key == "raw_coverage", 
                  quality_control.value_measured, 
                  NA_real_))

# effective_coverage_min
p <- ggplot(coverage_table, aes(x=donorID, y=quality_control.value_measured, 
                                group=quality_control_type.key, 
                                color=quality_control_type.key)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(dir,"/coverage_summary.png"), plot = p)
