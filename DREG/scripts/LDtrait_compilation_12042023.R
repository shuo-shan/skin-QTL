#!/usr/bin/env Rscript
# written by Crystal Shan 11/2023
# load all QTLs and compile LDtrait table

library(tidyverse)
library(magrittr)
library(LDlinkR)

# read in QTL table (total of 4884 QTLs)
dir <- "~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results"
outDir <- paste0(dir,"/LDtrait")
modelF <- paste0(dir,"/modeling_results_featureSelectedModel_phenotypeRankNormCPM.txt")
modelTab <- data.table::fread(modelF, fill=TRUE) %>% unique()
snpList <- unique(modelTab$SNP)

# loop through the table and get their LDtrait table
compile_LDtraitTable <- function(this.snp) {
  tryCatch({
  my_token = "908c3efbf915"
  this_res <- LDtrait(snps = this.snp,
                      pop = "ALL",
                      r2d = "r2",
                      r2d_threshold = 0.1,
                      win_size = 5e+05,
                      token = my_token,
                      genome_build = "grch38")
  this_res <- this_res %>% arrange(desc(R2)) 
  return(this_res)
  }, 
  error = function(err) {
  this_res <- data.frame(result="no associated trait was found at r2 > 0.1 in ALL population by LDtrait")
   return(this_res)
  })
}

time_log <- data.frame(Iteration = numeric(),
                       ExecutionTime = numeric(),
                       EndTime = numeric())

for (i in 450:length(snpList)) {
  start_time <- Sys.time()
  this_snp <- snpList[i]
  this_table <- compile_LDtraitTable(this_snp)
  data.table::fwrite(this_table, file=paste0(outDir,"/",this_snp,".txt"), 
                     quote=F, sep="\t", row.names=F, col.names=T)
  end_time <- Sys.time()
  
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  time_log <- rbind(time_log, data.frame(Iteration = i, ExecutionTime = execution_time, EndTime = end_time))
}


# organize the compiled table
resultF <- paste0(dir,"/QTLs_LDtrait_association_ALLpopulation_r20.1.txt")
resultTab <- data.table::fread(resultF, fill=TRUE) %>% unique()

traits <- resultTab %>% dplyr::filter(resultTab$R2 > 0.8) %>% pull(GWAS_Trait)
View(table(traits))
traitList <- c("Systemic lupus erythematosus",
               "Systemic lupus erythematosus (MTAG)",
                unique(traits[grepl("MHC class I",traits)]))
resultTab_key <- resultTab %>% dplyr::filter(GWAS_Trait %in% traitList)
unique(resultTab_key$Query)
GWAS_hits <- unique(resultTab_key$RS_Number)
resultTab_key2 <- resultTab_key %>% dplyr::filter(Query %in% GWAS_hits)
# "rs1131476"  "rs11755393" "rs12620999" "rs3869132"  "rs4681679"  "rs9270984"  "rs9271366"  "rs9469857" 