#!/usr/bin/env Rscript
# written by Crystal Shan 07/2024
# load all QTLs and compile LDproxy table

library(tidyverse)
library(magrittr)
library(LDlinkR)

# read in QTL table (separately by celltype)
dir <- "~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results"
outDir <- paste0(dir,"/LDproxy")

inDir <- "~/Downloads/nl/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/"
celltype <- "KRT"
modelF <- data.table::fread(paste0(inDir,"compiled_table_",celltype,"_reQTL.txt"))
snpList <- unique(c(modelF$QTL))
# exclude reQTLs that were also in MEL (for faster processing)
reQTL_MEL <- data.table::fread(paste0(inDir,"compiled_table_","MEL","_reQTL.txt")) %>% pull(QTL) %>% unique()
snpList <- modelF %>% 
  dplyr::filter(!QTL %in% reQTL_MEL) %>% 
  pull(QTL) %>% unique()


# function
compile_LDproxyTable <- function(this.snp, threshold) {
  tryCatch({
  my_token = "0f622f6ac228" #cornell token
  this_res <- LDproxy(snp = this.snp,
                      pop = "ALL",
                      r2d = "r2",
                      token = my_token,
                      genome_build = "grch38")
  this_res <- this_res %>% 
    dplyr::filter(R2>=threshold) %>% 
    arrange(desc(R2)) %>% 
    dplyr::mutate(query_snp=this.snp)
  return(this_res)
  }, 
  error = function(err) {
  this_res <- data.frame(result="")
   return(this_res)
  })
}

# loop through the table and get their LDtrait table
for (i in 5376:length(snpList)) {
  start_time <- Sys.time()
  this_snp <- snpList[i]
  this_table <- compile_LDproxyTable(this_snp, 0.6)
  data.table::fwrite(this_table, file=paste0(outDir,"/LDproxy_results_",celltype,"_reQTLs.txt"), 
                     quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  time_log <- data.frame(Iteration = i, ExecutionTime = execution_time, EndTime = end_time)
  data.table::fwrite(time_log, file=paste0(outDir,"/log_LDproxy_results_",celltype,"_reQTLs.txt"), 
                     quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
}

# after the table is fully compiled. add these colnames to the tables
result_colnames <- c("RS_Number","Coord","Alleles","MAF","Distance","Dprime","R2","Correlated_Alleles","FOREGEdb","RegulomeDB","Function","query_snp")
log_colnames <- c("Iteration","ExecutionTime","EndTime")