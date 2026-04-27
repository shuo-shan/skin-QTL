#!/usr/bin/env Rscript
# written by Crystal Shan 09/2023
library(tidyverse)
library(LDlinkR)

args = commandArgs(trailingOnly=TRUE)
this.snp=args[1]
outDir=args[2]


compile_LDtraitTable <- function(this.snp) {
  #tryCatch({
    my_token = "908c3efbf915"
    this.res <- LDtrait(snps = this.snp,
                        pop = "CEU",
                        r2d = "r2",
                        r2d_threshold = 0.6,
                        win_size = 5e+05,
                        token = my_token,
                        genome_build = "grch38")
    this.res <- this.res %>% arrange(desc(R2)) 
    return(this.res)
  #}, 
  #error = function(err) {
    #this.res <- data.frame(result="no associated trait was found")
   # return(this.res)
  #})
}

this.table <- compile_LDtraitTable(this.snp)

write.table(this.table, file=paste0(outDir,"/LDtrait_result_",this.snp,".txt"), quote=F, sep="\t", row.names = F)
