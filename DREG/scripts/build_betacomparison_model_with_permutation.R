#!/usr/bin/env Rscript
# written by Crystal Shan 03/2022. modified on 12/2022.
# do lm() with linear regression separately for PBS and IFN, then compare beta and s.e between the two models.
library(tidyverse)
library(caret)
library(glmnet)

# parsing input arguments
args = commandArgs(trailingOnly=TRUE)
Dir=args[1]
snp=args[2]
g=args[3]
PBSFile=args[4]
IFNFile=args[5]
genotypeFile=args[6]

# # # # #example of snp-gene pair
# Dir="~/Downloads/nl/human/skin/eQTLs/DREG/reQTL_betacomparison/KRT"
# snp="rs1036173"
# g="ADRB2"
# PBSFile=paste0(Dir,"/","modelingResult","/temp_",snp,"/PBS_",g,".txt")
# IFNFile=paste0(Dir,"/","modelingResult","/temp_",snp,"/IFN_",g,".txt")
# genotypeFile=paste0(Dir,"/","modelingResult","/temp_",snp,"/genotype.txt")

##########_________main script________________________________________________
##### compile data table ##### 
PBS=data.table::fread(PBSFile, header=TRUE)
IFN=data.table::fread(IFNFile, header=TRUE)
genotype=data.table::fread(genotypeFile, header=TRUE)
cat("compiled data matrix for modeling, building models now \n")

##### model 8. beta-comparison model. CPM. no covariates for both PBS and IFN. #####
model8.PBS <- function() {
  phenotype=PBS
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, genotype, by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    model.min <- lm(phenotype ~ ., as.data.frame(mat))
    p.nominal <- summary(model.min)$coefficients["genotype", "Pr(>|t|)"]
    p.adj <- NA
    beta <- summary(model.min)$coefficients["genotype", "Estimate"]
    stdErr <- summary(model.min)$coefficients["genotype", "Std. Error"]
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("PBSeQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("PBSeQTL: successfully built the model...")
    p.nominal <- summary(model.min)$coefficients["genotype", "Pr(>|t|)"]
    p.adj <- NA
    beta <- summary(model.min)$coefficients["genotype", "Estimate"]
    stdErr <- summary(model.min)$coefficients["genotype", "Std. Error"]
    this_result_table=c(p.nominal, p.adj, beta, stdErr)
  }
  return(this_result_table)
}
model8.IFN <- function() {
  phenotype=IFN
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, genotype, by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    model.min <- lm(phenotype ~ ., as.data.frame(mat))
    p.nominal <- summary(model.min)$coefficients["genotype", "Pr(>|t|)"]
    p.adj <- NA
    beta <- summary(model.min)$coefficients["genotype", "Estimate"]
    stdErr <- summary(model.min)$coefficients["genotype", "Std. Error"]
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("IFNeQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("IFNeQTL: successfully built the model...")
    p.nominal <- summary(model.min)$coefficients["genotype", "Pr(>|t|)"]
    p.adj <- NA
    beta <- summary(model.min)$coefficients["genotype", "Estimate"]
    stdErr <- summary(model.min)$coefficients["genotype", "Std. Error"]
    this_result_table=c(p.nominal, p.adj, beta, stdErr)
  }
  return(this_result_table)
}
model8.PBSandIFN.permute <- function() {
  phenotype=PBS
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, genotype, by="donor") %>% 
    dplyr::select(-donor)
  
  # permute donor genotype without replacement
  mat.shuffled <- mat
  mat.shuffled$genotype <- sample(mat$genotype, nrow(mat), replace = FALSE)
  # fit permuted matrix
  model.permutated <- lm(phenotype ~ ., mat.shuffled)
  pval_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Pr(>|t|)"]
  beta_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Estimate"]
  se_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Std. Error"]
  
  this_result_table.pbs = c(pval_genotype.permutated, beta_genotype.permutated, se_genotype.permutated)
  
  ###### now shuffle with IFN
  phenotype=IFN
  colnames(phenotype)=c("donor","phenotype")
  mat.ifn <- inner_join(phenotype, genotype, by="donor") %>% 
    dplyr::select(-donor)
  
  # permute donor genotype without replacement
  mat.shuffled.ifn <- mat.ifn[,c(1,2)]
  mat.shuffled.ifn$genotype <- mat.shuffled$genotype
  # fit permuted matrix
  model.permutated <- lm(phenotype ~ ., mat.shuffled.ifn)
  pval_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Pr(>|t|)"]
  beta_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Estimate"]
  se_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Std. Error"]
  
  this_result_table.ifn = c(pval_genotype.permutated, beta_genotype.permutated, se_genotype.permutated)
  
  this_result_table <- c(this_result_table.pbs, this_result_table.ifn)
  return(this_result_table)
}

# compute beta comparison ratio
results.PBS <- model8.PBS()
results.IFN <- model8.IFN()
beta_comparison <- (results.PBS[3] - results.IFN[3]) / (sqrt((results.PBS[4]^2) + (results.IFN[4]^2)))
p.pnorm <- 2 * pmin(pnorm(beta_comparison), 1 - pnorm(beta_comparison))
# permutation
cat("permuting for 10K times now... \n")
stats_from_permutation <- data.frame(pval.pbs=numeric(), beta.pbs=numeric(), se.pbs=numeric(), 
                                     pval.ifn=numeric(), beta.ifn=numeric(), se.ifn=numeric())
for (i in 1:10000) {
  stats_from_permutation[i,1:6] <- model8.PBSandIFN.permute()
}
stats_from_permutation$beta_comparison <- (stats_from_permutation$beta.pbs - stats_from_permutation$beta.ifn) / sqrt((stats_from_permutation$se.pbs)^2 + (stats_from_permutation$se.ifn)^2)
cat("done with permutation \n")
# calculate empirical p.val
num_of_beta_comparison_more_extreme = length( stats_from_permutation$beta_comparison[ which(abs(stats_from_permutation$beta_comparison) >= abs(beta_comparison))])
p.empirical <- num_of_beta_comparison_more_extreme / nrow(stats_from_permutation)

##### Write output table. ##### 
result_table=c(snp,g,
               beta_comparison,
               p.pnorm,
               p.empirical)
result_table[is.na(result_table)] <- "."

# # column names for result_table
# colnames(result_table) <- c("SNP","GENE","z.betaComp","p.betaCompPnorm","p.betaComp10KPermut")

output_fname=paste0(Dir,"/","modelingResult","/",snp,"_",g,".txt")

write.table(paste(result_table, collapse="\t"),file=output_fname,quote=F, sep="\t", row.names=F, col.names=F)
cat("Done! Cya!!!! =D \n")