#!/usr/bin/env Rscript
# written by Crystal Shan 03/2022. modified on 12/2022.
# do lm() with linear regression separately for PBS and IFN, then compare genotype between the two models.
library(tidyverse)
library(caret)
library(glmnet)

# parsing input arguments
args = commandArgs(trailingOnly=TRUE)
Dir=args[1]
snp=args[2]
g=args[3]
combination=args[4] # a=n_genotypePCs, b=n_latentVariables
PBSFile=args[5]
IFNFile=args[6]
genotypeFile=args[7]
covariatesFile=args[8]

##example of snp-gene pair
#dir="~/Downloads/nl/human/skin/eQTLs/DREG/edQTL_also_eQTL/"
#snp="rs11750025"
#g="ERAP2"
#combination="a1b1" # a=n_genotypePCs, b=n_latentVariables
#PBSFile=paste0(dir,"/","modelingResult_",combination,"/temp_",snp,"/PBS_",g,".txt")
#IFNFile=paste0(dir,"/","modelingResult_",combination,"/temp_",snp,"/IFN_",g,".txt")
#genotypeFile=paste0(dir,"/","modelingResult_",combination,"/temp_",snp,"/genotype.txt")
#covariatesFile=paste0(dir,"/covariates_",combination,".txt")

##########_________main script________________________________________________
### compile data table
PBS=data.table::fread(PBSFile, header=TRUE)
IFN=data.table::fread(IFNFile, header=TRUE)
genotype=data.table::fread(genotypeFile, header=TRUE)
covariates=data.table::fread(covariatesFile, header=TRUE)
cat("compiled data matrix for modeling, building models now \n")


### build linear regression model with regularization
build_linear_model_log2FC_nopermutation <- function() {
  fudge=10
  phenotype=inner_join(PBS,IFN,by="donor") %>% mutate(phenotype=log2((IFN+fudge) / (PBS+fudge))) %>% dplyr::select(-c(PBS,IFN))
  mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    cv.glmnet(mat[,-1], mat[,1], alpha = 1)},
    error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("failed to build to model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("successfully built the model \n")
    # Find the best lambda using cross-validation
    cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1)
    # Fit the final model on the training data
    model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
    # Display regression coefficients
    modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
    colnames(modelcoef) <- c("feature","coef")
    modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
    selected_features <- modelcoef$feature
    # If genotype is kept as a selected feature, fit the model
    if ('genotype' %in% selected_features){
      # Build the linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
      p.adj <- NA
      all_features <- data.frame(feature=colnames(mat)[2:length(colnames(mat))] )
      x <- as.data.frame(summary(model.fs)$coefficients)
      x$feature <- rownames(x)
      p_list <- data.frame(pval = x[-1,"Pr(>|t|)"])
      p_list$feature <- x[-1,"feature"] %>% gsub("`","",.)
      all_features <- left_join(all_features, p_list, by="feature")
      this_result_table=c(p.nominal,p.adj,beta,stdErr,all_features$pval)
    } else {
      cat("genotype did not pass feature selection! \n")
      this_result_table=c(NA,NA,NA,NA,rep(NA,ncol(mat)-1))}
  }
  return(this_result_table)
}
build_linear_model_PBS_nopermutation <- function() {
  phenotype=PBS
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    cv.glmnet(mat[,-1], mat[,1], alpha = 1)},
    error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("failed to build to model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("successfully built the model \n")
    # Find the best lambda using cross-validation
    cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5)
    # Fit the final model on the training data
    model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
    # Display regression coefficients
    modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
    colnames(modelcoef) <- c("feature","coef")
    modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
    selected_features <- modelcoef$feature
    # If genotype is kept as a selected feature, fit the model
    if ('genotype' %in% selected_features){
      # Build the linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
      p.adj <- NA
      all_features <- data.frame(feature=colnames(mat)[2:length(colnames(mat))] )
      x <- as.data.frame(summary(model.fs)$coefficients)
      x$feature <- rownames(x)
      p_list <- data.frame(pval = x[-1,"Pr(>|t|)"])
      p_list$feature <- x[-1,"feature"] %>% gsub("`","",.)
      all_features <- left_join(all_features, p_list, by="feature")
      this_result_table=c(p.nominal,p.adj,beta,stdErr,all_features$pval)
    } else {
      cat("genotype did not pass feature selection! \n")
      this_result_table=c(NA,NA,NA,NA,rep(NA,ncol(mat)-1))
      }
  }
  return(this_result_table)
}
build_linear_model_IFN_nopermutation <- function() {
  phenotype=IFN
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    cv.glmnet(mat[,-1], mat[,1], alpha = 1)},
    error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("failed to build to model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("successfully built the model \n")
    # Find the best lambda using cross-validation
    cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1)
    # Fit the final model on the training data
    model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
    # Display regression coefficients
    modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
    colnames(modelcoef) <- c("feature","coef")
    modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
    selected_features <- modelcoef$feature
    # If genotype is kept as a selected feature, fit the model
    if ('genotype' %in% selected_features){
      # Build the linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
      p.adj <- NA
      all_features <- data.frame(feature=colnames(mat)[2:length(colnames(mat))] )
      x <- as.data.frame(summary(model.fs)$coefficients)
      x$feature <- rownames(x)
      p_list <- data.frame(pval = x[-1,"Pr(>|t|)"])
      p_list$feature <- x[-1,"feature"] %>% gsub("`","",.)
      all_features <- left_join(all_features, p_list, by="feature")
      this_result_table=c(p.nominal,p.adj,beta,stdErr,all_features$pval)
    } else {
      cat("genotype did not pass feature selection! \n")
      this_result_table=c(NA,NA,NA,NA,rep(NA,ncol(mat)-1))
      }
  }
  return(this_result_table)
}


### Write output table.
res_log2FC=build_linear_model_log2FC_nopermutation()
res_PBS=build_linear_model_PBS_nopermutation()
res_IFN=build_linear_model_IFN_nopermutation()

result_table=c(snp,g,res_log2FC,res_PBS,res_IFN)

output_fname=paste0(Dir,"/","modelingResult","_",combination,"/",snp,"_",g,".txt")

write.table(paste(result_table, collapse="\t"),file=output_fname,quote=F, sep="\t", row.names=F, col.names=F)
cat("Done! Cya!!!! =D \n")



###### unused code
# build_linear_model_log2FC_permutation <- function(fudge) {
#   fudge=100
#   phenotype=inner_join(PBS,IFN,by="donor") %>% mutate(phenotype=log2((IFN+fudge) / (PBS+fudge))) %>% dplyr::select(-c(PBS,IFN))
#   mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor)
#   # first try if model building can be successful before proceeding
#   flag <- TRUE
#   tryCatch({
#     lm(phenotype ~ ., mat)},
#     error = function(e){flag<<-FALSE})
#   if (!flag) stop("failed to build model")
#   cat("successfully built the model \n")
#   # build model
#   model <- lm(phenotype ~ ., mat)
#   summary(model)
#   p.nominal <- summary(model)$coefficients["genotype", "Pr(>|t|)"]
#   beta <- summary(model)$coefficients["genotype", "Estimate"]
#   stdErr <- summary(model)$coefficients["genotype", "Std. Error"]
#   res <- data.frame(residual=summary(model)$residuals) %>% t()
#   ### permutation without resplacement on genotype 1000 times and perform linear regression each time.
#   if (p.nominal < 0.01) {
#     cat("permuting response variables... \n")
#     pval_from_permutation <- c()
#     for (i in 1:1000) {
#       # permute donor indices without replacement
#       shuffled_donors <- sample(rownames(mat), nrow(mat), replace = FALSE) %>% as.numeric
#       # permute data matrix genotype based on the shuffled donor indices
#       mat.permutated <- mat %>% mutate(genotype=mat$genotype[shuffled_donors])
#       # fit permuted matrix
#       model.permutated <- lm(phenotype ~ ., mat.permutated)
#       pval_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Pr(>|t|)"]
#       # write to variable
#       pval_from_permutation <- c(pval_from_permutation, pval_genotype.permutated)
#     }
#     # calculated p.adj
#     num_of_pval_more_extreme = length( pval_from_permutation[ which(pval_from_permutation <= p.nominal)])
#     p.adj <- num_of_pval_more_extreme / length(pval_from_permutation)
#     cat("Permutation testing completed. \n")
#   } else {
#     p.adj <- NA
#   }
#   this_result_table=as.character(c(p.nominal,p.adj,beta,stdErr))
#   return(this_result_table)
# }
