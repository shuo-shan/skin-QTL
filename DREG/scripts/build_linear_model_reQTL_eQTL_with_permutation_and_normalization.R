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

# #example of snp-gene pair
# dir="~/Downloads/nl/human/skin/eQTLs/DREG/KRT_new"
# snp="rs2648726"
# g="ITGA1"
# combination="a10b10" # a=n_genotypePCs, b=n_latentVariables
# PBSFile=paste0(dir,"/","modelingResult_",combination,"/temp_",snp,"/PBS_",g,".txt")
# IFNFile=paste0(dir,"/","modelingResult_",combination,"/temp_",snp,"/IFN_",g,".txt")
# genotypeFile=paste0(dir,"/","modelingResult_",combination,"/temp_",snp,"/genotype.txt")
# covariatesFile=paste0(dir,"/covariates_",combination,".txt")

##########_________main script________________________________________________
### compile data table
PBS=data.table::fread(PBSFile, header=TRUE)
IFN=data.table::fread(IFNFile, header=TRUE)
genotype=data.table::fread(genotypeFile, header=TRUE)
covariates=data.table::fread(covariatesFile, header=TRUE)
cat("compiled data matrix for modeling, building models now \n")


### inverse normal transformation
# INT stands for rank-based inverse normal transformation
# function
rankNorm <- function (y) {
  # input y: numeric vector of CPM across all genes per donor.
  k <- 0.375 # an offset to ensure the z-score is finite. from Blom transform.
  n <- length(y)
  # Ranks.
  r <- rank(y, ties.method = "average") # if same value, same rank
  # Apply transformation.
  r.prob <- (r - k) / (n - 2 * k + 1)
  y.rankNorm <- stats::qnorm(r.prob) # because qnorm(1) is infinite and qnorm(-1) is -Inf
  return(y.rankNorm)
}
PBS.norm <- data.frame(donor=PBS$donor, PBS=rankNorm(PBS$PBS))
IFN.norm <- data.frame(donor=IFN$donor, IFN=rankNorm(IFN$IFN))
# We inferred a significant eQTL when the nominal P-value was less than 5 x 10-8 , 
# which is a threshold commonly applied to genome-wide association studies (45-49), 
# and corresponded to a false discovery rate (50) of4% and 12% for the ANOVA and additive model, respectively
# Vicente CT, Revez JA, Ferreira MAR. Lessons from ten years of genome-wide association studies of asthma. Clin Transl Immunology.2017;6(12):e165.

### build linear regression model with lasso regularization (alpha=1) and permutation
# output table header for each function: genotype feature information on p.nominal, p.permutation, beta, stdError
build_linear_model_log2FCnorm_permutation <- function() {
  fudge=10
  phenotype=inner_join(PBS,IFN,by="donor") %>% mutate(phenotype=log2((IFN+fudge) / (PBS+fudge))) %>% dplyr::select(-c(PBS,IFN))
  # rankNorm transform phenotype
  phenotype <- data.frame(donor=phenotype$donor, phenotype=rankNorm(phenotype$phenotype))
  mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
    cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
    # Fit the final model on the data
    model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
    # Display regression coefficients
    modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
    colnames(modelcoef) <- c("feature","coef")
    modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
    selected_features <- modelcoef$feature
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("failed to build to model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("successfully built the model \n")
    # If genotype is kept as a selected feature, fit the model
    if ('genotype' %in% selected_features) {
      # Build the linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
      p.adj <- NA
      pval_from_permutation <- c()
      for (i in 1:1000) {
        # permute donor indices without replacement
        shuffled_donors <- sample(rownames(mat.fs), nrow(mat.fs), replace = FALSE) %>% as.numeric
        # permute data matrix genotype based on the shuffled donor indices
        mat.permutated <- mat.fs %>% mutate(genotype=mat.fs$genotype[shuffled_donors])
        # fit permuted matrix
        model.permutated <- lm(phenotype ~ ., mat.permutated)
        pval_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Pr(>|t|)"]
        # write to variable
        pval_from_permutation <- c(pval_from_permutation, pval_genotype.permutated)
      }
      # calculated p.adj
      num_of_pval_more_extreme = length( pval_from_permutation[ which(pval_from_permutation <= p.nominal)])
      p.adj <- num_of_pval_more_extreme / length(pval_from_permutation)
      cat("Permutation testing completed. \n")
      
      this_result_table=as.character(c(p.nominal,p.adj,beta,stdErr))
      return(this_result_table)
      # header: log2FCnorm_pnominal, log2FCnorm_pperm, log2FCnorm_beta, log2FCnorm_se
    }
  }
}
build_linear_model_PBSnorm_permutation <- function() {
  phenotype=PBS.norm
  # rankNorm transform phenotype
  phenotype <- data.frame(donor=phenotype$donor, phenotype=rankNorm(phenotype$PBS))
  mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
    cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
    # Fit the final model on the data
    model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
    # Display regression coefficients
    modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
    colnames(modelcoef) <- c("feature","coef")
    modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
    selected_features <- modelcoef$feature
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("failed to build to model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("successfully built the model \n")
    # If genotype is kept as a selected feature, fit the model
    if ('genotype' %in% selected_features) {
      # Build the linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
      p.adj <- NA
      pval_from_permutation <- c()
      for (i in 1:1000) {
        # permute donor indices without replacement
        shuffled_donors <- sample(rownames(mat.fs), nrow(mat.fs), replace = FALSE) %>% as.numeric
        # permute data matrix genotype based on the shuffled donor indices
        mat.permutated <- mat.fs %>% mutate(genotype=mat.fs$genotype[shuffled_donors])
        # fit permuted matrix
        model.permutated <- lm(phenotype ~ ., mat.permutated)
        pval_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Pr(>|t|)"]
        # write to variable
        pval_from_permutation <- c(pval_from_permutation, pval_genotype.permutated)
      }
      # calculated p.adj
      num_of_pval_more_extreme = length( pval_from_permutation[ which(pval_from_permutation <= p.nominal)])
      p.adj <- num_of_pval_more_extreme / length(pval_from_permutation)
      cat("Permutation testing completed. \n")
      
      this_result_table=as.character(c(p.nominal,p.adj,beta,stdErr))
      return(this_result_table)
      # header: PBSnorm_pnominal, PBSnorm_pperm, PBSnorm_beta, PBSnorm_se
    }
    else {
      this_result_table=c(NA,NA,NA,NA)
      return(this_result_table)
    }
  }
}
build_linear_model_IFNnorm_permutation <- function() {
  phenotype=IFN.norm
  # rankNorm transform phenotype
  phenotype <- data.frame(donor=phenotype$donor, phenotype=rankNorm(phenotype$IFN))
  mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
    cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
    # Fit the final model on the data
    model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
    # Display regression coefficients
    modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
    colnames(modelcoef) <- c("feature","coef")
    modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
    selected_features <- modelcoef$feature
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("failed to build to model \n")
    this_result_table=c(NA,NA,NA,NA)
  } else {
    cat("successfully built the model \n")
    # If genotype is kept as a selected feature, fit the model
    if ('genotype' %in% selected_features) {
      # Build the linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
      p.adj <- NA
      pval_from_permutation <- c()
      for (i in 1:1000) {
        # permute donor indices without replacement
        shuffled_donors <- sample(rownames(mat.fs), nrow(mat.fs), replace = FALSE) %>% as.numeric
        # permute data matrix genotype based on the shuffled donor indices
        mat.permutated <- mat.fs %>% mutate(genotype=mat.fs$genotype[shuffled_donors])
        # fit permuted matrix
        model.permutated <- lm(phenotype ~ ., mat.permutated)
        pval_genotype.permutated <- summary(model.permutated)$coefficients["genotype", "Pr(>|t|)"]
        # write to variable
        pval_from_permutation <- c(pval_from_permutation, pval_genotype.permutated)
      }
      # calculated p.adj
      num_of_pval_more_extreme = length( pval_from_permutation[ which(pval_from_permutation <= p.nominal)])
      p.adj <- num_of_pval_more_extreme / length(pval_from_permutation)
      cat("Permutation testing completed. \n")
      
      this_result_table=as.character(c(p.nominal,p.adj,beta,stdErr))
      return(this_result_table)
      # header: IFNnorm_pnominal, IFNnorm_pperm, IFNnorm_beta, IFNnorm_se
    } else {
      this_result_table=c(NA,NA,NA,NA)
      return(this_result_table)
    }
  }
}


### Write output table.
res_log2FC=build_linear_model_log2FCnorm_permutation()
res_PBS=build_linear_model_PBSnorm_permutation()
res_IFN=build_linear_model_IFNnorm_permutation()

result_table=c(snp,g,res_log2FC,res_PBS,res_IFN)

output_fname=paste0(Dir,"/","modelingResult","_",combination,"/",snp,"_",g,".txt")

write.table(paste(result_table, collapse="\t"),file=output_fname,quote=F, sep="\t", row.names=F, col.names=F)
cat("Done! Cya!!!! =D \n")



# ###### unused code
# ### build linear regression model with lasso regularization (alpha=1)
# build_linear_model_log2FC_nopermutation <- function() {
#   fudge=10
#   phenotype=inner_join(PBS,IFN,by="donor") %>% mutate(phenotype=log2((IFN+fudge) / (PBS+fudge))) %>% dplyr::select(-c(PBS,IFN))
#   mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
#   # first try if model building can be successful before proceeding
#   flag <- TRUE
#   tryCatch({
#     # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
#     cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
#     # Fit the final model on the data
#     model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
#     # Display regression coefficients
#     modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
#     colnames(modelcoef) <- c("feature","coef")
#     modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
#     selected_features <- modelcoef$feature
#   },
#   error = function(e){flag<<-FALSE})
#   if (!flag) {
#     cat("failed to build to model \n")
#     this_result_table=c(NA,NA,NA,NA)
#   } else {
#     cat("successfully built the model \n")
#     # If genotype is kept as a selected feature, fit the model
#     if ('genotype' %in% selected_features){
#       # Build the linear regression model after feature selection
#       mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
#       model.fs <- lm(phenotype ~ ., mat.fs)
#       p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
#       beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
#       stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
#       res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
#       p.adj <- NA
#       all_features <- data.frame(feature=colnames(mat)[2:length(colnames(mat))] )
#       x <- as.data.frame(summary(model.fs)$coefficients)
#       x$feature <- rownames(x)
#       p_list <- data.frame(pval = x[-1,"Pr(>|t|)"])
#       p_list$feature <- x[-1,"feature"] %>% gsub("`","",.)
#       all_features <- left_join(all_features, p_list, by="feature")
#       this_result_table=c(p.nominal,p.adj,beta,stdErr,all_features$pval)
#     } else {
#       cat("genotype did not pass feature selection! \n")
#       this_result_table=c(NA,NA,NA,NA,rep(NA,ncol(mat)-1))}
#   }
#   return(this_result_table)
# }
# build_linear_model_PBS_nopermutation <- function() {
#   phenotype=PBS
#   colnames(phenotype)=c("donor","phenotype")
#   mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
#   # first try if model building can be successful before proceeding
#   flag <- TRUE
#   tryCatch({
#     # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
#     cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
#     # Fit the final model on the data
#     model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
#     # Display regression coefficients
#     modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
#     colnames(modelcoef) <- c("feature","coef")
#     modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
#     selected_features <- modelcoef$feature
#   },
#   error = function(e){flag<<-FALSE})
#   if (!flag) {
#     cat("failed to build the model \n")
#     this_result_table=c(NA,NA,NA,NA)
#   } else {
#     cat("successfully built the model, testing genotype in feature selection \n")
#     # If genotype is kept as a selected feature, fit the model
#     if ('genotype' %in% selected_features){
#       # Build the linear regression model after feature selection
#       mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
#       model.fs <- lm(phenotype ~ ., mat.fs)
#       p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
#       beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
#       stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
#       res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
#       p.adj <- NA
#       all_features <- data.frame(feature=colnames(mat)[2:length(colnames(mat))] )
#       x <- as.data.frame(summary(model.fs)$coefficients)
#       x$feature <- rownames(x)
#       p_list <- data.frame(pval = x[-1,"Pr(>|t|)"])
#       p_list$feature <- x[-1,"feature"] %>% gsub("`","",.)
#       all_features <- left_join(all_features, p_list, by="feature")
#       this_result_table=c(p.nominal,p.adj,beta,stdErr,all_features$pval)
#     } else {
#       cat("genotype did not pass feature selection! \n")
#       this_result_table=c(NA,NA,NA,NA,rep(NA,ncol(mat)-1))
#     }
#   }
#   return(this_result_table)
# }
# build_linear_model_IFN_nopermutation <- function() {
#   phenotype=IFN
#   colnames(phenotype)=c("donor","phenotype")
#   mat <- inner_join(phenotype, inner_join(genotype, covariates, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
#   # first try if model building can be successful before proceeding
#   flag <- TRUE
#   tryCatch({
#     # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
#     cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
#     # Fit the final model on the data
#     model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
#     # Display regression coefficients
#     modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
#     colnames(modelcoef) <- c("feature","coef")
#     modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
#     selected_features <- modelcoef$feature
#   },
#   error = function(e){flag<<-FALSE})
#   if (!flag) {
#     cat("failed to build the model \n")
#     this_result_table=c(NA,NA,NA,NA)
#   } else {
#     cat("successfully built the model \n")
#     # If genotype is kept as a selected feature, fit the model
#     if ('genotype' %in% selected_features){
#       # Build the linear regression model after feature selection
#       mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
#       model.fs <- lm(phenotype ~ ., mat.fs)
#       p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
#       beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
#       stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
#       res <- data.frame(residual=summary(model.fs)$residuals) %>% t()
#       p.adj <- NA
#       all_features <- data.frame(feature=colnames(mat)[2:length(colnames(mat))] )
#       x <- as.data.frame(summary(model.fs)$coefficients)
#       x$feature <- rownames(x)
#       p_list <- data.frame(pval = x[-1,"Pr(>|t|)"])
#       p_list$feature <- x[-1,"feature"] %>% gsub("`","",.)
#       all_features <- left_join(all_features, p_list, by="feature")
#       this_result_table=c(p.nominal,p.adj,beta,stdErr,all_features$pval)
#     } else {
#       cat("genotype did not pass feature selection! \n")
#       this_result_table=c(NA,NA,NA,NA,rep(NA,ncol(mat)-1))
#     }
#   }
#   return(this_result_table)
# }
