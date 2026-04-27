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

# # # #example of snp-gene pair
# Dir="~/Downloads/nl/human/skin/eQTLs/DREG/transQTL/MEL"
# snp="rs229527"
# g="FLOT2"
# combination="a10b10" # a=n_genotypePCs, b=n_latentVariables
# PBSFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/PBS_",g,".txt")
# IFNFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/IFN_",g,".txt")
# genotypeFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/genotype.txt")

##########_________main script________________________________________________
##### compile data table ##### 
PBS=data.table::fread(PBSFile, header=TRUE)
IFN=data.table::fread(IFNFile, header=TRUE)
genotype=data.table::fread(genotypeFile, header=TRUE)
covariates_log2FC=data.table::fread(paste0(Dir,"/covariates_phenotype_Log2FCwithDummy10_",combination,".txt"), header=TRUE)
covariates_PBS=data.table::fread(paste0(Dir,"/covariates_phenotype_CPM_PBS_",combination,".txt"), header=TRUE)
covariates_IFN=data.table::fread(paste0(Dir,"/covariates_phenotype_CPM_IFN_",combination,".txt"), header=TRUE)

covariates_rankNormlog2FC=data.table::fread(paste0(Dir,"/covariates_phenotype_rankNormLog2FCwithDummy10_",combination,".txt"), header=TRUE)
covariates_rankNormPBS=data.table::fread(paste0(Dir,"/covariates_phenotype_rankNormCPM_PBS_",combination,".txt"), header=TRUE)
covariates_rankNormIFN=data.table::fread(paste0(Dir,"/covariates_phenotype_rankNormCPM_IFN_",combination,".txt"), header=TRUE)
cat("compiled data matrix for modeling, building models now \n")


##### inverse normal transformation ##### 
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


##### build linear regression model with CPM as independent variable, and lasso regularization (alpha=1) with permutation ##### 
build_linear_model_log2FC_with_permutation <- function() {
  fudge=10
  phenotype=inner_join(PBS,IFN,by="donor") %>% mutate(phenotype=log2((IFN+fudge) / (PBS+fudge))) %>% dplyr::select(-c(PBS,IFN))
  mat <- inner_join(phenotype, inner_join(genotype, covariates_log2FC, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # stability selection. only continue if genotype is selected feature in 90% of the 100 5-fold-cross-validation runs. 
    genotype_selected_table <- rep(FALSE,100)
    for (i in 1:100) {
      # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
      cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
      # not used due to variability. Find the best lambda using bootstrapping.
      #bestLambda <- bootstrap_and_find_lambda(mat[,-1], mat[,1], 100)
      # Fit the final model on the data
      model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
      # Display regression coefficients
      modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
      colnames(modelcoef) <- c("feature","coef")
      modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
      selected_features <- modelcoef$feature
      if ("genotype" %in% selected_features) {
        genotype_selected_table[i] <- TRUE
      }
      genotype_selection_freq <- length(which(genotype_selected_table==T))/100
    }
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("reQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
  } else {
    cat("reQTL: successfully built the model...")
    # If genotype is kept as a selected feature and stable, fit the model
    if ("genotype" %in% selected_features & genotype_selection_freq >= 0.90){
      cat("genotype is a stable selected feature \n")
      Sys.time()
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      
      Sys.time()
      this_result_table_minimal=c(p.nominal,p.adj,beta,stdErr)
      
      # Build the feature selected linear regression model
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]

      Sys.time()
      this_result_table_featureSelected=c(p.nominal,p.adj,beta,stdErr)
      
      # join two result tables
      this_result_table=c(this_result_table_minimal,this_result_table_featureSelected)
      
    } else {
      cat("genotype did not pass feature selection or isn't stable! \n")
      this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)}
  }
  return(this_result_table)
}
build_linear_model_PBS_with_permutation <- function() {
  phenotype=PBS
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, inner_join(genotype, covariates_PBS, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # stability selection. only continue if genotype is selected feature in 90% of the 100 5-fold-cross-validation runs. 
    genotype_selected_table <- rep(FALSE,100)
    for (i in 1:100) {
      # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
      cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
      # not used due to variability. Find the best lambda using bootstrapping.
      #bestLambda <- bootstrap_and_find_lambda(mat[,-1], mat[,1], 100)
      # Fit the final model on the data
      model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
      # Display regression coefficients
      modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
      colnames(modelcoef) <- c("feature","coef")
      modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
      selected_features <- modelcoef$feature
      if ("genotype" %in% selected_features) {
        genotype_selected_table[i] <- TRUE
      }
      genotype_selection_freq <- length(which(genotype_selected_table==T))/100
    }
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("PBSeQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
  } else {
    cat("PBSeQTL: successfully built the model...")
    # If genotype is kept as a selected feature in at least 90% of the cv runs, fit the minimal model
    if ("genotype" %in% selected_features & genotype_selection_freq >= 0.90){
      cat("genotype is a stable selected feature \n")
      Sys.time()
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]

      Sys.time()
      this_result_table_minimal=c(p.nominal,p.adj,beta,stdErr)
      
      # Build the feature selected linear regression model
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]

      Sys.time()
      this_result_table_featureSelected=c(p.nominal,p.adj,beta,stdErr)
      
      # join two result tables
      this_result_table=c(this_result_table_minimal,this_result_table_featureSelected)
      
    } else {
      cat("genotype did not pass feature selection or isn't stable! \n")
      this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
    }
  }
  return(this_result_table)
}
build_linear_model_IFN_with_permutation <- function() {
  phenotype=IFN
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, inner_join(genotype, covariates_IFN, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # stability selection. only continue if genotype is selected feature in 90% of the 100 5-fold-cross-validation runs. 
    genotype_selected_table <- rep(FALSE,100)
    for (i in 1:100) {
      # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
      cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
      # not used due to variability. Find the best lambda using bootstrapping.
      #bestLambda <- bootstrap_and_find_lambda(mat[,-1], mat[,1], 100)
      # Fit the final model on the data
      model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
      # Display regression coefficients
      modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
      colnames(modelcoef) <- c("feature","coef")
      modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
      selected_features <- modelcoef$feature
      if ("genotype" %in% selected_features) {
        genotype_selected_table[i] <- TRUE
      }
      genotype_selection_freq <- length(which(genotype_selected_table==T))/100
    }
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("IFNeQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
  } else {
    cat("IFNeQTL: successfully built the model...")
    # If genotype is kept as a selected feature in at least 90% of the cv runs, fit the minimal model
    if ("genotype" %in% selected_features & genotype_selection_freq >= 0.90){
      cat("genotype is a stable selected feature \n")
      Sys.time()
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
 
      Sys.time()
      this_result_table_minimal=c(p.nominal,p.adj,beta,stdErr)
      
      # Build the feature selected linear regression model
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]

      Sys.time()
      this_result_table_featureSelected=c(p.nominal,p.adj,beta,stdErr)
      
      # join two result tables
      this_result_table=c(this_result_table_minimal,this_result_table_featureSelected)
      
    } else {
      cat("genotype did not pass feature selection or isn't stable! \n")
      this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
      }
  }
  return(this_result_table)
}

##### build linear regression model with rankNormalizedValue as independent variable, and lasso regularization (alpha=1) with permutation ##### 
build_linear_model_rankNormlog2FC_with_permutation <- function() {
  fudge=10
  phenotype=inner_join(PBS,IFN,by="donor") %>% mutate(phenotype=log2((IFN+fudge) / (PBS+fudge))) %>% dplyr::select(-c(PBS,IFN))
  # rankNorm transform phenotype
  phenotype <- data.frame(donor=phenotype$donor, phenotype=rankNorm(phenotype$phenotype))
  mat <- inner_join(phenotype, inner_join(genotype, covariates_rankNormlog2FC, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # stability selection. only continue if genotype is selected feature in 90% of the 100 5-fold-cross-validation runs. 
    genotype_selected_table <- rep(FALSE,100)
    for (i in 1:100) {
      # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
      cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
      # not used due to variability. Find the best lambda using bootstrapping.
      #bestLambda <- bootstrap_and_find_lambda(mat[,-1], mat[,1], 100)
      # Fit the final model on the data
      model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
      # Display regression coefficients
      modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
      colnames(modelcoef) <- c("feature","coef")
      modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
      selected_features <- modelcoef$feature
      if ("genotype" %in% selected_features) {
        genotype_selected_table[i] <- TRUE
      }
      genotype_selection_freq <- length(which(genotype_selected_table==T))/100
    }
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("rankNormreQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
  } else {
    cat("rankNormreQTL: successfully built the model...")
    # If genotype is kept as a selected feature and stable, fit the model
    if ("genotype" %in% selected_features & genotype_selection_freq >= 0.90){
      cat("genotype is a stable selected feature \n"); Sys.time()
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      Sys.time()
      this_result_table_minimal=c(p.nominal,p.adj,beta,stdErr)
      
      # Build the feature selected linear regression model
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      Sys.time()
      this_result_table_featureSelected=c(p.nominal,p.adj,beta,stdErr)
      
      # join two result tables
      this_result_table=c(this_result_table_minimal,this_result_table_featureSelected)
      
    } else {
      cat("genotype did not pass feature selection or isn't stable! \n")
      this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)}
  }
  return(this_result_table)
}
build_linear_model_rankNormPBS_with_permutation <- function() {
  phenotype=PBS.norm
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, inner_join(genotype, covariates_rankNormPBS, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # stability selection. only continue if genotype is selected feature in 90% of the 100 5-fold-cross-validation runs. 
    genotype_selected_table <- rep(FALSE,100)
    for (i in 1:100) {
      # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
      cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
      # not used due to variability. Find the best lambda using bootstrapping.
      #bestLambda <- bootstrap_and_find_lambda(mat[,-1], mat[,1], 100)
      # Fit the final model on the data
      model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
      # Display regression coefficients
      modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
      colnames(modelcoef) <- c("feature","coef")
      modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
      selected_features <- modelcoef$feature
      if ("genotype" %in% selected_features) {
        genotype_selected_table[i] <- TRUE
      }
      genotype_selection_freq <- length(which(genotype_selected_table==T))/100
    }
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("rankNormPBSeQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
  } else {
    cat("rankNormPBSeQTL: successfully built the model...")
    # If genotype is kept as a selected feature in at least 90% of the cv runs, fit the minimal model
    if ("genotype" %in% selected_features & genotype_selection_freq >= 0.90){
      cat("genotype is a stable selected feature \n"); Sys.time()
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      Sys.time()
      this_result_table_minimal=c(p.nominal,p.adj,beta,stdErr)
      
      # Build the feature selected linear regression model
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      Sys.time()
      this_result_table_featureSelected=c(p.nominal,p.adj,beta,stdErr)
      
      # join two result tables
      this_result_table=c(this_result_table_minimal,this_result_table_featureSelected)
      
    } else {
      cat("genotype did not pass feature selection or isn't stable! \n")
      this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
    }
  }
  return(this_result_table)
}
build_linear_model_rankNormIFN_with_permutation <- function() {
  phenotype=IFN.norm
  colnames(phenotype)=c("donor","phenotype")
  mat <- inner_join(phenotype, inner_join(genotype, covariates_rankNormIFN, by="donor"), by="donor") %>% dplyr::select(-donor) %>% as.matrix()
  # first try if model building can be successful before proceeding
  flag <- TRUE
  tryCatch({
    # stability selection. only continue if genotype is selected feature in 90% of the 100 5-fold-cross-validation runs. 
    genotype_selected_table <- rep(FALSE,100)
    for (i in 1:100) {
      # Find the best lambda using cross-validation. lambda is a tuning parameter that controls the overall penalty strength.
      cv <- cv.glmnet(mat[,-1], mat[,1], alpha = 1, nfolds=5, standardize=TRUE)
      # not used due to variability. Find the best lambda using bootstrapping.
      #bestLambda <- bootstrap_and_find_lambda(mat[,-1], mat[,1], 100)
      # Fit the final model on the data
      model <- glmnet(mat[,-1], mat[,1], alpha = 1, lambda = cv$lambda.min)
      # Display regression coefficients
      modelcoef <- coef(model)[,] %>% as.data.frame() %>% rownames_to_column("feature")
      colnames(modelcoef) <- c("feature","coef")
      modelcoef <- modelcoef %>% filter(coef != 0 & feature != "(Intercept)")
      selected_features <- modelcoef$feature
      if ("genotype" %in% selected_features) {
        genotype_selected_table[i] <- TRUE
      }
      genotype_selection_freq <- length(which(genotype_selected_table==T))/100
    }
  },
  error = function(e){flag<<-FALSE})
  if (!flag) {
    cat("rankNormIFNeQTL: failed to build the model \n")
    this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
  } else {
    cat("rankNormIFNeQTL: successfully built the model...")
    # If genotype is kept as a selected feature in at least 90% of the cv runs, fit the minimal model
    if ("genotype" %in% selected_features & genotype_selection_freq >= 0.90){
      cat("genotype is a stable selected feature \n"); Sys.time()
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      Sys.time()
      this_result_table_minimal=c(p.nominal,p.adj,beta,stdErr)
      
      # Build the feature selected linear regression model
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]

      this_result_table_featureSelected=c(p.nominal,p.adj,beta,stdErr)
      
      # join two result tables
      this_result_table=c(this_result_table_minimal,this_result_table_featureSelected)
      
    } else {
      cat("genotype did not pass feature selection or isn't stable! \n")
      this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)
    }
  }
  return(this_result_table)
}


##### Write output table. ##### 
res_log2FC = build_linear_model_log2FC_with_permutation()
res_PBS = build_linear_model_PBS_with_permutation()
res_IFN = build_linear_model_IFN_with_permutation()
res_rankNormlog2FC = build_linear_model_rankNormlog2FC_with_permutation()
res_rankNormPBS = build_linear_model_rankNormPBS_with_permutation()
res_rankNormIFN = build_linear_model_rankNormIFN_with_permutation()

result_table=c(snp,g,
               res_log2FC,res_rankNormlog2FC, 
               res_PBS,res_rankNormPBS,
               res_IFN, res_rankNormIFN)
result_table[is.na(result_table)] <- "."

# # column names for result_table
# colnames(result_table)[c(1,2)] <- c("SNP","GENE")
# colnames(result_table)[c(3,4,5,6)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("log2FC_minimalModel_reQTL_genotype_",x)))
# colnames(result_table)[c(7,8,9,10)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("log2FC_featureSelected_reQTL_genotype_",x)))
# colnames(result_table)[c(11,12,13,14)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormlog2FC_minimalModel_reQTL_genotype_",x)))
# colnames(result_table)[c(15,16,17,18)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormlog2FC_featureSelected_reQTL_genotype_",x)))
# 
# colnames(result_table)[c(19,20,21,22)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_minimalModel_PBSeQTL_genotype_",x)))
# colnames(result_table)[c(23,24,25,26)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_featureSelected_PBSeQTL_genotype_",x)))
# colnames(result_table)[c(27,28,29,30)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_minimalModel_PBSeQTL_genotype_",x)))
# colnames(result_table)[c(31,32,33,34)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_featureSelected_PBSeQTL_genotype_",x)))
# 
# colnames(result_table)[c(35,36,37,38)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_minimalModel_IFNeQTL_genotype_",x)))
# colnames(result_table)[c(39,40,41,42)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("CPM_featureSelected_IFNeQTL_genotype_",x)))
# colnames(result_table)[c(43,44,45,46)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_minimalModel_IFNeQTL_genotype_",x)))
# colnames(result_table)[c(47,48,49,50)] <- unlist(lapply(c("pval","pperm","beta","se"),function(x) paste0("rankNormCPM_featureSelected_IFNeQTL_genotype_",x)))


output_fname=paste0(Dir,"/","modelingResult","_",combination,"/",snp,"_",g,".txt")

write.table(paste(result_table, collapse="\t"),file=output_fname,quote=F, sep="\t", row.names=F, col.names=F)
cat("Done! Cya!!!! =D \n")