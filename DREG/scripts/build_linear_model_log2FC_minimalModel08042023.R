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

# # #example of snp-gene pair
#Dir="~/Downloads/nl/human/skin/eQTLs/DREG/MEL_minimal"
#snp="rs2910686"
#g="CAST"
#combination="a10b10" # a=n_genotypePCs, b=n_latentVariables
#PBSFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/PBS_",g,".txt")
#IFNFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/IFN_",g,".txt")
#genotypeFile=paste0(Dir,"/","modelingResult_",combination,"/temp_",snp,"/genotype.txt")

##########_________main script________________________________________________
##### compile data table ##### 
PBS=data.table::fread(PBSFile, header=TRUE)
IFN=data.table::fread(IFNFile, header=TRUE)
genotype=data.table::fread(genotypeFile, header=TRUE)
covariates_log2FC=data.table::fread(paste0(Dir,"/covariates_phenotype_Log2FCwithDummy10_",combination,".txt"), header=TRUE)
covariates_PBS=data.table::fread(paste0(Dir,"/covariates_phenotype_CPM_PBS_",combination,".txt"), header=TRUE)
covariates_IFN=data.table::fread(paste0(Dir,"/covariates_phenotype_CPM_IFN_",combination,".txt"), header=TRUE)
cat("compiled data matrix for modeling, building models now \n")

##### build linear regression model with CPM as independent variable, and lasso regularization (alpha=1), without permutation  ##### 
build_linear_model_log2FC_nopermutation <- function() {
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
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
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
      this_result_table=c(NA,NA,NA,NA,NA,NA,NA,NA)}
  }
  return(this_result_table)
}
build_linear_model_PBS_nopermutation <- function() {
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
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
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
build_linear_model_IFN_nopermutation <- function() {
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
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
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
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      # Permutation for 100,000 times
      pval_from_permutation <- c()
      # Using 100,000 permutations reduces the uncertainty near p =0.05 to ±0.1% and allows p-values as small as 0.00001
      for (i in 1:100000) {
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
      this_result_table_minimal=c(p.nominal,p.adj,beta,stdErr)
      
      # Build the feature selected linear regression model
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
      # Permutation for 100,000 times
      pval_from_permutation <- c()
      # Using 100,000 permutations reduces the uncertainty near p =0.05 to ±0.1% and allows p-values as small as 0.00001
      for (i in 1:100000) {
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
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
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
      # Build the minimal linear regression model after feature selection
      mat.fs <- mat[,c(1,which(colnames(mat) %in% selected_features))] %>% as.data.frame()
      mat.fs <- mat.fs[,c(1,2)]
      model.fs <- lm(phenotype ~ ., mat.fs)
      p.nominal <- summary(model.fs)$coefficients["genotype", "Pr(>|t|)"]
      p.adj <- NA
      beta <- summary(model.fs)$coefficients["genotype", "Estimate"]
      stdErr <- summary(model.fs)$coefficients["genotype", "Std. Error"]
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
res_log2FC=build_linear_model_log2FC_nopermutation()
res_PBS=build_linear_model_PBS_nopermutation()
res_IFN=build_linear_model_IFN_nopermutation()

result_table=c(snp,g,res_log2FC,res_PBS,res_IFN)
result_table[is.na(result_table)] <- "."

output_fname=paste0(Dir,"/","modelingResult","_",combination,"/",snp,"_",g,".txt")

write.table(paste(result_table, collapse="\t"),file=output_fname,quote=F, sep="\t", row.names=F, col.names=F)
cat("Done! Cya!!!! =D \n")

##### unused code ##### 
### function to find lambda by bootstrapping and output selected features
bootstrap_and_find_lambda <- function(x, y, num_bootstrap) {
  # Number of observations in the dataset
  n <- nrow(x)
  
  # Number of features in the dataset
  p <- ncol(x)
  
  # Empty matrix to store bootstrap samples
  bootstrap_samples <- matrix(NA, nrow = n, ncol = num_bootstrap)
  
  # Empty vector to store optimal lambda for each bootstrap sample
  optimal_lambdas <- numeric(num_bootstrap)
  
  # Generate bootstrap samples
  for (i in 1:num_bootstrap) {
    bootstrap_indices <- sample(1:n, size = n, replace = TRUE)
    bootstrap_x <- x[bootstrap_indices, ]
    bootstrap_y <- y[bootstrap_indices]
    bootstrap_samples[, i] <- bootstrap_y
    
    # Fit the model using glmnet with alpha = 1 for the bootstrap sample
    model <- glmnet(bootstrap_x, bootstrap_y, alpha = 1)
    
    # Find the optimal lambda for the bootstrap sample based on MSE
    mse <- apply((predict(model, newx = bootstrap_x, s = model$lambda) - bootstrap_y)^2, 2, mean)
    optimal_lambda_idx <- which.min(mse)
    optimal_lambdas[i] <- model$lambda[optimal_lambda_idx]
  }
  
  # Calculate the average of optimal lambdas over all bootstrap samples
  avg_optimal_lambda <- mean(optimal_lambdas)
  
  return(avg_optimal_lambda)
}

