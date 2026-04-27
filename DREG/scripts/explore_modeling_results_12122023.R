#!/usr/bin/env Rscript
# written by Crystal Shan 12/2023

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)

##### model for predicting eQTL power calculation for variants with different MAF in keratinocytes
# 1. fetch keratinocyte model result table

# 2. from the model result, get the variants with MAF = c(0.1+0.9, 0.2+0.8, 0.3+0.7, 0.4+0.6, 0.5)
#.    sigma.x = sqrt(2* MAF * (1-MAF))
# 3. for each MAF-tier, calculate the following:
#    - lambda.a = eQTL's mean slope
#    - sigma.y = eQTL linked gene's standard deviation across all x values.
# 4. use power.SLR to calculate
#    - alpha = 0.05/nTests
#    - power.SLR(n, lambda.a, sigma.x, sigma.y, alpha=0.05, verbose=TRUE)


##### minimal model power calculation n=28 ##### 
# load modeling_results_minimalModel_phenotypeRankNormCPM
dir <- "~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results"
resF <- paste0(dir,"/modeling_results_minimalModel_phenotypeRankNormCPM.txt")
res <- data.table::fread(resF)
res[res == "."] <- NA
res <- res %>% 
  mutate(across(contains(c("reQTL_","PBSeQTL_","IFNeQTL_")), as.numeric)) %>%
  mutate(across(where(is.numeric), ~as.numeric(sprintf("%.3e", .))))

# set the significance cut-off to be Bonferroni adjusted pval to threshold.
# threshold denominator is the sum of all tests in FRB, KRT, MEL
# number of tests in FRB, KRT, MEL: 8038333 + 7941046 + 7961789
threshold = 0.05 / 8038333
critical_Zalpha <- qnorm( 1 - threshold/2)
critical_Zpower <- qnorm(0.80)

# histogram of each model's slope. 
hist(res[which(res$reQTL_pval < threshold),]$reQTL_beta , breaks = 100, main="reQTL beta")
hist(res[which(res$PBSeQTL_pval < threshold),]$PBSeQTL_beta , breaks = 100, main = "PBS eQTL beta")
hist(res[which(res$IFNeQTL_pval < threshold),]$IFNeQTL_beta , breaks = 100, main = "IFN eQTL beta")

# check out the effect size
abs(res[which(res$reQTL_pval < threshold),]$reQTL_beta) %>% sort()
summary(abs(res[which(res$reQTL_pval < threshold),]$reQTL_beta)) 
summary(abs(res[which(res$PBSeQTL_pval < threshold),]$PBSeQTL_beta)) 
summary(abs(res[which(res$IFNeQTL_pval < threshold),]$IFNeQTL_beta)) 
effect_size = 1.15

# standard deviation of the residuals
summary(res[which(res$reQTL_pval < threshold),]$reQTL_se)
summary(res[which(res$PBSeQTL_pval < threshold),]$PBSeQTL_se)
summary(res[which(res$IFNeQTL_pval < threshold),]$IFNeQTL_se)
# average standard error is 0.13. 
0.13 * sqrt(35) # result: 0.7690904
var_X = 0.13 * sqrt(35)

# calculate desired sample size
n = ( (critical_Zalpha + critical_Zpower)^2 * var_X )/( effect_size^2 ) # answer = 28



##### feature selected model power calculation n=21 ######## 
# load modeling_results_featureSelectedModel_phenotypeRankNormCPM
dir <- "~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results"
resF <- paste0(dir,"/modeling_results_featureSelectedModel_phenotypeRankNormCPM.txt")
res <- data.table::fread(resF)
res[res == "."] <- NA
res <- res %>% 
  mutate(across(contains(c("reQTL_","PBSeQTL_","IFNeQTL_")), as.numeric)) %>%
  mutate(across(where(is.numeric), ~as.numeric(sprintf("%.3e", .))))

# set the significance cut-off to be Bonferroni adjusted pval to threshold.
# threshold denominator is the sum of all tests in FRB, KRT, MEL
threshold = 0.05 / (8038333 + 7941046 + 7961789)
critical_Zalpha <- qnorm( 1 - threshold/2)
critical_Zpower <- qnorm(0.80)

# histogram of each model's slope. 
hist(res[which(res$reQTL_pval < threshold),]$reQTL_beta , breaks = 100, main="reQTL beta")
hist(res[which(res$PBSeQTL_pval < threshold),]$PBSeQTL_beta , breaks = 100, main = "PBS eQTL beta")
hist(res[which(res$IFNeQTL_pval < threshold),]$IFNeQTL_beta , breaks = 100, main = "IFN eQTL beta")

# check out the effect size
abs(res[which(res$reQTL_pval < threshold),]$reQTL_beta) %>% sort()
summary(abs(res[which(res$reQTL_pval < threshold),]$reQTL_beta)) 
summary(abs(res[which(res$PBSeQTL_pval < threshold),]$PBSeQTL_beta)) 
summary(abs(res[which(res$IFNeQTL_pval < threshold),]$IFNeQTL_beta)) 
effect_size = 1.15

# standard deviation of the residuals
summary(res[which(res$reQTL_pval < threshold),]$reQTL_se)
summary(res[which(res$PBSeQTL_pval < threshold),]$PBSeQTL_se)
summary(res[which(res$IFNeQTL_pval < threshold),]$IFNeQTL_se)
# average standard error is 0.10. 
0.10 * sqrt(35) # result: 0.591608
var_X = 0.10 * sqrt(35)

# calculate desired sample size
n = ( (critical_Zalpha + critical_Zpower)^2 * var_X )/( effect_size^2 ) # answer = 21



##### checking out top hits filtered by effect size #####
# load modeling_results_featureSelectedModel_phenotypeRankNormCPM
dir <- "~/Downloads/nl/human/skin/eQTLs/website/data/modeling_results"
resF <- paste0(dir,"/modeling_results_featureSelectedModel_phenotypeRankNormCPM.txt")
res <- data.table::fread(resF)
res[res == "."] <- NA
res <- res %>% 
  mutate(across(contains(c("reQTL_","PBSeQTL_","IFNeQTL_")), as.numeric)) %>%
  mutate(across(where(is.numeric), ~as.numeric(sprintf("%.3e", .))))

# set the significance cut-off to be Bonferroni adjusted pval to threshold.
# threshold denominator is the sum of all tests in FRB, KRT, MEL
threshold = 0.05 / (8038333 + 7941046 + 7961789)

# histogram of each model's slope. 
hist(res[which(res$reQTL_pval < threshold),]$reQTL_beta , breaks = 100, main="reQTL beta")
hist(res[which(res$PBSeQTL_pval < threshold),]$PBSeQTL_beta , breaks = 100, main = "PBS eQTL beta")
hist(res[which(res$IFNeQTL_pval < threshold),]$IFNeQTL_beta , breaks = 100, main = "IFN eQTL beta")

# how many target genes are related to antigen presentation?
res[which(res$reQTL_pval < threshold),]$GENE %>% unique()
res[which(res$PBSeQTL_pval < threshold),]$GENE %>% unique()
res[which(res$PBSeQTL_pval < threshold),]$GENE %>% unique()
res[which(res$reQTL_pval < 0.000001),]$GENE %>% unique()
######## power caclulation ######## 
library(pwr)  

# Parameters
effect_size <- 1.1 / 0.13 # Cohen's d = mean slope divided by standard error of slope
sig_level <- 0.05 / 8038333 
power <- 0.80

# Calculate sample size
sample_size_result <- pwr.t.test(d = effect_size, 
                                 sig.level = sig_level, 
                                 power = power, 
                                 type = "one.sample", 
                                 alternative = "two.sided")
# effect size (Cohen's d) - difference between the means divided by the pooled standard deviation

# Output the result
print(sample_size_result)
# n=10.1343 --> 11



######## end ######## 