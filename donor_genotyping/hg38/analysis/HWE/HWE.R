#!/usr/bin/env Rscript
# written by Crystal Shan 08/2023
# HWE analysis

library(tidyverse)
library(magrittr)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
library("HardyWeinberg")
#vignette("HardyWeinberg")

x<-c(MM=298,MN=489,NN=213)
HW.test<-HWChisq(x,verbose=TRUE)

# Bydefault,HWChisq applies a continuity correction.
# This is not recommended for low minor allele frequencies. In order to perform a 
# chi-square test without Yates’ continuity correction, it is necessary to set 
# the cc parameter to zero:
HW.test<-HWChisq(x,cc=0,verbose=TRUE)

# If the genotype counts aa, ab, bb are collected in a m x 3 matrix, with each row 
# representing a marker, then HWEtests can be run over each row in the matrix by
# the routines HWChisqMat. These routines return a list with the p values and test
# statistics for each marker.
set.seed(123)
X2 <- HWData(100,100)
colnames(X2) <- c("MM","MN","NN")
res <- HWChisqMat(X2, cc=0)
output <- cbind(X2, res$chisqvec, res$pvalvec) %>% 
  set_colnames(c("MM","MN","NN","chisq","pval")) %>%
  as.data.frame()
# output$pval < 0.05 are the markers that deviate from HWE.
output[which(output$pval < 0.05),]


# If I wish to perform all possible tests for a marker, I run:
HW.results<-HWAlltests(x,verbose=TRUE,include.permutation.test=TRUE)


# simulate data: HWData
# by default generate simulated genotype data with multinomial distribution around HWE
# if specified I could generate data with exact HWE. (exactequilibrium=TRUE)
set.seed(123)
n <- 100 # 100 donors
m <- 100 # 100 markers
X1 <- HWData(m,n,p=rep(0.5,m))
X2 <- HWData(m,n) # 100 markers under HWE with a random uniform allele frequency 
X3 <- HWData(m,n,exactequilibrium=TRUE) # 100 markers under exact HWE with a random uniform allele frequency

opar<-par(mfrow=c(2,2),mar=c(1,0,3,0)+0.1)
par(mfg=c(1,1))
HWTernaryPlot(X1,main="(a)",vbounds=FALSE)
par(mfg=c(1,2))
HWTernaryPlot(X2,main="(b)",vbounds=FALSE)
par(mfg=c(2,1))
HWTernaryPlot(X3,main="(c)",vbounds=FALSE)
par(opar)




