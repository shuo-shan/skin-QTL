#!/bin/bash
# goal: make covariate table for modeling based on specified combination of covariates. donor is the first column
# step1. input arguments
n_gPC=$1
n_phenoPC=$2
phenotype=$3 # example: CPM_PBS, Log2FCwithDummy10, rankNormCPM_PBS, rankNormLog2FCwithDummy10

# step2. specificy the file locations
datadir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis
genotypePCs=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/pca/pca_eigenvec_table.txt
phenotypePCs=${datadir}/PCs_${phenotype}_expressedGenes.txt

# step 3. slice covariate table
# 1. genotype PCs
cat $genotypePCs | cut -f1  > temp_genotypePC_donors.txt
cat $genotypePCs | cut -f4- > temp_genotypePC_PCs.txt
awk -v j=${n_gPC} '{for (i=1;i<=j;i++) printf $i"\t"; print ""}' temp_genotypePC_PCs.txt | paste temp_genotypePC_donors.txt - | sort -k1 > genotypePCs.txt
rm temp*.txt
# 2. phenotype PCs
cat $phenotypePCs | cut -f1 > temp_phenotypePC_donors.txt
cat $phenotypePCs | cut -f2- > temp_phenotypePC_PCs.txt
awk -v j=${n_phenoPC} '{for (i=1;i<=j;i++) printf $i"\t"; print ""}' temp_phenotypePC_PCs.txt | paste temp_phenotypePC_donors.txt - | sort -k1 > phenotypePCs.txt
rm temp*.txt

# step 4. compile covariate table
join -j 1 genotypePCs.txt phenotypePCs.txt > covariates_phenotype_${phenotype}_a${n_gPC}b${n_phenoPC}.txt
rm genotypePCs.txt phenotypePCs.txt















