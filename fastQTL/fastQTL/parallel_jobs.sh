#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/%J%I.peer_log2FC.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/%J%I.peer_log2FC.err"
### script to interface Gencove to retrieve, upload, and process genotyping data

### crystal shan 09/2021
### bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash

############################
export PATH=$PATH:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz
phenotypeF1=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_cpm.MEL_PBS.bed.gz
phenotypeF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz
phenotypeF3=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log2FC_p0.05.MEL.bed.gz
covF1=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_cpm_MEL_PBS_covariates.txt
covF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt
covF3=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log2FC_p0.05_MEL_covariates.txt

###########################
# TMM-normalized log2 CPM from Melanocyte-PBS samples
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF2} --permute 1000 --normal --cov ${covF2} --out permutations_tmm_normalized_log2cpm.MEL_PBS.txt.gz --commands 40 commands.40.txt
while read c; do
     echo "
        export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin;
	genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz;
	phenotypeF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz;
	covF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt;
        ${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10000] -R "select[rh=6]" -q long \
	-o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.out" |\
     echo "submitted job"
done < commands.40.txt
# make sure to rename files

###########################
# log2 fold-change to detect responseQTL
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF3} --permute 1000 --normal --cov ${covF3} --out permutations_log2FC_p0.05.MEL.txt.gz --commands 40 commands.log2FC.40.txt
while read c; do
     echo "
        export PATH=${PATH}:/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/FastQTL/bin;
        genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/filtered.vcf.gz;
        phenotypeF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm.MEL_PBS.bed.gz;
        covF2=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/tmm_normalized_log2cpm_MEL_PBS_covariates.txt;
        ${c}" |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.log2FC.out" |\
     echo "submitted job"
done < commands.log2FC.40.txt

###########################
# log2 fold-change to detect responseQTL
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL
fastQTL --vcf ${genotypeF} --bed ${phenotypeF1} --permute 1000 --normal --cov ${covF1} --out permutations_tmm_normalized_cpm.MEL_PBS.txt.gz --commands 40 commands.cpm.40.txt
while read c; do
     echo ${c} |\
     bsub -W 128:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long \
        -o "/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/fastQTL/log/%J%I.permutation.cpm.out" |\
     echo "submitted job"
done < commands.cpm.40.txt
