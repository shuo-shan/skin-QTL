#!/bin/bash
#BSUB -n 20
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=22500]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./query.%J%I.err"
#BSUB -o "./query.%J%I.out"

module load bcftools/1.9

cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/

cat /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/01.24.2022/skin_disease_GWAS_SNPs_06112022.txt > snp.list

bcftools view --include ID==@snp.list bd841628-fcc2-487a-8460-f5428237f0c9.merged.vcf.gz -Oz -o skin_disease_GWAS_SNPs_all_donors_merged.vcf.gz –-threads 19

