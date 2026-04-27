#!/bin/bash
#BSUB -n 26
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=18000]
#BSUB -q long
#BSUB -W 121:00
### script to read genotype data into plink
#### bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash

#### packages & env
module load condas/2018-05-11
source activate sshan_isoform
#conda install -c bioconda plink #install ver1.90

#############################
# PLINK tutorial
cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink
wget https://zzz.bwh.harvard.edu/plink/hapmap1.zip # genotype .ped, .map file, 2 phenotype files .phe
# download command-line help
plink --help > plink-help.txt
# checking input file
plink --file hapmap1 --noweb
plink --ped hapmap1.ped --map hapmap1.map --noweb

##############################
# all donors
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
### make sampleList for the 28 donors to keep
ls /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/fastq | grep -v 'VB' | grep -v 'CB' | grep -v 'F47' > all_donors.txt
cat all_donors.txt | awk '{OFS="\t"}{print $1,$1,1}' > all_donors_info.txt
grep -f /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/samples.txt all_donors_info.txt > 28_donors_info.txt 
sampleList=all_donors_info.txt
### remove non-biallelic variants & duplicated variants
bcftools view -H ${genotypeF} | cut -f3 | sort | uniq -d > duplicated
plink --vcf ${genotypeF} --biallelic-only strict list --update-sex ${sampleList} --recode --snps-only --maf 0.10 --exclude duplicated --out 37donors
plink --vcf ${genotypeF} --recode --keep-fam /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/samples.txt --update-sex ${sampleList} --exclude duplicated --snps-only --biallelic-only strict list --maf 0.10 --out 28donors
### obtain skin-disease GWAS SNPs
pref=28donors_snps_with_ld
gwas_snps=/nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/skin_disease_GWAS_SNPs.txt
plink --file 28donors --r2 --ld-snp-list ${gwas_snps} --ld-window-kb 1000 --ld-window-r2 0.6 --out ${pref}
# SNPs with high LD info are stored in: 28donors_snps_with_ld.ld
### organize GWAS SNPs and all those in high LD into a single column
cat 28donors_snps_with_ld.ld | awk '{OFS=";"}NR>1{print $3,$6}' | tr ';' '\n' | sort | uniq > 28donors_snps_with_ld_lst.txt


##############################
this=rs706779
grep $this 28donors_snps_with_ld.ld

###############################
# 11/02/2021
# get 600 SNPs surrounding rs706779 in 37 donors file
cd /nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink
f=37donors.map
cat ${f} | awk '{if ($1==10) print $2,$1,$4}' > temp.37donors.chr10.txt
cat temp.37donors.chr10.txt | grep -w -n rs706779 #result: 14630
head -14630 temp.37donors.chr10.txt | tail -300 | cut -d' ' -f1 > temp.rs706779.300SNPs.upstream.txt
head -14929 temp.37donors.chr10.txt | tail -300 | cut -d' ' -f1 > temp.rs706779.300SNPs.downstream.txt
# run ldmatrix with these lists 
# 
