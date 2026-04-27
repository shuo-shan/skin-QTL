#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=50040]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# script written by Crystal SHan 03/2022
# goal: filter donor vcf file with GWAS SNPs

# In this folder:
cd /nl/umw_manuel_garber/human/skin/eQTLs/fastQTL/GWAS_SNPs/03.08.2022
cat skin_disease_snps.txt | sed 's/chr//g' | sort | uniq > skin_disease_GWAS_SNPs.txt
### get all GWAS SNPs for skin disease
module load condas/2018-05-11
source activate sshan_isoform
module load bcftools/1.9  
temp=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/plink/37donors_snps_plus_high_ld_lst.txt
genotypeF=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
bcftools view -h ${genotypeF} > header
bcftools view -H ${genotypeF} > body
split -l 100000 body temp.split.
rm commands.txt
for f in temp.split.*;do
  echo "grep -w -f ${temp} ${f} > grepped.${f}" >> commands.txt
done
while read c;do
  echo ${c} | bsub -J abdcf -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long -e "./log/%J%I.err" -o "./log/%J%I.out"
done < commands.txt
while [[ $(bjobs | grep 'abdcf' | wc -l) != 0 ]] ; do echo $(bjobs | grep 'abdcf' | wc -l) "jobs remaining"; sleep 5; done
cat grepped.temp.split.* > tempbody
cat header tempbody > 37donors_skin_disease_GWAS_SNPs_plus_high_ld.vcf
bgzip 37donors_skin_disease_GWAS_SNPs_plus_high_ld.vcf && tabix -p vcf 37donors_skin_disease_GWAS_SNPs_plus_high_ld.vcf.gz 

rm grepped.temp.split.*
rm temp.split.*
rm temp header body tempbody filtered_skin_disease_GWAS_SNPs.vcf
rm commands.txt

