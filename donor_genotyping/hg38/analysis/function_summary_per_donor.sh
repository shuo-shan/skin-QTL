#!/bin/bash
module load bcftools/1.9
dir=/nl/umw_manuel_garber/human/skin/eQTLs/donor_genotyping/hg38/vcf/organized
donor=$1

cd ${dir}
f=${dir}/${donor}.vcf.gz
echo "step 1: extract all snps" ; date;
bcftools view -v snps ${f} -Oz -o ${dir}/temp/${donor}.snps.vcf.gz 
echo "step 2: extract all snps with direct observed genotype" ; date;
bcftools query -f "%ID[\t%GP\t%GT\t%DS]\n" ${dir}/temp/${donor}.snps.vcf.gz | awk '{if ( $2~/1,1e-10,1e-10/ || $2~/1e-10,1,1e-10/ || $2~/1e-10,1e-10,1/ ) print $0}' > ${dir}/temp/${donor}.snps.obs.vcf
echo "step 3: extract all snps with imputed observed genotype"; date;
bcftools query -f "%ID\t%FILTER[\t%GP\t%DS]\n" ${dir}/temp/${donor}.snps.vcf.gz | awk '{if ( !( $3~/1,1e-10,1e-10/ || $3~/1e-10,1,1e-10/ || $3~/1e-10,1e-10,1/ ) ) print $0}' > ${dir}/temp/${donor}.snps.imputed.vcf
echo "step 4: extract all snps with imputed observed genotype that failed quality filter"; date;
cat ${dir}/temp/${donor}.snps.imputed.vcf | awk '{if ($2~/LOW/) print $0}' > ${dir}/temp/${donor}.snps.imputed.lowconf.txt


all=$(bcftools view -H ${dir}/temp/${donor}.snps.vcf.gz | wc -l | cut -d' ' -f1) 
obs=$(wc -l ${dir}/temp/${donor}.snps.obs.vcf |  cut -d' ' -f1)
impfail=$(wc -l ${dir}/temp/${donor}.snps.imputed.lowconf.txt | cut -d' ' -f1)

echo -e ${donor}"\t"${all}"\t"${obs}"\t"${impfail} >> ${dir}/summary.txt
rm ${dir}/temp/${donor}.snps.vcf.gz ${dir}/temp/${donor}.snps.obs.vcf ${dir}/temp/${donor}.snps.imputed.vcf ${dir}/temp/${donor}.snps.imputed.lowconf.txt

