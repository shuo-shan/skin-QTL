#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=100000]
#BSUB -q long
#BSUB -W 124:00
#BSUB -e "./vcf2bed.%J%I.err"
#BSUB -o "./vcf2bed.%J%I.out"

dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG
jobname=vcf2bed
master_vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz

module load bcftools/1.16
source activate fastQTL

cd ${dir}

# option 1: filter the vcf file by list of snps
#ln -s /pi/manuel.garber-umw/human/skin/eQTLs/GWAS_SNPs/snponly-GWAS-catalog-all-associations-autosomal-snp.txt
#bcftools view --include ID==@snponly-GWAS-catalog-all-associations-autosomal-snp.txt ${master_vcf} -Oz -o ${dir}/snponly-GWAS-catalog-all-associations-autosomal-snp.vcf.gz
#genotype_table=${dir}/snponly-GWAS-catalog-all-associations-autosomal-snp.vcf.gz

# option 2: no-filtering the vcf file by list of snps
genotype_table=${master_vcf}


# obtain the genotype and info from vcf file 
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" -H ${genotype_table} | head -1 > header
cat header | tr '\t' '\n' | cut -d']' -f2 | tr '\n' '\t'  | sed '$s/\t$/\n/' > header2
cat header2 | sed 's/POS\tPOS/START\tEND/g'  > header3
bcftools query -f "%CHROM\t%POS\t%POS\t%ID\t%REF\t%ALT{0}[\t%GT]\n" ${genotype_table} -o temp.filtered.genotype.vcf
cat temp.filtered.genotype.vcf | awk '{OFS="\t"}{print $1,$2,$2+1,$4,$5,$6}' > temp1
cat temp.filtered.genotype.vcf | cut -f7- > temp2
paste temp1 temp2 > temp3
cat temp3 > temp.filtered.genotype.bed
cat header3 > temp.header
rm header header2 header3 temp.filtered.genotype.vcf temp1 temp2 temp3
echo "there are " $(wc -l temp.filtered.genotype.bed) " GWAS SNPs in our data with at least 3 ALT counts".

### further filtering ( snps -->  snps) 
# this ensures that all the SNPs we are looking at have at least 3 donors in each of the 3 genotypes.
# split bed file into 50K entries per file.
split -l 30000 temp.filtered.genotype.bed temp.split.
sleep 10
rm commands.txt
for f in temp.split.*;do
  echo "date;Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/scripts/further_filtering_SNPs.R ${dir} ${f} temp.header; date" >> commands.txt
done
while read c; do
  echo ${c} | bsub -W 72:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -q long -e "./log/further_filtering.%J%I.err" -o "./log/further_filtering.%J%I.out"
done < commands.txt
while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; sleep 5; done
# all outputs are named new.xxxx. without a header. join them and then add a header.
cat temp.header new.temp.split.* > filtered.genotype.bed
#mv filtered.genotype.bed gwas-snps_filtered.genotype.bed
mv filtered.genotype.bed master_filtered_genotype.bed
rm temp.header temp.split.* new.temp.split.* temp.filtered.genotype.bed commands.txt
echo "there are " $(wc -l gwas-snps_filtered.genotype.bed) " GWAS SNPs with all 3 genotypes and at least 3 donors per genotype".
# this results in a filtered.genotype.bed that contain SNPs with all 3 genotypes, and at least 3 donors per genotype.

