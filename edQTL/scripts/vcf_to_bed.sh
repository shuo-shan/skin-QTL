#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=100000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -e "./vcf2bed.%J%I.err"
#BSUB -o "./vcf2bed.%J%I.out"

dir=/pi/manuel.garber-umw/human/skin/eQTLs/edQTL
jobname=vcf2bed
genotype_table=${dir}/output/snps.vcf.gz

module load bcftools/1.16
conda activate fastQTL

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
cat temp.header temp.filtered.genotype.bed > genotype.bed
rm header header2 header3 temp.filtered.genotype.vcf temp1 temp2 temp3 temp.header temp.filtered.genotype.bed

