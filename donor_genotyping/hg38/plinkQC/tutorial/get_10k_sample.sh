#!/bin/bash

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/plinkQC
outdir=${dir}/sampled_by_chr
chr=$1
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz

### get test.vcf
module load bcftools
cd ${dir}
echo "Processing $chr"
bcftools view -r $chr $vcf | bcftools view -H | shuf -n 10000 | cut -f1,2 > ${outdir}/${chr}_random10000.txt
