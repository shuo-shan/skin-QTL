#!/bin/bash

f="$1"
dir=$(dirname "$f")
prefix=$(basename "$f" .vcf.gz)

echo "${f}"
echo "${dir}"
echo "${prefix}"
rm ${dir}/${prefix}_lowGPresetGT.vcf.gz* 

module load bcftools
echo "resetting now for ${prefix}"; date
bcftools +setGT ${f} -- -t q -i 'max(GP)<0.90' -n "./." | bcftools view -Oz -o ${dir}/${prefix}_lowGPresetGT.vcf.gz
echo "resetting GT done, indexing now"; date
bcftools index ${dir}/${prefix}_lowGPresetGT.vcf.gz
