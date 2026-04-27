#!/bin/bash
# parse gwas stats file by chromosome per trait

dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/analysis/fine-mapping/coloc
cd $dir
trait=$1 # vitiligo_verma2024
chr=$2	 # 1...22, X

# 
echo $trait
cd ${dir}/${trait}
f=$(ls standardized_*.h.tsv.gz)
echo "Input: ${f}"

# grab header once
header="$(zcat "${f}" | head -n1)"

# split into chr1..chr22 + chrX (adjust if you need Y/MT)
out_gz="${dir}/${trait}/standardized_${trait}.chr${chr}.tsv.gz"
echo "  -> chr${chr}"

{
  echo "${header}"
  zcat "${f}" | awk -F'\t' -v OFS='\t' -v c="${chr}" 'NR>1 && $6==c'
} | gzip -c > "${out_gz}"
