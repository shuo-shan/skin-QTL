#!/bin/bash
snp=$1

dir=/pi/manuel.garber-umw/human/skin/eQTLs/edQTL
cd /pi/manuel.garber-umw/human/skin/eQTLs/edQTL/output

module load plink/1.90b6.27
plinkFname=bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.plink
plink --file ${plinkFname} --r2 --ld-snp ${snp} --ld-window-kb 1000 --ld-window-r2 0.8 --out high_ld_variants_of_${snp}
cat high_ld_variants_of_${snp}.ld | awk '{OFS="\t"}NR>1{if ($3 != $6) print $3,$6,$7}' >> ${dir}/output/gwas_snps_and_LD_tags.txt
rm high_ld_variants_of_${snp}*; date;

