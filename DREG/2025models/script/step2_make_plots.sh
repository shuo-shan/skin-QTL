#!/usr/bin/env bash
set -euo pipefail
module load bcftools

# Set-up
DIR="/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL"
RES_FILE=${DIR}/results/reQTL_IFNB_pvalE-09.txt
OUT_NAME="reQTL_IFNB_pvalE-09"
PAIR_FILE="${DIR}/plots/pairs_${OUT_NAME}.txt"
OUT_PDF="${DIR}/plots/${OUT_NAME}.pdf"

cd ${DIR}/results
for f in result_*.tsv; do
  awk -F'\t' 'NR>1 && $10 != "" && $10+0 < 1e-9' "$f" >> ${RES_FILE}
done
cat ${RES_FILE} | cut -f1,2 | sort -u > "${PAIR_FILE}"

# Resources
CPM_FILE="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/CPM.sampleFiltered.metaConverted.txt"
META_FILE="/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/metadata.sampleFiltered.txt"
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz

# Intermediate files
snplist="${DIR}/plots/snps_${OUT_NAME}.txt"
vcf_subset="${DIR}/plots/snps_${OUT_NAME}.vcf.gz"
samples="${DIR}/plots/temp.${OUT_NAME}.samples.txt"
header="${DIR}/plots/temp.${OUT_NAME}.header.txt"
body="${DIR}/plots/temp.${OUT_NAME}.body.tsv"
geno="${DIR}/plots/temp.${OUT_NAME}.genotype.tsv"
geno_num="${DIR}/plots/genotype_${OUT_NAME}.tsv"

mkdir -p "${DIR}/plots"
cd ${DIR}/plots

# -------- Subset genotype VCF file to SNP(s) of interest --------------------- #
echo "==== Subsetting genotype VCF file ====";date
# 1) collect SNP IDs from this chunk
cat ${RES_FILE} | cut -f1 | sort -u > "$snplist"
# 2) extract VCF records
bcftools view --threads 5 --include ID==@"$snplist" "$vcf" -Oz -o "$vcf_subset"
# 3) clean sample names
bcftools query -l "$vcf_subset" | cut -d'_' -f1 > "$samples"
# 4) build header
printf "CHROM\tPOS\tID\tREF\tALT\t" > "$header"
paste -sd '\t' "$samples" >> "$header"
# 5) extract GT matrix
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf_subset" > "$body"
cat "$header" "$body" > "$geno"
# 6) convert genotypes to numeric
sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' "$geno" > "$geno_num"
# 7) cleanup intermediate files (leave only numeric genotype)
rm "$snplist" "$vcf_subset" "$samples" "$header" "$body" "$geno"

# -------- Set-up Singularity and Dependencies for Rscript ------ #
echo "==== Making 3x3 grid plots for SNP:gene pair(s) ====";date

export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
Rscript plot_reQTL_grid.R \
  --pairs="${PAIRS_FILE}" \
  --geno="${GENO_FILE}" \
  --out="${OUT_PDF}"
