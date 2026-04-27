#!/bin/bash
#BSUB -J concat_impute
#BSUB -R "rusage[mem=10000]"
#BSUB -o step6_concat_all_chr_hg19.out
#BSUB -e step6_concat_all_chr_hg19.err
#BSUB -n 8
#BSUB -q short
#BSUB -W 02:00

module load bcftools
module load plink2/alpha6.1amd
module load htslib

# ---- PARAMETERS  ----
VCF_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output
OUT_DIR=${VCF_DIR}/final_output
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

# Chromosome list (space-separated).
CHROM_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

# VCF filename prefix pattern (must include "chr" placeholder).
# E.g., for files like "chr1_imputed_combined_hg19.vcf.gz", use:
HG19_PREFIX="chr{CHR}_imputed_combined_hg19_rsID.vcf.gz"

# Output file bases
HG19_OUT_BASE=skineQTL_imputed_hg19

# ---- Function: Expand file list for bcftools concat ----
build_vcf_list() {
  local pattern="$1"
  local list=""
  for CHR in $CHROM_LIST; do
    file="${VCF_DIR}/$(echo ${pattern} | sed "s/{CHR}/${CHR}/")"
    if [[ -f ${file} ]]; then
      list+="$file "
    else
      echo "Warning: Missing file for chr$CHR: $file" >&2
    fi
  done
  echo "$list"
}

# ---- Step 1: Concatenate hg19 VCFs ----
HG19_VCF_LIST=$(build_vcf_list ${HG19_PREFIX})
HG19_OUT=${OUT_DIR}/${HG19_OUT_BASE}.vcf
echo "Concatenating hg19 VCFs..."
bcftools concat -a -O v ${HG19_VCF_LIST} -o ${HG19_OUT}
echo "Compressing...";date
bgzip -@ 8 ${HG19_OUT}
echo "Indexing...";date
tabix ${HG19_OUT}.gz
