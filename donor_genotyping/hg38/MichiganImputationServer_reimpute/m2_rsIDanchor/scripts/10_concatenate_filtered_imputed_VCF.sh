#!/bin/bash
# shuo.shan@umassmed.edu
#BSUB -J combine_imputed_vcfs
#BSUB -R "rusage[mem=120000]"
#BSUB -o step10_concatFilteredVCF_%I.out
#BSUB -e step10_concatFilteredVCF_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 8

# Load necessary modules
module load bcftools
module load htslib

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor
VCF_DIR=${DIR}/MIS_results/combined_output
OUT_DIR=${VCF_DIR}/final_output
cd ${DIR}/MIS_results/combined_output/final_output


# ----------- concatenate genotype data
# Chromosome list (space-separated).
CHROM_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

# VCF filename prefix pattern (must include "chr" placeholder).
HG38_PREFIX="chr{CHR}_imputed_combined_hg38_rsID_filtered.vcf.gz"

# Output file bases
HG38_OUT_BASE=skineQTL_imputed_hg38_filtered

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

# ---- Step 2: Concatenate hg38 VCFs ----
echo "Concatenating hg38 VCFs..."
HG38_VCF_LIST=$(build_vcf_list ${HG38_PREFIX})
HG38_OUT=${OUT_DIR}/${HG38_OUT_BASE}.vcf

bcftools concat -a -O v ${HG38_VCF_LIST} -o ${HG38_OUT}
bgzip -@ 8 ${HG38_OUT}
tabix ${HG38_OUT}.gz
