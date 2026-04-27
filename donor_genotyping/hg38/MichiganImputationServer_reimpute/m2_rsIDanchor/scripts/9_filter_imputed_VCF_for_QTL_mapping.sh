#!/bin/bash
# shuo.shan@umassmed.edu
#BSUB -J combine_imputed_vcfs[1-23]
#BSUB -R "rusage[mem=120000]"
#BSUB -o step9_filterVCF_%I.out
#BSUB -e step9_filterVCF_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load bcftools
module load htslib

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor
VCF_DIR=${DIR}/MIS_results/combined_output
CHROM_LIST=({1..22} X)
chr_name="chr${CHROM_LIST[$((LSB_JOBINDEX-1))]}"
cd ${DIR}/MIS_results/combined_output/final_output

# ----------- filter genotype data --------------
echo "Processing ${chr_name}..."

VCF1=${VCF_DIR}/${chr_name}_imputed_combined_hg38_rsID.vcf.gz
VCF2=${VCF_DIR}/${chr_name}_imputed_combined_hg38_rsID_filtered.vcf.gz
bcftools +fill-tags ${VCF1} -- -t MAF,F_MISSING |\
        bcftools view -i 'MAF>=0.08 && F_MISSING<=0.1' ${VCF1} -Oz -o ${VCF2}
tabix ${VCF2}

n1=$(bcftools view -H ${VCF1} | wc -l | cut -d' ' -f1)
n2=$(bcftools view -H ${VCF2} | wc -l | cut -d' ' -f1)
pct=$(echo "${n2}/${n1}*100" | bc -l)
pct=$(printf "%.2f" ${pct})

echo -e "${chr_name}\t${n1}\t${n2}\t${pct}" >> ${DIR}/log/step9_summary.txt
