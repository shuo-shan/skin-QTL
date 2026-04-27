#!/bin/bash
#BSUB -J convert_to_rsID[1-23]
#BSUB -R "rusage[mem=120000]"
#BSUB -o step4_rename_to_rsID_hg19_%I.out
#BSUB -e step4_rename_to_rsID_hg19_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load plink
module load plink2/alpha6.1amd
module load bcftools
module load htslib

#### -------- convert to rsID -------------
# ---- PARAMETERS  ----
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute
VCF_DIR=${DIR}/m2_rsIDanchor/MIS_results/combined_output
cd ${DIR}/data
chr_name=$(awk -v idx=${LSB_JOBINDEX} 'NR==idx {print $2}' ${DIR}/data/refseq_to_chr_hg19.txt)
refseq_id=$(awk -v idx=${LSB_JOBINDEX} 'NR==idx {print $1}' ${DIR}/data/refseq_to_chr_hg19.txt)

# ---- download resources ---- #
# dbSNP index files for hg19
#cd ${DIR}/data
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

# -----------------------------------------------------------
# Step 1: Extract relevant fields from dbSNP VCF
# -----------------------------------------------------------
mkdir -p ${DIR}/data/dbSNP_hg19
cd ${DIR}/data/dbSNP_hg19

## Subset dbSNP by RefSeq accession
#echo "Processing $refseq_id -> $chr_name"; date
#bcftools view -r ${refseq_id} GCF_000001405.25.gz -Oz -o dbSNP_hg19_${chr_name}.vcf.gz
#tabix -p vcf dbSNP_hg19_${chr_name}.vcf.gz
#
## Rename to UCSC-style chromosome names
#echo "renaming to USCS-style chromosome names"; date
#bcftools annotate --rename-chrs refseq_to_chr_hg19.txt \
#dbSNP_hg19_${chr_name}.vcf.gz -Oz -o dbSNP_hg19_${chr_name}_renamed.vcf.gz
#tabix -p vcf dbSNP_hg19_${chr_name}_renamed.vcf.gz

# -----------------------------------------------------------
# Step 2: Build ID mapping file
# -----------------------------------------------------------
# Compile map file (2 columns: chr:pos:ref:alt rsID)
# features: single alt allele per row; considers allele swap
# I had created two files: vitiligo_snps.txt, and refseq_to_chr.txt
echo "starting Step 2: Build ID mapping file"; date
cd ${DIR}/data/dbSNP_hg19
bcftools view -H dbSNP_hg19_${chr_name}_renamed.vcf.gz | grep -v '^NT_' | grep -v '^NW_' | cut -f1-5 | sort -k3 > dbSNP_hg19_${chr_name}_snps.txt

# For each record:
#   - Split multiple ALT alleles
#   - For each allele:
#       * Print "chr:pos:ref:alt   rsID"  (canonical orientation)
#       * Print "chr:pos:alt:ref   rsID"  (swapped orientation, in case REF/ALT flipped)
#   - Remove "chr" prefix (so matches your imputed VCF IDs)
awk '{
  split($5, alts, ",");
  for (i in alts) {
    alt = alts[i];
    if (alt != ".") {
      # canonical orientation
      print $1":"$2":"$4":"alt"\t"$3;
      # swapped alleles orientation
      print $1":"$2":"alt":"$4"\t"$3;
    }
  }
}' dbSNP_hg19_${chr_name}_snps.txt | sed 's/chr//g' > idmap_hg19_${chr_name}.txt

# -----------------------------------------------------------
# Step 3: Apply ID mapping to imputed VCF
# -----------------------------------------------------------
echo "starting Step 3: Apply ID mapping to imputed VCF"; date

cd ${DIR}/data/dbSNP_hg19
IDMAP=${DIR}/data/dbSNP_hg19/idmap_hg19_${chr_name}.txt
VCF_IN=${VCF_DIR}/${chr_name}_imputed_combined_hg19.vcf.gz
VCF_OUT=${VCF_DIR}/${chr_name}_imputed_combined_hg19_rsID.vcf.gz

# Extract VCF header
bcftools view -h "${VCF_IN}" > temp.header.${chr_name}

# Replace IDs in the VCF body
bcftools view -H ${VCF_IN} | awk -v mapfile="$IDMAP" '
BEGIN {
    FS=OFS="\t"
    # load mapping file into array
    while ((getline < mapfile) > 0) {
        old=$1; new=$2;
        idmap[old]=new;
    }
}
# header lines
/^#/ { print; next }

{
    # if ID is in map, replace, if not in map, keep hg19 chr:pos:ref:alt
    if ($3 in idmap) {
        $3=idmap[$3];
    }
    print
}' | awk '!seen[$3]++' > temp.body.${chr_name}

# Recombine header + body and bgzip compress, index
cat temp.header.${chr_name} temp.body.${chr_name} | bgzip -c > ${VCF_OUT}
tabix -p vcf ${VCF_OUT}
rm temp.header.${chr_name} temp.body.${chr_name}

# -----------------------------------------------------------
# Step 4: Convert VCF to PLINK bed/bim/fam
# -----------------------------------------------------------
echo "starting step 4: convert VCF to PLINK files"; date

cd ${VCF_DIR}
plink --vcf ${VCF_OUT} \
      --make-bed \
      --out hg19_${chr_name}_rsID \
      --double-id

