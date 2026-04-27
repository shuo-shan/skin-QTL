#!/bin/bash
# convert genotype VCF hg38->hg19 by first rsID anchor then liftover variants without rsID
# author: shuo.shan@umassmed.edu
#BSUB -J rsIDanchor_then_liftover_hg38tohg19_chr6[1-23]
#BSUB -R "rusage[mem=120000]"
#BSUB -o hg38tohg19_%I.out
#BSUB -e hg38tohg19_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load plink
module load plink2/alpha6.1amd
module load bcftools
module load htslib
module load picard/3.1.1
module load samtools

# ------------------------
# Global Parameters
# ------------------------
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor
DATA_DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/data
chr_name=$(awk -v idx=${LSB_JOBINDEX} 'NR==idx {print $2}' ${DATA_DIR}/refseq_to_chr_hg19.txt)
refseq_id=$(awk -v idx=${LSB_JOBINDEX} 'NR==idx {print $1}' ${DATA_DIR}/refseq_to_chr_hg19.txt)
mkdir -p ${DIR}/data
mkdir -p ${DIR}/log

# -------------------------------------------------------------------------- #
# ---- STEP 0: Create plink file set of hg38 preMIS VCF file            ---- #
# -------------------------------------------------------------------------- #
echo "step0: create plink2 file set for preMIS hg38 VCF";date

cd ${DATA_DIR}
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz

# Define the 2 donors to remove: they failed plinkQC heterogeneity test
echo "skineQTL-F086,skineQTL-F104" | tr "," "\n" > samples_to_remove.txt
cat samples_to_remove.txt | awk '{OFS=FS="\t"}{print "0",$0}' > samples_to_remove_plink.txt

# Make a 3-column sex file: FID  IID  SEX (1=male, 2=female, 0=unknown), required by plink2
bcftools query -l ${vcf} \
  | awk 'BEGIN{OFS="\t"} {print 0, $1, 1}' > sex.update

# Create PLINK2 file set
BFILE_HG38="${DATA_DIR}/cohort.hg38.preMIS"
plink2 \
	--vcf ${vcf} \
	--split-par b38 \
	--remove samples_to_remove_plink.txt \
	--update-sex sex.update \
	--make-bed \
	--out ${BFILE_HG38}

# -------------------------------------------------------------------------- #
# ---- STEP 1. CONFIRM INPUT hg38 GENOTYPE VARIANTS HAVE rsID,-------------- # 
# ----         SUBSET TO THIS CHR                             -------------- #
# -------------------------------------------------------------------------- #
echo "step1: confirm input hg38 variants have rsID in ${chr_name}";date
#cat ${BFILE_HG38}.bim | awk '$2 ~ /^rs/' | wc -l
#cat ${BFILE_HG38}.bim | awk '$2 !~ /^rs/' | wc -l
# confirmed: 98.5% of variants in input VCF have rsID.

# Convert PLINK to VCF
cd ${DIR}/data
# subset hg38 vcf to this chromosome
VCF_HG38_PREFIX=preMIS_hg38
bcftools view -S ^${DATA_DIR}/samples_to_remove.txt -r ${chr_name} ${vcf} -Oz -o ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}.vcf.gz
tabix -p vcf ${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}.vcf.gz

# --------------------------------------------------- #
# ---- STEP 2. HG38->HG19 VIA rsID ANCHOR METHOD ---- #
# --------------------------------------------------- #
echo "step2";date
DBSNP_VCF=${DATA_DIR}/dbSNP_hg19/dbSNP_hg19_${chr_name}_renamed.vcf.gz
IDMAP_TXT=${DATA_DIR}/dbSNP_hg19/idmap_hg19_${chr_name}.txt
# 2.1) Download resources
# dbSNP index files for hg19
#cd ${DATA_DIR}
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
#
## 2.2) Extract relevant fields from dbSNP VCF
## Subset dbSNP by RefSeq accession corresponding to chromosome of interest
#cd ${DATA_DIR}/dbSNP_hg19
#echo "Processing $refseq_id -> $chr_name"; date
#bcftools view -r ${refseq_id} ${DATA_DIR}/GCF_000001405.25.gz -Oz -o ${DATA_DIR}/dbSNP_hg19/dbSNP_hg19_${chr_name}.vcf.gz
#tabix -p vcf ${DATA_DIR}/dbSNP_hg19/dbSNP_hg19_${chr_name}.vcf.gz
#
## Rename dbSNP chromosome to UCSC-style from refseq-style
#echo "renaming to UCSC-style chromosome name"; date
#bcftools annotate --rename-chrs ${DATA_DIR}/refseq_to_chr_hg19.txt \
#${DATA_DIR}/dbSNP_hg19/dbSNP_hg19_${chr_name}.vcf.gz -Oz -o ${DBSNP_VCF}
#tabix -p vcf ${DBSNP_VCF}
#
# 2.3) Build rsID -> hg19 coordinate mapping file
# Compile map file (3 columns: rsID chr pos)
echo "Building rsID→hg19 mapping for ${chr_name}"; date
bcftools view -H ${DBSNP_VCF} \
       | cut -f1-3 \
       | awk '{OFS=FS="\t"}{print $3,$1,$2}' \
       | sed 's/chr//g' \
       > ${IDMAP_TXT}
echo "Done: dbSNP_hg19_${chr_name}_renamed.vcf.gz + ${IDMAP_TXT}"

# Step 2.4) For variants with rsID in mapping file: 
# Reannotate genome positions to hg19 from dbSNP mapping made in previous step
cd ${DIR}/data
THIS_VCF_HG38=${DIR}/data/${VCF_HG38_PREFIX}_${chr_name}.vcf.gz
echo "Annotating with hg19 coordinates now using rsID anchor now"; date

# 2.4.1. Build new header
# Take old header, drop hg38 contigs
bcftools view -h ${THIS_VCF_HG38} | grep -v "^##contig=" > header_${chr_name}.tmp
# Split old header: meta-info vs column header
head -n -1 header_${chr_name}.tmp > header_${chr_name}_part1.tmp
tail -n 1 header_${chr_name}.tmp > header_${chr_name}_part2.tmp
# Extract hg19 contigs from dbSNP file (only chr-prefixed contigs), drop "chr" prefix to match hg19 reference
bcftools view -h ${DBSNP_VCF} | grep "^##contig=<ID=chr" | sed 's/chr//g' >> header_${chr_name}_part1.tmp
cat header_${chr_name}_part1.tmp header_${chr_name}_part2.tmp > header_${chr_name}.tmp

# 2.4.2. Replace chr/pos based only on rsID
# if variant has corresponding rsID, update their position to hg19 position
# if variant doesn't have rsID, keep it for later liftOver
bcftools view -H ${THIS_VCF_HG38} | \
	sed 's/^chr//' | \
awk -v OFS="\t" -v mapfile=${IDMAP_TXT} -v chr_name=${chr_name} '
BEGIN {
  # load mapping file
  while ((getline < mapfile) > 0) {
    id=$1; chr=$2; pos=$3;
    coord[id]=chr"\t"pos;
  }
  missingFile=chr_name"_no_rsID_variants.tmp"
  rsIDanchoredFile="body_"chr_name".tmp"
  system("rm -f "out) # clean old files
}
{
  split($3, ids, ";")   # split on ";". col#3 of VCF = ID field.
  foundAny=0
  for (i in ids) {
    vcf_id=ids[i]

    if (vcf_id in coord) {
	    # Found in mapping --> update CHR & POS
      split(coord[vcf_id], f, "\t")
      $1=f[1]   # new CHR
      $2=f[2]   # new POS
      $3=vcf_id # replace ID with this rsID
      print > rsIDanchoredFile
      foundAny=1
    }
  }
  if (foundAny==0) {
    # no matching rsIDs in coord → print original line with all IDs.
    # not in mapping --> write entire line of hg38 position to missing file
    print > missingFile
  }
}
'
nVarTotal=$(bcftools view -H ${THIS_VCF_HG38} | cut -f3 | tr ";" "\n" | wc -l) # get total number of variants (also count merged variants)
nVarHasrsID=$(wc -l body_${chr_name}.tmp | cut -d' ' -f1)
nVarNorsID=$(wc -l ${chr_name}_no_rsID_variants.tmp | cut -d' ' -f1)
pct=$(echo "${nVarHasrsID}/(${nVarHasrsID} + ${nVarNorsID})*100" | bc -l)
pct=$(printf "%.2f" ${pct})
echo "Step2: Before MIS imputation, among ${nVarTotal} variants (hg38, ${chr_name}): ${nVarHasrsID} variants have rsID in hg19 (${pct}%). ${nVarNorsID} do not. Don't 100% add up because of merged variants." > ${DIR}/log/script0_${chr_name}.log
cat ${DIR}/log/script0_${chr_name}.log

# -----------------------------------------------------------
# Step 3: Generate hg19 PLINK set from rsID-anchored VCF
# -----------------------------------------------------------
echo "step3";date
cd ${DIR}/data
rsIDanchored_hg19_vcf=${DIR}/data/${chr_name}_varWithRsID_hg19.vcf.gz

# Recombine header and body, strip "chr" prefix (hg19 reference uses bare 1..22, X, Y, MT)
cat body_${chr_name}.tmp | awk '{OFS=FS="\t"}{sub(/^chr/, "", $1); print}' > body_${chr_name}.tmp.tmp
cat header_${chr_name}.tmp body_${chr_name}.tmp.tmp | bgzip -c > ${rsIDanchored_hg19_vcf}

# Sort and index the VCF for downstream tools
bcftools sort -Oz -o ${rsIDanchored_hg19_vcf}.sorted.gz ${rsIDanchored_hg19_vcf}
tabix -f -p vcf ${rsIDanchored_hg19_vcf}.sorted.gz

# Creating new PLINK1.9 or PLINK2 file sets will lose the sex information. Add it here.
awk '{print $1, $2, $5}' ${BFILE_HG38}.fam > sex_info_${chr_name}.txt

# Convert VCF to PLINK2 format (pgen/pvar/psam), keeping only ACGT SNPs
OUTPAR=${chr_name}.hg19.acgt
plink2 \
	--vcf ${rsIDanchored_hg19_vcf}.sorted.gz \
	--snps-only just-acgt \
	--make-pgen \
	--update-sex sex_info_${chr_name}.txt \
	--split-par b37 \
	--out ${OUTPAR}

# Log the number of variants remaining after applying just-acgt filter:
txt=$(cat ${OUTPAR}.log | grep "variants loaded from")
top=$(echo ${txt} | cut -d' ' -f1)
bottom=$(echo ${txt} | cut -d' ' -f4)
pct=$(echo "${top}/${bottom}*100" | bc -l)
pct=$(printf "%.2f" ${pct})
echo "Step3: After applying just-acgt filter, ${txt}, (${pct}%)">> ${DIR}/log/script0_${chr_name}.log

# Convert PLINK2 set back to PLINK1.9 bed/bim/fam for Rayner tools
plink2 \
        --pfile ${OUTPAR} \
        --make-bed \
        --out ${OUTPAR}

# Standardize chromosome codes to match hg19/1kG conventions:
# X=23, Y=24, MT=26. PAR1/PAR2/XY merged into 23.
awk 'BEGIN{OFS="\t"}
{
  c=$1
  if (c=="PAR1" || c=="PAR2" || c=="XY") c=23;   # treat PAR records as X (coords are on X)
  else if (c=="X")  c=23;
  else if (c=="Y")  c=24;
  else if (c=="MT" || c=="M") c=26;
  $1=c; print
}' "${OUTPAR}.bim" > ${OUTPAR}.chr_renamed.bim

# Overwrite original BIM file with corrected chromosome codes
mv ${OUTPAR}.chr_renamed.bim ${OUTPAR}.bim

# -----------------------------------------------------------
# Step 4: Prepare dataset for Rayner's HRC-1000G-check-bim
# -----------------------------------------------------------

echo "step4";date
# Apply basic QC (remove variants failing HWE) and make clean PLINK set
plink --bfile ${OUTPAR} \
  --hwe 0.001 \
  --make-bed \
  --out panelprep_${chr_name}

# Log
top=$(cat panelprep_${chr_name}.log | grep "pass filters and QC" | cut -d' ' -f1)
bottom=$(cat panelprep_${chr_name}.log | grep "variants loaded from" | cut -d' ' -f1)
pct=$(echo "${top}/${bottom}*100" | bc -l)
pct=$(printf "%.2f" ${pct})
echo "Step4: After filtering by HWE 0.001, ${top} out of ${bottom} variants remaining, (${pct}%)" >> ${DIR}/log/script0_${chr_name}.log

# Generate frequency file required by Rayner check-bim
plink --bfile panelprep_${chr_name} --freq --out panelprep_${chr_name}

# -----------------------------------------------------------
# Step 5: Run Rayner check-bim against HRC reference panel
# -----------------------------------------------------------

echo "step5: running Rayner against HRC ref panel";date
# Create working directory and symlink panelprep files
mkdir -p ${DIR}/pre_imputation_vcf_HRC
mkdir -p ${DIR}/pre_imputation_vcf_HRC/${chr_name}
cd ${DIR}/pre_imputation_vcf_HRC/${chr_name}
ln -s ${DIR}/data/panelprep_${chr_name}.bed .
ln -s ${DIR}/data/panelprep_${chr_name}.bim .
ln -s ${DIR}/data/panelprep_${chr_name}.fam .
ln -s ${DIR}/data/panelprep_${chr_name}.frq .
ln -s ${DIR}/data/panelprep_${chr_name}.log .

# Path to HRC reference sites file (hg19/GRCh37 coordinates)
REF_HRC="${DATA_DIR}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"

# Run the check-bim script (generates Run-plink.sh with correction commands)
perl ${DIR}/../bin/HRC-1000G-check-bim.pl \
  -b panelprep_${chr_name}.bim \
  -f panelprep_${chr_name}.frq \
  -r ${REF_HRC} \
  -o ${DIR}/pre_imputation_vcf_HRC/${chr_name} \
  -h

# Apply Rayner's suggested fixes for HRC
bash Run-plink.sh

# Log number of variants remaining after Rayner tool
mv LOG-panelprep_${chr_name}-HRC.txt ${DIR}/log/script0_${chr_name}_Rayner_HRC.log
log1=${DIR}/pre_imputation_vcf_HRC/${chr_name}/panelprep_${chr_name}-updated.log
log2=${DIR}/log/script0_${chr_name}_Rayner_HRC.log
top=$(cat ${log1} | grep "pass filters and QC" | cut -d' ' -f1)
bottom=$(cat ${log2} | grep "Total processed" | cut -d' ' -f3)
pct=$(echo "${top}/${bottom}*100" | bc -l)
pct=$(printf "%.2f" ${pct})
idx=$(cat ${log2} | grep -n "Matching to HRC" | cut -d':' -f1)
echo "Step5: After applying the Rayner tool against HRC ref panel, ${top} out of ${bottom} variants remaining, (${pct}%)" >> ${DIR}/log/script0_${chr_name}.log
echo -e "------Step5 details-------\n $(awk -v idx=${idx} 'NR>=idx' ${log2}) \n -------------" >> ${DIR}/log/script0_${chr_name}.log

# Gather all VCF files
OUTDIR=${DIR}/pre_imputation_vcf_HRC/vcf
mkdir -p ${OUTDIR}
f=panelprep_${chr_name}-updated-chr${LSB_JOBINDEX}.vcf
mv ${f} ${OUTDIR}
cd ${OUTDIR}
bcftools view --threads 8 -Oz -o ${f}.gz ${f}
tabix -f -p vcf ${f}.gz
rm ${f}

# -----------------------------------------------------------
# Step 6: Run Rayner check-bim against 1000 Genomes Phase 3 panel
# -----------------------------------------------------------

echo "step6";date
# Create working directory and symlink panelprep files
mkdir -p ${DIR}/pre_imputation_vcf_1000G
mkdir -p ${DIR}/pre_imputation_vcf_1000G/${chr_name}
cd ${DIR}/pre_imputation_vcf_1000G/${chr_name}
ln -s ${DIR}/data/panelprep_${chr_name}.bed .
ln -s ${DIR}/data/panelprep_${chr_name}.bim .
ln -s ${DIR}/data/panelprep_${chr_name}.fam .
ln -s ${DIR}/data/panelprep_${chr_name}.frq .
ln -s ${DIR}/data/panelprep_${chr_name}.log .

# Path to 1000G Phase3 reference panel (legend file with strand/orientation info)
REF_1000G="${DATA_DIR}/1000GP_Phase3_combined.legend.gz"

# Run the check-bim script against 1000G reference
perl ${DIR}/../bin/HRC-1000G-check-bim.pl \
  -b panelprep_${chr_name}.bim \
  -f panelprep_${chr_name}.frq \
  -r ${REF_1000G} \
  -o ${DIR}/pre_imputation_vcf_1000G/${chr_name} \
  -g -p ALL

# Apply Rayner's suggested fixes for 1000G
bash Run-plink.sh

# Log number of variants remaining after Rayner tool
mv LOG-panelprep_${chr_name}-1000G.txt ${DIR}/log/script0_${chr_name}_Rayner_1000G.log
log1=${DIR}/pre_imputation_vcf_1000G/${chr_name}/panelprep_${chr_name}-updated.log
log2=${DIR}/log/script0_${chr_name}_Rayner_1000G.log
top=$(cat ${log1} | grep "pass filters and QC" | cut -d' ' -f1)
bottom=$(cat ${log2} | grep "Total processed" | cut -d' ' -f3)
pct=$(echo "${top}/${bottom}*100" | bc -l)
pct=$(printf "%.2f" ${pct})
idx=$(cat ${log2} | grep -n "Matching to 1000G" | cut -d':' -f1)
echo "Step6: After applying the Rayner tool against 1000G ref panel, ${top} out of ${bottom} variants remaining, (${pct}%)" >> ${DIR}/log/script0_${chr_name}.log
echo -e "------Step6 details-------\n $(awk -v idx=${idx} 'NR>=idx' ${log2}) \n -------------" >> ${DIR}/log/script0_${chr_name}.log


# Gather all VCF files
OUTDIR=${DIR}/pre_imputation_vcf_1000G/vcf
mkdir -p ${OUTDIR}
f=panelprep_${chr_name}-updated-chr${LSB_JOBINDEX}.vcf
mv ${f} ${OUTDIR}
cd ${OUTDIR}
bcftools view --threads 8 -Oz -o ${f}.gz ${f}
tabix -f -p vcf ${f}.gz 
rm ${f}

# -----------------------------------------------------------
# Step 7: Extract chr6 region (28–34 Mb) for HLA imputation
# -----------------------------------------------------------

echo "step7"; date
if [ "${chr_name}" == "chr6" ]; then
    echo "Extracting chr6 (28–34 Mb) for HLA imputation..."

    mkdir -p ${DIR}/pre_imputation_HLA
    cd ${DIR}/pre_imputation_HLA

    # Use the 1000G-corrected dataset from Step 6 (after Rayner)
    # Replace <OUTNAME> with the basename created by Run-plink.sh
    BFILE=${DIR}/pre_imputation_vcf_1000G/${chr_name}/panelprep_${chr_name}
    OUT=${DIR}/pre_imputation_HLA/chr6_28_34Mb_hg19

    plink2 --bfile ${BFILE} \
      --chr 6 \
      --from-bp 28000000 --to-bp 34000000 \
      --hwe 0.001 \
      --snps-only just-acgt \
      --max-alleles 2 \
      --recode vcf bgz \
      --out ${OUT}

    tabix -f -p vcf ${OUT}.vcf.gz

    echo "Submit ${OUT}.vcf.gz to the Michigan Imputation Server HLA pipeline."


    # -------------------------------------------------------
    # Optional: handle MIS-reported bad SNPs
    # -------------------------------------------------------
    SNPS_TO_REMOVE=${DIR}/pre_imputation_HLA/statistics/snps-excluded.txt
    if [ -f "$SNPS_TO_REMOVE" ]; then
        echo "Removing SNPs flagged by MIS and regenerating chr6 VCF..."
        OUTPUT_PREFIX=${OUT}_updated

        plink2 --bfile ${BFILE} \
          --chr 6 \
          --from-bp 28000000 --to-bp 34000000 \
          --hwe 0.001 \
          --exclude ${SNPS_TO_REMOVE} \
          --snps-only just-acgt \
          --max-alleles 2 \
          --recode vcf bgz \
          --out ${OUTPUT_PREFIX}

        tabix -f -p vcf ${OUTPUT_PREFIX}.vcf.gz
    else
        echo "No SNPs-to-exclude.txt found yet. Skipping exclusion step."
    fi
fi


# -----------------------------------------------------------
# Step 8: Clean-up
# -----------------------------------------------------------
if [ ${chr_name} != "chr6" ]; then
  cd ${DIR}/data
  #rm dbSNP_hg19_${chr_name}_renamed.vcf.gz*
  #rm dbSNP_hg19_${chr_name}.vcf.gz*
  #rm preMIS_hg38_${chr_name}.vcf.gz* 
  rm ${chr_name}_no_rsID_variants.tmp
  rm header_${chr_name}_part1.tmp
  rm header_${chr_name}_part2.tmp
  rm header_${chr_name}.tmp
  rm body_${chr_name}.tmp
  rm body_${chr_name}.tmp.tmp
  rm ${DATA_DIR}/dbSNP_hg19/idmap_hg19_${chr_name}.txt
  rm ${chr_name}_varWithRsID_hg19.vcf.gz*
  rm ${chr_name}.hg19.acgt*
  rm sex_info_${chr_name}.txt
  rm panelprep_${chr_name}.bed
  rm panelprep_${chr_name}.bim
  rm panelprep_${chr_name}.f*
  rm -r ${DIR}/pre_imputation_vcf_HRC/${chr_name}
  rm -r ${DIR}/pre_imputation_vcf_1000G/${chr_name}
fi
