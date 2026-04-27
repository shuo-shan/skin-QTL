#!/bin/bash
#BSUB -J rsID_to_hg38[1-23]
#BSUB -R "rusage[mem=120000]"
#BSUB -o step5_get_hg38_coordinates_via_rsIDanchor_%I.out
#BSUB -e step5_get_hg38_coordinates_via_rsIDanchor_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load plink
module load plink2/alpha6.1amd
module load bcftools
module load htslib
module load picard/3.1.1

# ------------------------
# Parameters
# ------------------------
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute
DATA_DIR=${DIR}/data
IN_DIR=${DIR}/m2_rsIDanchor/MIS_results/combined_output
OUT_DIR=${DIR}/m2_rsIDanchor/MIS_results/combined_output
VCF_DIR=${DIR}/m2_rsIDanchor/MIS_results/combined_output

chr_name=$(awk -v idx=${LSB_JOBINDEX} 'NR==idx {print $2}' ${DATA_DIR}/refseq_to_chr_hg38.txt)
refseq_id=$(awk -v idx=${LSB_JOBINDEX} 'NR==idx {print $1}' ${DATA_DIR}/refseq_to_chr_hg38.txt)

IN_VCF=${IN_DIR}/${chr_name}_imputed_combined_hg19_rsID.vcf.gz
OUT_VCF=${OUT_DIR}/${chr_name}_imputed_combined_hg38_rsID.vcf.gz

mkdir -p ${DATA_DIR}/dbSNP_hg38
DBSNP_VCF=${DATA_DIR}/dbSNP_hg38/dbSNP_${chr_name}_hg38_renamed.vcf.gz
IDMAP_TXT=${DATA_DIR}/dbSNP_hg38/idmap_${chr_name}_hg38.txt

### ---- download resources ---- #
#cd ${DIR}/data
# dbSNP GRCh38
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
#wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi

# ------------------------
# Step 0: Create isolated working directory
# ------------------------
TMPID=${chr_name}_${LSB_JOBID:-$$}
WORKDIR=${OUT_DIR}/tmp/${TMPID}
mkdir -p "${WORKDIR}"
cd "${WORKDIR}" || { echo "Cannot enter ${WORKDIR}"; exit 1; }

echo "Working in isolated tmp dir: ${WORKDIR}"

## -----------------------------------------------------------
## Step 1: Extract relevant fields from dbSNP VCF
## -----------------------------------------------------------
## Subset dbSNP by RefSeq accession
#echo "Processing $refseq_id -> $chr_name"; date
#bcftools view -r ${refseq_id} GCF_000001405.40.gz -Oz -o dbSNP_${chr_name}_hg38.vcf.gz
#tabix -p vcf dbSNP_${chr_name}_hg38.vcf.gz
#
## Rename to UCSC-style chromosome names
#echo "renaming to USCS-style chromosome names"; date
#bcftools annotate --rename-chrs ${DIR}/data/refseq_to_chr_hg38.txt \
#dbSNP_${chr_name}_hg38.vcf.gz -Oz -o dbSNP_${chr_name}_hg38_renamed.vcf.gz
#tabix -p vcf dbSNP_${chr_name}_hg38_renamed.vcf.gz
#
# -----------------------------------------------------------
# Step 2: Build rsID → hg38 coordinate mapping file
# -----------------------------------------------------------
cd ${DIR}/data
echo "Building rsID→hg38 mapping for $chr_name"; date
bcftools view -H ${DIR}/data/dbSNP_hg38/dbSNP_${chr_name}_hg38_renamed.vcf.gz \
	| cut -f1-3 \
	| awk '{OFS=FS="\t"}{print $3,$1,$2}' \
	> $IDMAP_TXT
echo "Done: dbSNP_${chr_name}_hg38_renamed.vcf.gz + $IDMAP_TXT"

# -----------------------------------------------------------
# Step 3: For variants with rsID in mapping file: 
#         Reannotate genome positions to hg38 from dbSNP mapping made in previous step
# -----------------------------------------------------------
echo "Annotating with hg38 coordinates now"; date

cd ${DIR}/data/dbSNP_hg38
## Build new header
# Take old header, drop hg19 contigs
bcftools view -h ${IN_VCF} | grep -v "^##contig=" > header_${chr_name}_hg19.tmp
# Split: meta-info vs column header
head -n -1 header_${chr_name}_hg19.tmp > header_${chr_name}_part1.tmp
tail -n 1 header_${chr_name}_hg19.tmp > header_${chr_name}_part2.tmp
# Extract hg38 contigs (only chr-prefixed contigs)
bcftools view -h ${DBSNP_VCF} | grep "^##contig=<ID=chr" >> header_${chr_name}_part1.tmp
cat header_${chr_name}_part1.tmp header_${chr_name}_part2.tmp > header_${chr_name}.tmp
# Clean-up 
rm header_${chr_name}_hg19.tmp header_${chr_name}_part1.tmp header_${chr_name}_part2.tmp

# Replace chr/pos based only on rsID
# if variant has corresponding rsID, update their position to hg38 position
# if variant doesn't have rsID, keep it for later liftOver
bcftools view -H ${IN_VCF} | \
awk -v OFS="\t" -v mapfile="$IDMAP_TXT" -v chr_name="$chr_name" '
BEGIN {
  # load mapping file
  while ((getline < mapfile) > 0) {
    id=$1; chr=$2; pos=$3; 
    coord[id]=chr"\t"pos;
  }
  # open file handles for 2 cases: variants without rsID, and variants with rsID
  missingFile = chr_name"_no_rsID_variants.tmp"
  out="body_"chr_name".tmp"
  system("rm -f "out) # clean old files
}
{
  vcf_id=$3   # column 3 of VCF = ID field
  if (vcf_id in coord) {
    # Found in mapping → update CHR & POS
    split(coord[vcf_id],f,"\t");
    $1=f[1];   # CHR (hg38)
    $2=f[2];   # POS (hg38)
    # ID stays as rsID
    print > "body_"chr_name".tmp"
  } else {
    # Not in mapping → write entire line to missing file
    print > missingFile
  }
}
'
total=$(bcftools view -H ${IN_VCF} | wc -l | cut -d' ' -f1)
noRSID=$(wc -l ${chr_name}_no_rsID_variants.tmp | cut -d' ' -f1)
pct=$(echo "scale=3; ${noRSID}/${total}*100" | bc)
echo "there are ${noRSID} variants out of ${total} without hg38 annotation in ${chr_name} data (${pct}%)" > ${DIR}/m2_rsIDanchor/log/hg19tohg38_rsIDanchor_${chr_name}.txt
cat ${DIR}/m2_rsIDanchor/log/hg19tohg38_rsIDanchor_${chr_name}.txt
mv body_${chr_name}.tmp body_all_${chr_name}.tmp
# head ${chr_name}_no_rsID_variants.tmp | cut -f1-12
# head body_${chr_name}.tmp | cut -f-12

# -----------------------------------------------------------
# Step 5: Finalize VCF (sort + compress + index)
# -----------------------------------------------------------
# Split header into meta-info vs column header
echo "Finalizing hg38 VCF for ${chr_name}"; date
TMPID=${chr_name}_${LSB_JOBID:-$$}
head -n -1 header_${chr_name}.tmp > header_meta_${TMPID}.tmp
tail -n 1  header_${chr_name}.tmp > header_cols_${TMPID}.tmp

# Ensure missing INFO fields are declared
grep -q "ID=ReverseComplementedAlleles" header_meta_${TMPID}.tmp || \
echo '##INFO=<ID=ReverseComplementedAlleles,Number=0,Type=Flag,Description="Alleles were reverse complemented during liftover">' >> header_meta_${TMPID}.tmp

grep -q "ID=SwappedAlleles" header_meta_${TMPID}.tmp || \
echo '##INFO=<ID=SwappedAlleles,Number=0,Type=Flag,Description="REF/ALT were swapped during liftover">' >> header_meta_${TMPID}.tmp

# Reassemble clean header
cat header_meta_${TMPID}.tmp header_cols_${TMPID}.tmp | grep -v "bcftools_viewCommand" | grep -v "bcftools_annotateCommand" > header_final_${TMPID}.tmp

# Build final VCF for rsID-mapped variants
cat header_final_${TMPID}.tmp body_all_${chr_name}.tmp | bgzip -c > ${OUT_VCF}.${TMPID}.unsorted

# Use a job-specific tmp dir for bcftools sort
TMPDIR_JOB=${OUT_DIR}/tmp/${TMPID}
mkdir -p $TMPDIR_JOB
bcftools sort -T ${TMPDIR_JOB} ${OUT_VCF}.${TMPID}.unsorted -Oz -o ${OUT_VCF}
tabix -f -p vcf ${OUT_VCF}

# ------------------------
# Cleanup
# ------------------------
rm -f header_*_${TMPID}.tmp body_*_${chr_name}.tmp ${OUT_VCF}.${TMPID}.unsorted
cd "${OUT_DIR}/tmp"
rm -rf "${WORKDIR}"
echo "Finished remapping $IN_VCF → $OUT_VCF"





















