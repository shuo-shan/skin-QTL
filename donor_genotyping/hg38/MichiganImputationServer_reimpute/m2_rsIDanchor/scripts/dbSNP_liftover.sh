#!/bin/bash
# use liftover to liftover dbSNP variants from hg38 to hg19, then compare with dbSNP hg19 positions
# author: shuo.shan@umassmed.edu
#BSUB -J liftover_test[1-2]
#BSUB -R "rusage[mem=120000]"
#BSUB -o liftoverTest_%I.out
#BSUB -e liftoverTest_%I.err
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
LIFTOVER="/pi/manuel.garber-umw/human/VIGOR/crystal/bin/liftOver"        # UCSC liftOver binary
CHAIN="/pi/manuel.garber-umw/human/VIGOR/crystal/data/hg38ToHg19.over.chain.gz"    # UCSC chain
HG19_FA="/pi/manuel.garber-umw/human/VIGOR/crystal/data/human_g1k_v37.fasta" # hg19/GRCh37 fasta (no 'chr' in names)
mkdir -p ${DATA_DIR}/liftover_testing
cd ${DATA_DIR}/liftover_testing

# --------------------------------------------------- #
# ---- STEP 1. get hg19 rsID positions ---- #
# --------------------------------------------------- #
DBSNP_VCF_hg19=${DATA_DIR}/dbSNP_hg19/dbSNP_hg19_${chr_name}_renamed.vcf.gz
IDMAP_TXT_hg19=${DATA_DIR}/dbSNP_hg19/idmap_hg19_${chr_name}.txt
DBSNP_VCF_hg38=${DATA_DIR}/dbSNP_hg38/dbSNP_${chr_name}_hg38_renamed.vcf.gz
IDMAP_TXT_hg38=${DATA_DIR}/dbSNP_hg38/idmap_${chr_name}_hg38.txt

cat $IDMAP_TXT_hg19 | tr ':' '\t' | awk '{OFS="\t"}{print $5,"chr"$1,$2}' | sort -k1,1 > temp_hg19_${chr_name}.txt
cat $IDMAP_TXT_hg38 | awk '{OFS="\t"}{print $1,$2,$3}' | sort -k1,1 > temp_hg38_${chr_name}.txt

echo -e "#rsID\thg19_chr\thg19_pos\thg38_chr\thg38_pos" > merged_${chr_name}_hg19_hg38.txt
join -t $'\t' -1 1 -2 1 temp_hg19_${chr_name}.txt temp_hg38_${chr_name}.txt \
| awk '{OFS="\t"}{print $1, $2, $3, $4, $5}' \
| awk '!seen[$0]++' >> merged_${chr_name}_hg19_hg38.txt


# --------------------------------------------------- #
# ---- STEP 2. liftover hg38 --> hg19 ---- #
# --------------------------------------------------- #
# make hg38 bed file
awk 'NR>1 {OFS="\t"; print $4, $5-1, $5, $1}' merged_${chr_name}_hg19_hg38.txt > hg38_coords_${chr_name}.bed

# run UCSC liftOver
${LIFTOVER} hg38_coords_${chr_name}.bed ${CHAIN} lifted_hg19_${chr_name}.bed unmapped_${chr_name}.bed


# ----------------------------------------------------------------- #
# ---- STEP 3. compare liftOver hg19 coord vs. dbSNP hg19 coord---- #
# ----------------------------------------------------------------- #

# compile joined list
cat unmapped_${chr_name}.bed | grep -v "Deleted" | awk '{OFS="\t"}{print "NA","NA","NA",$4}' > failed_LiftOver_${chr_name}.txt
cat lifted_hg19_${chr_name}.bed failed_LiftOver_${chr_name}.txt > lifted_hg19_all_${chr_name}.txt

head -n1 merged_${chr_name}_hg19_hg38.txt \
  | awk '{print $0 "\tliftover_hg19_chr\tliftover_hg19_pos"}' \
  > merged_${chr_name}_hg19_hg38_with_liftover.txt

awk 'NR==FNR{lift_chr[$4]=$1; lift_pos[$4]=$3; next}
     NR>1 && $1 !~ /^#/{
         OFS="\t";
         chr = ($1 in lift_chr) ? lift_chr[$1] : "NA";
         pos = ($1 in lift_pos) ? lift_pos[$1] : "NA";
         print $0, chr, pos;
     }' \
     lifted_hg19_all_${chr_name}.txt merged_${chr_name}_hg19_hg38.txt \
     >> merged_${chr_name}_hg19_hg38_with_liftover.txt

# summarize
#echo "done, summarizing"; date
#N_total=$(cat merged_${chr_name}_hg19_hg38_with_liftover.txt | grep -v "#" | cut -f1 | sort | uniq | wc -l | cut -d' ' -f1)
#N_failedLift=$(cat merged_${chr_name}_hg19_hg38_with_liftover.txt | grep -v "#" | grep "NA" | cut -f1 | sort | uniq | wc -l | cut -d' ' -f1)
#N_differentChr=$(cat merged_${chr_name}_hg19_hg38_with_liftover.txt | grep -v "#" | grep -v "NA" | awk '{if ($2!=$6) print $1}' | sort | uniq | wc -l | cut -d' ' -f1)
#N_sameChr_differentPos=$(cat merged_${chr_name}_hg19_hg38_with_liftover.txt | grep -v "#" | grep -v "NA" | awk '{if ($2==$6) print $0}' | awk '{if ($3!=$7) print $1}' | sort | uniq | wc -l | cut -d' ' -f1)
#echo -e "${N_total}\t${N_failedLift}\t${N_differentChr}\t${N_sameChr_differentPos}" > liftOver_test_summary_${chr_name}.txt

echo "done, summarizing"; date

echo -e "N_total\tN_failedLift\tN_differentChr\tN_sameChr_differentPos" \
    > liftOver_test_summary_${chr_name}.txt

awk -v OFS='\t' '
    BEGIN {N_total=0; N_failedLift=0; N_diffChr=0; N_sameChr_diffPos=0}
    NR>1 && $1!~/#/ {
        N_total++
        if ($6=="NA" || $7=="NA") {
            N_failedLift++
        } else if ($2!=$6) {
            N_diffChr++
        } else if ($2==$6 && $3!=$7) {
            N_sameChr_diffPos++
        }
    }
    END {
        print N_total, N_failedLift, N_diffChr, N_sameChr_diffPos
    }
' merged_${chr_name}_hg19_hg38_with_liftover.txt >> liftOver_test_summary_${chr_name}.txt

date

# clean-up
#rm temp_hg19_${chr_name}.txt temp_hg38_${chr_name}.txt merged_${chr_name}_hg19_hg38.txt 
#rm hg38_coords_${chr_name}.bed failed_LiftOver_${chr_name}.txt
#rm lifted_hg19_${chr_name}.bed unmapped_${chr_name}.bed lifted_hg19_all_${chr_name}.txt
#rm merged_${chr_name}_hg19_hg38_with_liftover.txt

# archive
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
## 2.3) Build rsID -> hg19 coordinate mapping file
## Compile map file (3 columns: rsID chr pos)
#echo "Building rsID→hg19 mapping for ${chr_name}"; date
#bcftools view -H ${DBSNP_VCF} \
#       | cut -f1-3 \
#       | awk '{OFS=FS="\t"}{print $3,$1,$2}' \
#       | sed 's/chr//g' \
#       > ${IDMAP_TXT}
#echo "Done: dbSNP_hg19_${chr_name}_renamed.vcf.gz + ${IDMAP_TXT}"

## ======= download hg19 fasta and liftover tools =======
#cd ${dir}/data
## Download (GRCh37 / b37 without 'chr' prefix)
#wget -O human_g1k_v37.fasta.gz \
#  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
#gunzip human_g1k_v37.fasta.gz
#REF37=${dir}/data/human_g1k_v37.fasta
#
## Index for tools
#samtools faidx ${REF37}
#
## Create sequence dictionary (.dict) with GATK inside the container
#bsub -q interactive -Is singularity exec /share/pkg/containers/gatk/4.6.0.0/gatk-4.6.0.0.sif \
#  gatk CreateSequenceDictionary -R $REF37 -O ${REF37%.fasta}.dict
#
## get hg38tohg19 liftover chain
#wget -O hg38ToHg19.over.chain.gz \
#  'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
#
## get UCSC liftover binary
#cd ${dir}/scripts
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/liftOver
#chmod +x liftOver
## =================================
