#!/bin/bash
#BSUB -J combine_imputed_vcfs[1-23]
#BSUB -R "rusage[mem=120000]"
#BSUB -o step3_combine_imputed_vcfs_%I.out
#BSUB -e step3_combine_imputed_vcfs_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load plink
module load plink2/alpha6.1amd
module load bcftools
module load htslib

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor
cd ${DIR}/MIS_results

mkdir -p logs
OUTPUT_DIR="./combined_output"
mkdir -p "$OUTPUT_DIR"

TMP_DIR="./imputation_tmp"
mkdir -p "$TMP_DIR"

INPUT_DIRS=("1000G" "HLA" "HRC")
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz
bcftools query -l ${vcf} | awk '{print $1, $1, 0, 0, 1, -9}' > ${DIR}/data/MASTER.fam
bcftools query -l ${vcf} | awk '{print $1, $1, 0, 1, -9}' > ${DIR}/data/PHENO.txt
MASTER_FAM=${DIR}/data/MASTER.fam # FID   IID   FatherID   MotherID   Sex   Phenotype
PHENOTYPE=${DIR}/data/PHENO.txt # FID   IID   Age   Sex   Status

CHROM_LIST=({1..22} X)
chr="${CHROM_LIST[$((LSB_JOBINDEX-1))]}"

echo "Processing chromosome $chr..."

PGEN_BASENAMES=()
R2_SUMMARY_FILES=()

#### -------------------------- VCF -> plink2 ---------------------------- ####
# Step 1: Convert VCFs to plink2 pgen format with conditional allele renaming
for src in "${INPUT_DIRS[@]}"; do
    INPUT_VCF="${src}/chr${chr}.dose.vcf.gz"
    OUT_PREFIX="${src}_chr${chr}"
    if [[ ! -f "$INPUT_VCF" ]]; then
        echo "Warning: missing $INPUT_VCF"
	continue
    fi
    echo "working on chr${chr} in source: ${src}"


    # create plink file for VCF file: handle case of HLA variants
    if [[ "$src" == "VIGOR_HLA" && "$chr" == "6" ]]; then
        echo "creating VCF and plink2 files for HLA variants"
        HLA_VCF="${src}/${src}_chr${chr}_HLA_only.vcf.gz"

	# creat sex and ID map file for PLINK
        SEX_MAP="${TMP_DIR}/${OUT_PREFIX}_sex.map"
        cat $PHENOTYPE | awk '{print 0, $1"_"$2,$4}' > "${SEX_MAP}"
    
        MAP_FILE="${TMP_DIR}/${OUT_PREFIX}_samples.map"
        bcftools query -l "$INPUT_VCF" | awk -F'_' '{
            old_fid = 0;
            old_iid = $0;
            new_fid = $1"_"$2;
            new_iid = "";
            for (i = 3; i <= NF; i++) {
                new_iid = (new_iid == "") ? $i : new_iid"_"$i;
            }
            print old_fid, old_iid, new_fid, new_iid;
        }' > "${MAP_FILE}"


	# reheader to get rid of the weird "0_" in id names
	bcftools query -l "${HLA_VCF}" | sed 's/^0_//' > "${TMP_DIR}/${OUT_PREFIX}_samples.list"

	# Filter to HLA variants, reheader, bgzip + index
        bcftools view -i 'ID ~ "^HLA"' "${INPUT_VCF}" \
          | bcftools reheader -s "${TMP_DIR}/${OUT_PREFIX}_samples.list" \
          | bgzip -c > "${HLA_VCF}"
        
        tabix -f -p vcf "${HLA_VCF}"

	# create plink2 file set
	plink2 --vcf "$HLA_VCF" \
		--update-sex "${SEX_MAP}" \
		--make-pgen \
		--out "${TMP_DIR}/${OUT_PREFIX}.step1"

	plink2 --pfile "${TMP_DIR}/${OUT_PREFIX}.step1" \
		--extract-if-info "R2 >= 0.3" \
		--new-id-max-allele-len 10 missing \
		--update-ids "${MAP_FILE}" \
		--make-pgen --out "${TMP_DIR}/${OUT_PREFIX}"

    else
	echo "creating VCF and plink2 files"
	# reheader VCF to get rid of the weird "0_" in id names
	bcftools query -l ${INPUT_VCF} | sed 's/^0_//' > "${TMP_DIR}/${OUT_PREFIX}_samples.list"
	bcftools reheader -s "${TMP_DIR}/${OUT_PREFIX}_samples.list" ${INPUT_VCF} -o ${INPUT_VCF}.reheader.vcf

        # creat sex and ID map file for PLINK
        SEX_MAP="${TMP_DIR}/${OUT_PREFIX}_sex.map"
	bcftools query -l ${INPUT_VCF}.reheader.vcf | \
		awk '{print 0, $1, 1}' > "${SEX_MAP}"
        
        MAP_FILE="${TMP_DIR}/${OUT_PREFIX}_samples.map"
	bcftools query -l ${INPUT_VCF}.reheader.vcf | \
		awk '{print 0, $1, $1, $1}' > "${MAP_FILE}"
        
	# create plink file for VCF file
        # STEP1: VCF -> PGEN with sex info
        plink2 --vcf ${INPUT_VCF}.reheader.vcf \
                --split-par b37 \
                --make-pgen \
		--update-sex ${SEX_MAP} \
                --out "${TMP_DIR}/${OUT_PREFIX}.step1"


        # STEP2: Apply R2 filtering + rename variants + rename FID/IID
        plink2 --pfile "${TMP_DIR}/${OUT_PREFIX}.step1" \
                --extract-if-info "R2 >= 0.3" \
                --set-all-var-ids '@:#:$r:$a' \
                --new-id-max-allele-len 10 missing \
		--update-ids "${MAP_FILE}" \
                --make-pgen --out "${TMP_DIR}/${OUT_PREFIX}"
    fi    
    PGEN_BASENAMES+=("${TMP_DIR}/${OUT_PREFIX}")

    # extract R2 info from .pvar file; only keep variant if: R2 >= 0.3 and ID is not missing.
    # output columns are: rsID, R2, Source
    awk -v src="$src" '
    BEGIN { OFS = "\t" }
    !/^#/ {
        info = $6
        split(info, fields, ";")
        r2 = -1 
        for (i in fields) {
            if (fields[i] ~ /^R2=/) {
                split(fields[i], kv, "=")
                r2 = kv[2] + 0
            }
        }
        if (r2 >= 0.3 && $3 != ".") {
            print $3, r2, src
        }
    }' "${TMP_DIR}/${OUT_PREFIX}.pvar" > "${TMP_DIR}/${OUT_PREFIX}_variant_r2.tsv"
done


#### -------------------------- Combine per-source R2 summaries & select best variant ---------------------------- ####
# Step 2: Identify the highest R2 imputed variant for each variant across all sources (e.g., 1000G, HRC, HLA)

# Collect all *_variant_r2.tsv files for this chromosome (generated per source)
BEST_R2_FILE="${TMP_DIR}/chr${chr}_best_r2.tsv"
mapfile -t R2_SUMMARY_FILES < <(find "${TMP_DIR}" -name '*_variant_r2.tsv')

# Concatenate all R2 summary files and keep the row with the highest R2 per rsID
cat "${R2_SUMMARY_FILES[@]}" | \
awk '
BEGIN { OFS = "\t" }
{
    key = $1
    if (key in best_r2) {
        if ($2 > best_r2[key]) {
            best_r2[key] = $2
            best_variant[key] = $0
        }
    } else {
        best_r2[key] = $2
        best_variant[key] = $0
    }
}
END {
    for (k in best_variant) {
        print best_variant[k]
    }
}' > "${BEST_R2_FILE}"
sort "${BEST_R2_FILE}" -o "${BEST_R2_FILE}"


#### ------------------------- Step 3: Filter variants + intersect samples ------------------------- ####
# For each imputation source:
#   - Filter variants to those retained in the BEST_R2 list
#   - Identify overlapping samples across all sources (for fair comparison)
#   - Create filtered PLINK BFILES and .pvar files for downstream analysis

FILTERED_BFILES=()
INTERSECT_FILE="${TMP_DIR}/chr${chr}_intersect_samples.txt"
FIRST=true

# First loop: gather variant and sample sets, compute intersection
for src in "${INPUT_DIRS[@]}"; do
    OUT_PREFIX="${src}_chr${chr}"

    # proceed only if PGEN file exists
    if [[ -f "${TMP_DIR}/${OUT_PREFIX}.pgen" ]]; then

	    # extract variants from best-R2 file specific to this source
            VARIANT_LIST="${TMP_DIR}/${OUT_PREFIX}_filtered_variants.txt"
            grep "${src}" "${BEST_R2_FILE}" > "$VARIANT_LIST" 

	    # Extract sample IDs (FID IID) from the .psam
            SAMPLE_FILE="${TMP_DIR}/${OUT_PREFIX}_samples.txt"
            tail -n +2 "${TMP_DIR}/${OUT_PREFIX}.psam" | awk '{print $1,$2}' > "$SAMPLE_FILE"

	    # Compute sample intersection across sources
            if [ "$FIRST" = true ]; then
                cp "$SAMPLE_FILE" "$INTERSECT_FILE"
                FIRST=false
            else
                grep -F -f "$SAMPLE_FILE" "$INTERSECT_FILE" > "${INTERSECT_FILE}.tmp"
                mv "${INTERSECT_FILE}.tmp" "$INTERSECT_FILE"
            fi
    fi
done


# Second loop: create filtered PLINK files using variant + sample intersection
for src in "${INPUT_DIRS[@]}"; do
    OUT_PREFIX="${src}_chr${chr}"

    if [[ -f "${TMP_DIR}/${OUT_PREFIX}.pgen" ]]; then
        VARIANT_LIST="${TMP_DIR}/${OUT_PREFIX}_filtered_variants.txt"

	# Export as full BFILE (bed/bim/bam) filtered by sample + variant
        plink2 --pfile "${TMP_DIR}/${OUT_PREFIX}" \
            --extract "$VARIANT_LIST" \
            --keep "$INTERSECT_FILE" \
            --make-bed --out "${TMP_DIR}/${OUT_PREFIX}_filtered"
        
	# Also output updated .pvar only (can be used for QC or merges)
        plink2 --pfile "${TMP_DIR}/${OUT_PREFIX}" \
            --extract "$VARIANT_LIST" \
            --make-just-pvar --out "${TMP_DIR}/${OUT_PREFIX}_filtered"

	# Tract output BFILES
        FILTERED_BFILES+=("${TMP_DIR}/${OUT_PREFIX}_filtered")
    else
        echo "Warning: missing ${TMP_DIR}/${OUT_PREFIX}.pgen"
    fi
done


#### ----------------------------- Step 4: Merge filtered BFILES across sources ----------------------------- ####
# Goal: Combine the best-R² variants (after filtering) from all sources into a single BFILE set for this chromosome

# Create a merge list of BFILES (excluding the first one, which is the base)
MERGE_LIST="${TMP_DIR}/chr${chr}_merge_list.txt"
> "$MERGE_LIST"  # Clear or create the merge list file

# If chrX, fix PAR1/PAR2 chromosome labels in .bim files to 'X' (PLINK1.9 doesn't support 'PAR1')
if [[ "${chr}" == "X" ]]; then
    for f in "${FILTERED_BFILES[@]}"; do
        BIM="${f}.bim"
        cp "$BIM" "${BIM}.bak"
        awk 'BEGIN{OFS="\t"} { if ($1 == "PAR1" || $1 == "PAR2") $1 = "X"; print }' "${BIM}.bak" > "$BIM"
    done
fi

# Always populate the merge list (regardless of chromosome)
for f in "${FILTERED_BFILES[@]:1}"; do
    echo "$f" >> "$MERGE_LIST"
done

MERGED_PREFIX="${TMP_DIR}/chr${chr}_merged"

# Merge using PLINK1.9 (only PLINK1.9 supports --merge-list)
plink --bfile "${FILTERED_BFILES[0]}" \
      --merge-list "$MERGE_LIST" \
      --make-bed \
      --out "$MERGED_PREFIX"

# After merging, PLINK1.9 erased the sex and phenotype information
# so we re-update it here and convert to PLINK2 format.
MERGED_FAM="${MERGED_PREFIX}.fam"
UPDATED_SEX_FILE="${TMP_DIR}/chr${chr}_sex_update.txt"
UPDATED_PHENO_FILE="${TMP_DIR}/chr${chr}_pheno_update.txt"

awk '{print $1, $2}' "$MERGED_FAM" > "${TMP_DIR}/chr${chr}_samples_to_update.txt"

grep -Ff "${TMP_DIR}/chr${chr}_samples_to_update.txt" "$MASTER_FAM" | awk '{print $1, $2, $5}' > "$UPDATED_SEX_FILE"
grep -Ff "${TMP_DIR}/chr${chr}_samples_to_update.txt" "$MASTER_FAM" | awk '{print $1, $2, $6}' > "$UPDATED_PHENO_FILE"

plink2 \
  --bfile "$MERGED_PREFIX" \
  --split-par b37 \
  --update-sex "$UPDATED_SEX_FILE" \
  --pheno "$UPDATED_PHENO_FILE" \
  --make-pgen \
  --out "${OUTPUT_DIR}/chr${chr}_imputed_combined"



#### ----------------------------- Step 5: Final annotation of merged imputed dataset ----------------------- ####
# 1. Prepare annotation for .pvar variants, annotate variants with R2 and SOURCE
FILTERED_PVAR_FILES=()
for f in "${FILTERED_BFILES[@]}"; do
    FILTERED_PVAR_FILES+=("${f}.pvar")
done

awk -v r2_file="${BEST_R2_FILE}" '
    BEGIN {
        FS = OFS = "\t"
        while ((getline < r2_file) > 0) {
            id = $1
            r2 = $2
            src = $3
            info[id] = "R2=" r2 ";SOURCE=" src
        }
    }
    !/^#/ {
        id = $3
        if (id in info) {
            print $3, info[id]
        }
    }
' "${FILTERED_PVAR_FILES[@]}" | sort -k1,1V -k2,2n | uniq > "${MERGED_PREFIX}.filtered_combined.anno"


# 2. Attach INFO field to .pvar as a new 6th column; if a variant has no annotation, fills with R2=NA; SOURCE=NA
awk '
BEGIN { FS = OFS = "\t" }
NR==FNR { anno[$1] = $2; next }
/^#/ {
    if ($1 == "#CHROM") {
        print $1, $2, $3, $4, $5, "INFO"
    } else {
        print
    }
    next
}
{
    info = ($3 in anno) ? anno[$3] : "R2=NA;SOURCE=NA"
    print $1, $2, $3, $4, $5, info
}
' "${MERGED_PREFIX}.filtered_combined.anno" ${OUTPUT_DIR}/chr${chr}_imputed_combined.pvar > "${OUTPUT_DIR}/chr${chr}_imputed_combined_withanno.pvar"

# 3. Prepand INFO field header lines, makes it VCF-style compatible
pvar_file="${OUTPUT_DIR}/chr${chr}_imputed_combined_withanno.pvar"
{
    echo '##INFO=<ID=R2,Number=1,Type=Float,Description="Imputation INFO score from reference panel">'
    echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Imputation reference panel: HRC, 1KG, MHC_alleles">'
    cat "$pvar_file"
} > tmp && mv tmp "$pvar_file"

#### ----------------------------- Step 6: Create and post-process VCF ----------------------------- ####
# Export the PGEN-format PLINK2 data into a VCF.gz file
FINAL_PLINK_OUT_PREFIX="${TMP_DIR}/chr${chr}_imputed_combined"
plink2 --pfile "${OUTPUT_DIR}/chr${chr}_imputed_combined" \
    --export vcf bgz \
    --out "$FINAL_PLINK_OUT_PREFIX"

CLEANED_VCF="${TMP_DIR}/tmp_chr${chr}_imputed_combined_cleaned.vcf.gz"
FINAL_VCF="${OUTPUT_DIR}/chr${chr}_imputed_combined_hg19.vcf.gz"

# Step 6.1: Remove non-standard symbolic alleles that PLINK might include (e.g., REF = NON_REF)
bcftools view -e 'REF ~ "^<"' "${FINAL_PLINK_OUT_PREFIX}.vcf.gz" -O z -o "$CLEANED_VCF"
bcftools index -f "$CLEANED_VCF"

# Step 6.2: Rename chromosomes to UCSC-style (e.g. 1 -> chr1, X -> chrX)
bcftools annotate \
--rename-chrs <(echo -e "${chr}\tchr${chr}") \
-O z -o "$FINAL_VCF" "$CLEANED_VCF"
bcftools index -f "$FINAL_VCF"


##### ------------------------------ variant count summary file -------------------------------- ####
#./imputation_tmp/VIGOR_HRC_chrX.step1.log
#for src in "${INPUT_DIRS[@]}"; do
#	n1=$(cat ${TMP_DIR}/${src}_chr${chr}.log | grep "variants loaded from" | cut -d' ' -f4)
#	n2=$(cat ${TMP_DIR}/${src}_chr${chr}.log | grep "variants loaded from" | cut -d' ' -f1)
#	n3=$(cat ${TMP_DIR}/${src}_chr${chr}_filtered.log | grep "variants remaining after main filters"  | cut -d' ' -f1)
#	echo -e "MIS result: ${n1}, after R2 filter: ${n2}, after ${n3} ${src}" >> $TMP_DIR/variant_count_summary.txt
#done
#cat ${TMP_DIR}/chr${chr}_merged.bim | grep "variants loaded from"
#cat ${TMP_DIR}/chr${chr}_imputed_combined.log | grep "variants loaded from"

###### -------------------------------- chrX non-PAR variant haploidity sanity check ------------- ####
## in PAR1 and PAR2: should all be diploid, heterozygous is fine
## in non-PAR, females (XX) diploid okay, 0/0, 0/1, 1/1; males (XY) haploid only. should be all single allele genotype
#awk 'NR>1 { combined=$1"_"$2; print combined, $3 }' chrX_imputed_combined.psam > sex.map
#bcftools query -r chrX:2700157 -f '[%SAMPLE\t%GT\n]' chrX_imputed_combined_hg19.vcf.gz | awk 'FNR==NR{sex[$1]=$2; next} { gt=$2; hap=(gt ~ /^(\.|[01])$/); print $1, gt, "SEX="sex[$1], (hap?"HAP":"DIP") }' sex.map - | column -t
## passed the sanity check
