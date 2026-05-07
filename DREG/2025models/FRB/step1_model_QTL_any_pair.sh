#!/bin/bash
# make a PDF for gene snp pair of interest given the SNP is in the gene's cis-window.
# test QTL modeling result for any SNP gene pair of interest

module load bcftools

ct=FRB

g=$1   #g=IL1B
snp=$2 #snp=rs4790797
cond=$3
QTLtype=$4

## for testing only
#g=RASIP1
#snp=rs2287921
#cond=PBS
#QTLtype=eQTL

DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/results/temp
mkdir -p ${DIR}/results/temp/${g}_${snp}
cd ${DIR}/results/temp/${g}_${snp}

echo "starting the script for ${g}:${snp}"; date

# ------------------------------------
# specify here for nPEER and nGPC
# ------------------------------------
nPEER=10
nGPC=2

# ---------- get gene, SNP, metadata, and modeling stats info ----------- #
VST_FILE=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/VST.sampleFiltered.metaConverted.txt
head -1 ${VST_FILE} > vst.txt
awk -v gene=${g} '{OFS=FS="\t"}{if ($1==gene) print $0}' ${VST_FILE} >> vst.txt
echo "got gene VST"

# get SNP genotype
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz
vcf_subset=snp.vcf.gz
samples=samples.txt
header=header
body=body
geno=temp.snp.txt
geno_num=genotype.txt

echo "==== Subsetting genotype VCF for ${snp} ===="; date
# 1) collect SNP IDs from this chunk
echo "${snp}" > snplist.txt

# 2) extract VCF records
bcftools view --include ID==@snplist.txt ${vcf} -Oz -o ${vcf_subset}

# 3) clean sample names
bcftools query -l "$vcf_subset" | cut -d'_' -f1 | sed 's/F0/F/g' | sed 's/skineQTL-//g' > "$samples"

# 4) build header
printf "CHROM\tPOS\tID\tREF\tALT\t" > "$header"
paste -sd '\t' "$samples" >> "$header"

# 5) extract GT matrix
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf_subset" > "$body"
cat "$header" "$body" > "$geno"

# 6) convert genotypes to numeric
sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' "$geno" > "$geno_num"

# 7) cleanup intermediate files (leave only numeric genotype)
rm snplist.txt "$vcf_subset" "$samples" "$header" "$body" "$geno"

echo "got SNP genotype"


# ------------ call Rscript to plot --------------- #
echo "modeling SNP gene pair for QTL association now..."
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"
singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
	Rscript /pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/step1_model_QTL_any_pair.R ${g} ${snp} ${nPEER} ${nGPC} ${ct} ${cond} ${QTLtype}

echo "done!"; date










