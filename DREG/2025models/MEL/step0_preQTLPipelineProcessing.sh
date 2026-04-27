#!/bin/bash
# shuo.shan@umassmed.edu


# ----------------------------------------- #
# -------- Global Variables        -------- #
# ----------------------------------------- #
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL
metagene=/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt
genemodel=/pi/manuel.garber-umw/human/skin/eQTLs/literature/gencode/gencode_v34_allExpressedGenes_TSS.bed
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz

mkdir -p ${DIR}/data
mkdir -p ${DIR}/logs
mkdir -p ${DIR}/results
cd ${DIR}
# ----------------------------------------- #
# -------- Make SNP-gene pair list -------- #
# ----------------------------------------- #
module load bcftools
module load bedtools
module load htslib

# Convert VCF → BED format (1-based → 0-based)
bcftools view -H ${vcf} | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' > ${DIR}/data/SNPs.bed

# get every SNP within +/- 500kb of a gene's TSS
pairs=${DIR}/data/SNPs_near_TSS.bed
bedtools window -w 500000 -a ${genemodel} -b ${DIR}/data/SNPs.bed > ${pairs}

# add header
echo -e "gene_chr\tgene_start\tgene_end\tgene_name\t0\tstrand\tgene_id\tSNP_chr\tSNP_start\tSNP_end\tSNP_ID" \
  | cat - ${pairs} > ${pairs}.tmp && mv ${pairs}.tmp ${pairs}

# split into chunks
echo "Splitting into 200k-line chunks..."
mkdir -p ${DIR}/chunks
header=$(head -n1 ${pairs})
tail -n +2 ${pairs} | split --numeric-suffixes=0 --suffix-length=3 -l 200000 - ${DIR}/chunks/pairs_chunk_
# Re-add header to each chunk
for f in ${DIR}/chunks/pairs_chunk_*; do
  mv "$f" "$f.tsv"
  { echo "$header"; cat "$f.tsv"; } > "${f}.tmp" && mv "${f}.tmp" "${f}.tsv"
  rm -f "$f"
done

echo "Done. Chunks in ${DIR}/chunks/"


# detect total number of chunks automatically
NCHUNKS=$(ls ${DIR}/chunks/pairs_chunk_*.tsv | wc -l)
echo "Total chunks: ${NCHUNKS}"

# generate the right job-array
echo "#BSUB -J QTL[1-${NCHUNKS}]%70"
# #BSUB -J QTL[1-135]%70 # create 135 jobs but run only 70 at a time concurrently, queue the rest
