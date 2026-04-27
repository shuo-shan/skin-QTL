#!/bin/bash
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=5000]"
#BSUB -W 72:00
#BSUB -q long
#BSUB -J chunk[1-120]%120 # run 120 jobs concurrently
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/new_round1_adaptive_permut_%J_%I.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/MEL/logs/new_round1_adaptive_permut_%J_%I.err"

sleep $((RANDOM % 20)) # reduces heavy I/O thundering-herd effect in parallel jobs
module load bcftools
ct=MEL
DIR=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
mkdir -p ${DIR}/permutation
mkdir -p ${DIR}/permutation/chunk
mkdir -p ${DIR}/permutation/model_stats

all_genes=${DIR}/results_QC/all_master_pairs_${ct}_genes.txt
all_pairs=${DIR}/results_QC/all_master_pairs_${ct}.txt
all_modelstats=${DIR}/permutation/model_stats/all_model_results.txt
vcf="/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute/m2_rsIDanchor/MIS_results/combined_output/final_output/skineQTL_imputed_hg38_filtered.vcf.gz"

# Convert LSB_JOBINDEX (1–120) → chunk_id (000–119)
chunk_id=$(printf "%03d" $((LSB_JOBINDEX - 1)))

# ----------------------------------------------------------------------  #
#### summarize all SNP gene pairs with any QTL statistics for this celltype
#cd ${DIR}/results_QC
#for chunk in $(seq -w 0 134); do
#	echo ${chunk}
#  cat modeling_stats_postQC_*_"${chunk}".txt \
#    | grep -v 'snp' \
#    | awk '{print $1"_"$2}' \
#    | sort -u \
#    | tr '_' '\t' \
#    > "master_pairs_${ct}_${chunk}.txt"
#done
#
#cat master_pairs_${ct}*.txt | sort -t $'\t' -k2,2 -k1,1 > ${all_pairs}
#cat all_master_pairs_${ct}.txt | cut -f2 | sort -u > ${all_genes}
#rm master_pairs_${ct}*.txt
#
## ----------------------------------------------------------------------  #
## concatenate all model stats into one file
#rm -f ${DIR}/permutation/model_stats/model_results.txt
#cd ${DIR}/results
#awk 'NR==1' result_000.tsv > ${DIR}/permutation/model_stats/all_model_results.txt
#for f in result_*.tsv;do
#	echo ${f}
#	awk 'NR>1' ${f} >> ${DIR}/permutation/model_stats/all_model_results.txt
#done


# ----------------------------------------------------------------------  #
#### chunk gene_SNP pairs for every 100 genes
#cd ${DIR}/permutation/chunk
#rm gene_chunk_*
#NGENES=100
#split --numeric-suffixes=0 --suffix-length=3 -l ${NGENES} ${all_genes} ${DIR}/permutation/chunk/gene_chunk_
#echo "Num gene chunks: $(ls -1 gene_chunk_* | wc -l)"
#
#cd ${DIR}/permutation/chunk
#f=gene_chunk_${chunk_id}
#out_pair=pair_chunk_${chunk_id}.txt
#out_genotype=genotype_pair_chunk_${chunk_id}.txt
#out_modelstats=${DIR}/permutation/model_stats/model_stats_${chunk_id}.txt
#echo "Chunk ${chunk_id}: starting now"; date
#
### ---------------------------------------------------------- #
#### Load genes from gf into a hash; keep rows in all_pairs whose gene (col2) is in that hash
#awk -v FS="\t" -v OFS="\t" '
#NR==FNR {g[$1]=1; next}
#($2 in g) {print $1, $2}' ${f} ${all_pairs} > ${out_pair}


## ---------------------------------------------------------- #
#### Extract genotype of SNPs in this pair chunk file
#snplist=temp.${chunk_id}.snps.txt
#vcf_subset=temp.${chunk_id}.vcf.gz
#samples=temp.${chunk_id}.samples.txt
#header=temp.${chunk_id}.header.txt
#body=temp.${chunk_id}.body.txt
#geno=temp.${chunk_id}.genotype.txt

## 1) collect SNP IDs from this chunk
#cut -f1 ${out_pair} | sort -u > ${snplist}
#
## 2) extract VCF records
#bcftools view --threads 5 --include ID==@${snplist} ${vcf} -Oz -o ${vcf_subset}
#
## 3) clean sample names
#bcftools query -l ${vcf_subset} | cut -d'_' -f1 | sed 's/skineQTL-//' | sed 's/F0/F/' > ${samples}
#
## 4) build header
#printf "CHROM\tPOS\tID\tREF\tALT\t" > ${header}
#paste -sd '\t' ${samples} >> ${header}
#
## 5) extract GT matrix
#bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' ${vcf_subset} > ${body}
#cat ${header} ${body} > ${geno}
#
## 6) convert genotypes to numeric
#sed -E 's/0\/0/0/g; s/0\/1/1/g; s/1\/0/1/g; s/1\/1/2/g; s/\.\/\./NA/g' ${geno} > ${out_genotype}
#
## 7) cleanup intermediate files (leave only numeric genotype)
#rm "$snplist" "$vcf_subset" "$samples" "$header" "$body" "$geno"
#
## ---------------------------------------------------------- #
#### extract modeling results for the gene SNP pair
#rm -f ${out_modelstats}
#
#head -1 ${all_modelstats} > ${out_modelstats}
#
#cd ${DIR}/permutation/chunk
#awk -v FS="\t" -v OFS="\t" '
#NR==FNR {g[$1]=1; next}
#($2 in g) {print $0}' ${f} ${all_modelstats} >> ${out_modelstats}
#
## ---------------------------------------------------------- #
# Submit chunk to adaptive permutation Rscript
echo "Performing Round1 adaptive permutation for ${chunk_id}";date
cd ${DIR}

# Tell R to use my library inside the container:
export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.2"

singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_adaptive_perm_gene_eQTL_perGeneChunk.R MEL PBS eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_adaptive_perm_gene_eQTL_perGeneChunk.R MEL IFNG eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_adaptive_perm_gene_eQTL_perGeneChunk.R MEL IFNB eQTL ${chunk_id}
singularity exec /share/pkg/containers/rstudio_example/r_4.2.2.sif \
	Rscript ${DIR}/step4_adaptive_perm_gene_eQTL_perGeneChunk.R MEL TNF eQTL ${chunk_id}
                                                             








