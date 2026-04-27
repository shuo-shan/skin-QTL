#!/bin/bash
# written by shuo.shan@umassmed.edu, 03/2024

dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/QTL_region_TF_enrichment
cd ${dir}

#fetch and prune list of melanocyte reQTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk 'NR>1{if ($7 != "." && $7 < 0.000001) print $1}' | sort | uniq > reQTL.txt

genotype_table=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
module load bcftools
bcftools view --include ID==@reQTL.txt ${genotype_table} -Oz -o reQTL.vcf.gz
bcftools +prune -m 0.2 -w 1000 reQTL.vcf.gz -Oz -o pruned_reQTL.vcf.gz
bcftools index pruned_reQTL.vcf.gz
bcftools view -H pruned_reQTL.vcf.gz | cut -f3 > pruned_reQTL.txt

#fetch 100bp region centered on the reQTL
bed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/master_filtered_genotype.bed
grep -w -f pruned_reQTL.txt ${bed} > pruned_reQTL.bed
module load bedtools
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
bedtools slop -b 50 -i pruned_reQTL.bed -g ${genome} > pruned_reQTL_100bp.bed

#TF motif enrichment analysis using XSTREME
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme
export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH
bedtools getfasta -fi /share/GHPCC/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa -bed pruned_reQTL_100bp.bed > temp.fa
cat temp.fa | tr ':' '-' | tr '-' '_' > pruned_reQTL_100bp.fa
rm temp.fa
xstreme --mea-only --seed 42 --p pruned_reQTL_100bp.fa --m ${memeF}

# clean-up
mv xstreme.html ../
rm -r xstreme_out/
