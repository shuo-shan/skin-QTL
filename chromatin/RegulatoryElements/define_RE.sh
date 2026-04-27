#!/bin/bash
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# script overview: define promoters and enhancers based on ATACseq peaks
# method: ELisz 2018 paper. For details see Crystal Dropbox log on promoter_enhancer_annotation.docx

module load bedtools/2.30.0
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements
cd ${dir}

##### step 1. merge ATACseq peaks that are within 200bp of each other
krt=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/KRT_merged/report5430/chip/merged.bed
mel=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/MEL/report5356/chip/merged.bed
frb=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/merged.bed

geneExprs_krt=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/analysis/expressedGenes.txt
geneExprs_frb=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/analysis/expressedGenes.txt
geneExprs_mel=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/analysis/expressedGenes.txt

atac_krt=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_KRT_merged.bed
atac_mel=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_MEL_merged.bed
atac_frb=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_FRB_merged.bed

expressedGenes=${geneExprs_krt}
atac=${atac_krt}
h3k27acPeaks=${krt} # pick a celltype to focus on.
celltype=KRT # change this, too!

##### step 2. calculate peak distance to nearest TSS
# DEFINE TSS from Ensembl GTF file, by /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/log.sh
# columns are: TSS position, strand, Ensembl gene ID + transcript ID, Gene name + transcript ID, gene name, gene biotype
# TSS bed file is created from the start of all transcripts. script: /pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/ensembl2tss.sh
tss=/pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_transcription_start_sites.bed

# annotate each ATACseq peak with proximity to closest TSS
# columns: RE chromosome, RE start, RE end (half open, half closed), RE id, TSS chromosome, TSS start, TSS end, TSS strand, TSS_ensembl_gene_id, TSS_Ensembl_geneName, TSS_Ensembl_gene_biotype, distance 
bedtools closest -a ${atac} -b ${tss} -d > ${celltype}_with_TSS_distance.txt 

##### step 3. define promoters to be ATACseq peak < 500bp of TSS
# promoter name: promoter_geneName_promoterLocation_distance-to-TSS
cat ${celltype}_with_TSS_distance.txt | awk '{OFS="\t"}{if ($9 < 500) print $0,"promoter_"$8"_"$4"_"$9}' > ${celltype}_promoters.txt
##### step 4. define enhancers to be H3K27ac peaks with closest expressed gene within 300kbp of TSS
cat ${celltype}_with_TSS_distance.txt | awk '{if ($9 >= 500) print $0}' > temp.not_promoters.txt
# only pick TSSs of genes that are expressed from our RNA-seq data (avg > 10 cpm in either PBS or IFN across samples)
bash /pi/manuel.garber-umw/sshan/scripts/function_rapid_fgrep.sh ${expressedGenes} ${tss} TSS_of_expressed_genes.txt rapidGrep 500
bedtools closest -a temp.not_promoters.txt -b ${dir}/TSS_of_expressed_genes.txt -d | awk '{OFS="\t"}{if ($9 < 300000) print $0,"enhancer_"$8"_"$4"_"$9}' > ${celltype}_enhancers.txt
rm temp.not_promoters.txt
rm ${celltype}_with_TSS_distance.txt
rm -r rapid_fgrep_temp

##### step 5. annotate promoters/enhancers to be open/closed based on ATAC-seq overlap.
# annotate promoters
cat ${celltype}_promoters.txt | cut -f1,2,3,10 | bedtools intersect -a stdin -b ${atac} -u | awk '{OFS="\t"}{print $0,"open"}' > temp.promoters.open.txt
cat ${celltype}_promoters.txt | cut -f1,2,3,10 | bedtools intersect -a stdin -b ${atac} -v | awk '{OFS="\t"}{print $0,"closed"}' > temp.promoters.closed.txt
cat temp.promoters.open.txt temp.promoters.closed.txt | bedtools sort -i stdin > ${celltype}_promoters_annotated.txt
rm temp.promoters.open.txt temp.promoters.closed.txt ${celltype}_promoters.txt

# annotate enhancers with ATACseq info (closed/accessible)
cat ${celltype}_enhancers.txt | cut -f1,2,3,10 | bedtools intersect -a stdin -b ${atac} -u | awk '{OFS="\t"}{print $0,"open"}' > temp.enhancers.open.txt
cat ${celltype}_enhancers.txt | cut -f1,2,3,10 | bedtools intersect -a stdin -b ${atac} -v | awk '{OFS="\t"}{print $0,"closed"}' > temp.enhancers.closed.txt
cat temp.enhancers.open.txt temp.enhancers.closed.txt | bedtools sort -i stdin > ${celltype}_enhancers_annotated.txt
rm temp.enhancers.open.txt temp.enhancers.closed.txt

##### step 6. annotate promoters/enhancers to merge into one table for SPRITE usage
cat ${celltype}_promoters_annotated.txt | sed 's/_/\t/g' | awk '{OFS="\t"}{print $1,$2,$3,$1";"$2";"$3";"$5";"$4}' >> ${celltype}_k27acPeaks.bed
cat ${celltype}_enhancers_annotated.txt | awk '{OFS="\t"}{print $1,$2,$3,$1";"$2";"$3";enhancer;enhancer"}' >> ${celltype}_k27acPeaks.bed
bedtools sort -i ${celltype}_k27acPeaks.bed | awk '!seen[$4]++' > ${celltype}_k27acPeaks_sorted.bed
rm ${celltype}_k27acPeaks.bed

##### step 7. merge all promoters across 3 cts, and all enhancers across 3 cts
cat *_promoters_annotated.txt | bedtools sort -i stdin | bedtools merge -i stdin | awk '{OFS="\t"}{print $1,$2,$3,"promoter_"$1"_"$2"_"$3}' | bedtools sort -i stdin > allcts_promoters.bed
cat *_enhancers_annotated.txt | bedtools sort -i stdin | bedtools merge -i stdin | awk '{OFS="\t"}{print $1,$2,$3,"enhancer_"$1"_"$2"_"$3}' | bedtools sort -i stdin > allcts_enhancers.bed

##### step 8. compute ATACseq genome coverage of 200bp window flanking the center of each promoter/enhancer
bash multicov.sh

##### step 9. link promoter/enhancer to a gene
bedtools closest -a allcts_enhancers.bed -b ${dir}/TSS_of_expressed_genes.txt -d > allcts_enhancers_with_tss_distance.bed
bedtools closest -a allcts_promoters.bed -b ${dir}/TSS_of_expressed_genes.txt -d > allcts_promoters_with_tss_distance.bed






