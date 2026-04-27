#!/bin/bash
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# script overview: define promoters and enhancers based on H3K27ac peaks
# method: ELisz 2018 paper. For details see Crystal Dropbox log on promoter_enhancer_annotation.docx

module load bedtools/2.30.0
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method4
cd ${dir}

##### step 1. call summits in ATACseq peaks for all merged bam files, extend by 300bp.
summit=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits.bed
genomeFile=/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.chrom.sizes
atac_window=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_300bp_flanking_window.bed
bedtools slop -b 150 -i ${summit} -g ${genomeFile} | sort -k 1,1 -k2,2n > ${atac_window}

##### step2. define krt, mel, and frb enhancers to be ATACseq peaks that overlap with a chipseq peak in PBS or IFN for that celltype
chip_krt=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/KRT_merged/report5430/chip/merged.bed
chip_mel=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/MEL/report5356/chip/MEL_merged.bed
chip_frb=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/FRB_merged.bed
chip_krt_pbs=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/KRT_merged/report5430/chip/PBS_merged.bed
chip_krt_ifn=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/KRT_merged/report5430/chip/IFN_merged.bed
chip_mel_pbs=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/MEL/report5356/chip/PBS_merged.bed
chip_mel_ifn=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/MEL/report5356/chip/IFN_merged.bed
chip_frb_pbs=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/PBS_merged.bed
chip_frb_ifn=/pi/manuel.garber-umw/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/IFN_merged.bed

atac=${atac_window}
expressedGenes=${geneExprs_frb}
h3k27acPeaks=${chip_frb} # pick a celltype to focus on.
celltype=FRB # change this, too!

bedtools intersect -u -a ${atac_window} -b ${chip_krt} > temp_cRE_KRT.bed
bedtools intersect -u -a ${atac_window} -b ${chip_frb} > temp_cRE_FRB.bed
bedtools intersect -u -a ${atac_window} -b ${chip_mel} > temp_cRE_MEL.bed

##### step 3. take all cREs within 300kbp of a gene's TSS.
# this step is ran with this script: command_calling_function_fetch_cRE_surround_gene.sh
bsub -J fetchcRE < command_calling_function_fetch_cRE_surround_gene.sh

##### step 4. annotate ATAC-seq peaks with the H3K27ac activity in each celltype and condition. filter out peaks that are inactive in any ct and condition.
# annotate active/inactive atac peaks by overlapping KRT_PBS, KRT_IFN, MEL_PBS, etc
cat ${atac_window} | bedtools intersect -a stdin -b ${chip_krt_pbs} -u | awk '{OFS="\t"}{print $4}' > temp_atac_active_in_krt_pbs.txt
cat ${atac_window} | bedtools intersect -a stdin -b ${chip_krt_ifn} -u | awk '{OFS="\t"}{print $4}' > temp_atac_active_in_krt_ifn.txt
cat ${atac_window} | bedtools intersect -a stdin -b ${chip_mel_pbs} -u | awk '{OFS="\t"}{print $4}' > temp_atac_active_in_mel_pbs.txt
cat ${atac_window} | bedtools intersect -a stdin -b ${chip_mel_ifn} -u | awk '{OFS="\t"}{print $4}' > temp_atac_active_in_mel_ifn.txt
cat ${atac_window} | bedtools intersect -a stdin -b ${chip_frb_pbs} -u | awk '{OFS="\t"}{print $4}' > temp_atac_active_in_frb_pbs.txt
cat ${atac_window} | bedtools intersect -a stdin -b ${chip_frb_ifn} -u | awk '{OFS="\t"}{print $4}' > temp_atac_active_in_frb_ifn.txt
# next I want to create a dictionary of all atac peaks, if it's overlapping krt_pbs chip peak, annotate active in its column. create 6 columns for each ctxcondition.
# use R scirpt here: Rscript_annotate_H3K27ac_activity.R, R environment data here: myEnvironment_annotate_H3K27ac_activity_for_atac.RData
# this generates a dictionary_atac_peaks_annotated_H3K27ac_activity.txt file. around 100K atac-seq windows.


##### step 5. annotate ATAC-seq peaks with the neighboring genes upstream and downstream w.r.t. the peak.
# then I run this script: /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/command_join_table.sh 
# ^ this script produces three tables, one for each celltype. In each table the first column is the ATACseq peak name, second column is the closest gene, third column is the closest gene in the other direction. Table name is: FRB_regulatory_elements.txt
# table name is: MEL_atac_peak_and_neighboring_genes.txt FRB_atac_peak_and_neighboring_genes.txt KRT_atac_peak_and_neighboring_genes.txt

##### step 6. create gene-promoter dictionary and gene-enhancer dictionary (negative distance means gene tss is upstream of promoter.)
# header: gene, promoter name, distance, celltype
# ended up not using------> script: command_link_gene_cRE.sh. in this version, both expressed genes up and downstream of the enhancer are kept.
# this is the one I kept------> script: command_link_gene_cRE_closest.sh. in this version, only closest expressed gene to the enhancer is kept.
# ------> script: command_dictionary2bed.R. turn dictionary to bed file for IGV viewing.

##### step 7. calculate the genome coverage of all filtered ATACseq peaks in each ATACseq and ChIP-seq bam file.
# ------> script: command_multicov.sh


##### DIDN"T DO THIS: step 8. annotate promoters/enhancers to merge into one table for SPRITE usage
rm ${celltype}_atacPeaks.bed
cat ${celltype}_promoters.txt | sed 's/_/\t/g' | awk '{OFS="\t"}{print $1,$2,$3,$1";"$2";"$3";"$5";"$4}' >> ${celltype}_atacPeaks.bed
cat ${celltype}_enhancers.txt | awk '{OFS="\t"}{print $1,$2,$3,$1";"$2";"$3";enhancer;enhancer"}' >> ${celltype}_atacPeaks.bed
bedtools sort -i ${celltype}_atacPeaks.bed | awk '!seen[$4]++' > ${celltype}_atacPeaks_sorted.bed
rm ${celltype}_atacPeaks.bed





