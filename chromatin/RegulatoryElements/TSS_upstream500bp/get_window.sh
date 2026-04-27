#!/bin/bash
# written by crystal shan 02/2024
# get 500bp upstream to TSS of every gene and link it to gene

module load bedtools/2.30.0
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/TSS_upstream500bp
cd ${dir}

# TSS
# DEFINE TSS from Ensembl GTF file, by /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/log.sh
# columns are: TSS position, strand, Ensembl gene ID + transcript ID, Gene name + transcript ID, gene name, gene biotype
# TSS bed file is created from the start of all transcripts. script: /pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/ensembl2tss.sh
tss=/pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_transcription_start_sites.bed
cat ${tss} | awk '$4="+"' > temp_plusstrand.bed
cat ${tss} | awk '$4="-"' > temp_minusstrand.bed
cat temp_plusstrand.bed | awk '{OFS="\t"}{print $1,$2-500,$3,$4,$5}' > temp_plusstrand_window.bed
cat temp_minusstrand.bed | awk '{OFS="\t"}{print $1,$2,$3+500,$4,$5}' > temp_minusstrand_window.bed
cat temp_plusstrand_window.bed temp_minusstrand_window.bed > temp_TSS_upstream500window.bed

# change format as promoter bed file
cat temp_TSS_upstream500window.bed | awk '{OFS="\t"}{print $1,$2,$3,"window_"NR"_"$5,"TSS_upstream_500bp_window_"NR,$5,"celltype-agnostic","promoter"}' > TSS_upstream500window_alltranscripts.bed
rm temp_*

#### CONCERN!!! #####
# If some genes have multiple transcripts whose TSSs overlap one another, then this gene will be over-represented by having multiple TSS windows that are the same genomic region


