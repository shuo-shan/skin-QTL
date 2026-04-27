#!/bin/bash

#BSUB -n 1
#BSUB -R rusage[mem=50000]
#BSUB -W 124:00
#BSUB -q long # which queue we want to run in
#BSUB -R span[hosts=1]


### 09/07/2021
### script for ChIPseq data processing
### working in /nl/umw_manuel_garber/human/skin/eQTLs/ChIP-seq
### set-up
bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash
module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/ChIP-seq

##############################################################
### calculate distance from each merged peak to the nearest TSS of the highest expressed isoform for each gene
module load bedtools/2.29.2
cd ${Dir}
ct=MEL # cell type
f=${Dir}/peaks/${ct}_merged.bed # merged bed file
# get highest expressed isoform for each gene

# get TSS of highest expressed isoform for each gene
bedtools closest -D ref -t all -a -b ${f} 

