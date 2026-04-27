#!/bin/bash
module load bedtools
module load ucsc_utilities/20240312

workingdir=$1
SNPbed=$2
prefix=$3

########## fetch QTLs that overlap ChIPseq peaks
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
phyloP=/pi/manuel.garber-umw/human/skin/eQTLs/literature/phyloP/hg38.phyloP20way.bw

cd ${workingdir}
cut -f1-4 ${SNPbed} > ${workingdir}/temp_${prefix}.bed
bigWigAverageOverBed ${phyloP} ${workingdir}/temp_${prefix}.bed ${workingdir}/phyloP_scores_${prefix}.txt
rm temp_${prefix}.bed

echo "Wrote output to ${workingdir}/phyloP_scores_${prefix}.txt"; date

