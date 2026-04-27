#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=408000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/%J%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial/%J%I.out"
# preprocess ChromBPNet tutorial data

dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_tutorial
cd ${dir}

### Mapping isogenic reps to yield considated reads.
module load samtools
echo "merging now";date
samtools merge -f ${dir}/data/downloads/merged_unsorted.bam ${dir}/data/downloads/rep1.bam ${dir}/data/downloads/rep2.bam ${dir}/data/downloads/rep3.bam
echo "done merging, sorting now";date
samtools sort -@4 ${dir}/data/downloads/merged_unsorted.bam -o ${dir}/data/downloads/merged.bam
echo "sorting done, indexing now"; date
samtools index ${dir}/data/downloads/merged.bam
echo "indexing done, downloading and processing peak data";date
# for my own data, remove duplicates

### download peak data (they recommend using relaxed peak-calls)
# MACS2 with p-val = 0.01
# download overlap peaks (default peaks on ENCODE)
wget https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz -O ${dir}/data/downloads/overlap.bed.gz
# ensure that the peaks don't overlap blacklist
module load bedtools
bedtools slop -i ${dir}/data/downloads/blacklist.bed.gz -g ${dir}/data/downloads/hg38.chrom.sizes -b 1057 > ${dir}/data/downloads/temp.bed
bedtools intersect -v -a ${dir}/data/downloads/overlap.bed.gz -b ${dir}/data/downloads/temp.bed > ${dir}/data/downloads/peaks_no_blacklist.bed
echo "done with preprocessing step";date

### start singularity image to run chrombpnet interactively: Option 2
# start the gpu session
bsub -Is -q gpu -W 8:00 -R rusage[mem=100G] /bin/bash
# start the singularity image
singularity shell --nv /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/chrombpnet/chrombpnet_latest.sif

### define train, validation and test chr splits
head -n 24 ${dir}/data/downloads/hg38.chrom.sizes > ${dir}/data/downloads/hg38.chrom.subset.sizes
mkdir ${dir}/data/splits
chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ${dir}/data/splits/fold_0
chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr2 chr8 chr9 chr16 -vcr chr12 chr17 -op ${dir}/data/splits/fold_1
chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr4 chr11 chr12 chr15 chrY -vcr chr22 chr7 -op ${dir}/data/splits/fold_2
chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr5 chr10 chr14 chr18 chr20 chr22 -vcr chr6 chr21 -op ${dir}/data/splits/fold_3
chrombpnet prep splits -c ${dir}/data/downloads/hg38.chrom.subset.sizes -tcr chr7 chr13 chr17 chr19 chr21 chrX -vcr chr10 chr18 -op ${dir}/data/splits/fold_4 

### generate non-peaks (background regions)
chrombpnet prep nonpeaks -g ${dir}/data/downloads/hg38.fa -p ${dir}/data/downloads/peaks_no_blacklist.bed -c ${dir}/data/downloads/hg38.chrom.sizes -fl ${dir}/data/splits/fold_0.json -br ${dir}/data/downloads/blacklist.bed.gz -o ${dir}/data/output
