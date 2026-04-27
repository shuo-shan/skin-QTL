#!/bin/bash
#BSUB -n 11
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=10000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -e "./step2.%J%I.err"
#BSUB -o "./step2.%J%I.out"
# script from: https://raw.githubusercontent.com/wiki/GuttmanLab/sprite-pipeline/2.-Alignment.md

module load bowtie2/2.4.1
module load samtools/1.9

prefix=SPRITE-F37-KRT-PBS
bowtie2_index_dir=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/Bowtie2Index/genome
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cd ${dir}

bowtie2 -p 10 -t --phred33 -x ${bowtie2_index_dir} -U ${dir}/${prefix}.read1.barcoded.fastq.gz |\
samtools view -bq 20 -F 4 -F 256 - > ${dir}/${prefix}.DNA.bowtie2.mapq20.bam

