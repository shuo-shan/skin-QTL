#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=100000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

module load condas/2018-05-11
source activate sshan_isoform

module load bowtie/1.3.0
module load samtools/1.9
module load perl/5.28.1
module load python3/3.5.0
module load python3/3.5.0_packages/numpy/1.18.5


dir=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/HiCorr/HiCorr
bamdir=/nl/umw_manuel_garber/human/skin/eQTLs/literature/fibroblast/dixon_lung_fibroblast_imr90_hic/
name_of_data=IMR90
bamf=${bamdir}/${name_of_data}.bam
read_length=36
genome=/share/data/umw_biocore/genome_data/human/hg19/hg19.fa

# run HiCorr preprocess to turn bam file into loop files
#/HiCorr Bam-process-HindIII <bam_file> <name_of_your_data> <mapped_read_length_in_your_bam_file> <genome> HindIII
${dir}/HiCorr Bam-process-HindIII ${bamf} ${name_of_data} ${read_length} ${genome} HindIII
wait
echo "done turning bam file to loop files"; date

# run HiCorr bias correction using two *frag_loop files
${dir}/HiCorr HindIII ${name_of_data}.cis.frag_loop ${name_of_data}.trans.frag_loop ${name_of_data} hg19 
wait
echo "done running HiCorr";date

