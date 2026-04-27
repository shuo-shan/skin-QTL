#!/bin/bash
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# TUTORIAL: http://homer.ucsd.edu/homer/ngs/peakMotifs.html
dir=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin
outDir=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/HOMER

# INITIAL SET-UP:  generate configured genome file for hg38
#cd ${dir}/HOMER
#module load condas/2018-05-11
#source activate fastQTL #--> change this to whatever conda environment you have.
#conda install -c bioconda homer # ver4.11
#environment location: /home/ss65w/.conda/envs/fastQTL #--> change this to whatever conda environment you have.
#perl /home/ss65w/.conda/envs/fastQTL/share/homer/configureHomer.pl -install hg38 #ver6.3 #--> find your configureHomer.pl location and run this to configure the genome

# Method 1: running Motif enrichment and discovery for a given set of peaks
cd ${dir}/HOMER
module load condas/2018-05-11
source activate fastQTL
PATH=$PATH:/home/ss65w/.conda/envs/fastQTL/share/homer//bin/ #--> change this to conda install notes in the end.
pos_file=${dir}/ATAC_H3K27acChIP_intersect_FRB.bed
findMotifsGenome.pl ${pos_file} hg38 ${outDir} -size given

# Method 2: run Motif enrichment and discover for flanking regions around ATACseq summits that lie inside a putative regulatory element
pos_file=____ # bed file of ATACseq peak summits!!!
findMotifsGenome.pl ${pos_file} hg38 ${outDir} -size 200
