#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -e "./step4.%J%I.err"
#BSUB -o "./step4.%J%I.out"
# script from https://raw.githubusercontent.com/wiki/GuttmanLab/sprite-pipeline/4.-Clustering.md
# group all reads with the same barcodes into the same cluster. Also discard PCR duplicates.

module load condas/2018-05-11
source activate sprite_env1
# python 3.9.12
# pysam 0.19.0

prefix=SPRITE-F37-KRT-PBS
python_dir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/python
#python_dir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite2.0-pipeline-master/scripts/python
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cd ${dir}

python ${python_dir}/get_clusters.py --input ${dir}/${prefix}.DNA.bowtie2.mapq20.masked.sorted.bam --output ${dir}/${prefix}.clusters --num_tags 5


