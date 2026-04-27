#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -e "./step5.5.%J%I.err"
#BSUB -o "./step5.5.%J%I.out"
# script from https://github.com/GuttmanLab/sprite-pipeline/blob/master/Snakefile
# plot heatmap of the clusters

module load condas/2018-05-11
source activate fastQTL
# have R and gplots,  running

prefix=SPRITE-F37-KRT-PBS
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
r_dir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/r
cd ${dir}

Rscript ${r_dir}/plot_heatmap.R -i ${dir}/${prefix}.final.txt -m 255
















