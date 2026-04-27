#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20000]
#BSUB -q long
#BSUB -W 1:00
#BSUB -e "./step4.5.%J%I.err"
#BSUB -o "./step4.5.%J%I.out"
# script from https://github.com/GuttmanLab/sprite-pipeline/blob/master/Snakefile
# plot bar graph of cluster sizes

module load condas/2018-05-11
source activate fastQTL
# this env has required R packages

prefix=SPRITE-F37-KRT-PBS
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
r_dir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/r
cd ${dir}

# plot cluster size distribution
Rscript ${r_dir}/get_cluster_size_distribution.r ${dir} ${prefix}.clusters

mv cluster_sizes.pdf ${dir}/${prefix}.cluster_sizes.pdf
mv cluster_sizes.png ${dir}/${prefix}.cluster_sizes.png

















