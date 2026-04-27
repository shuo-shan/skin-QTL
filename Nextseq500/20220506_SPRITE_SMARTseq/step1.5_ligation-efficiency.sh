#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=80400]
#BSUB -q long
#BSUB -W 8:00
#BSUB -e "./step1.5.%J%I.err"
#BSUB -o "./step1.5.%J%I.out"

module load condas/2018-05-11
source activate sprite_env1
# python 3.9.12
# pysam 0.19.0

prefix=SPRITE-F37-KRT-PBS
scriptDir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/python
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cd ${dir}

#Get ligation efficiency
python ${scriptDir}/get_ligation_efficiency.py ${dir}/${prefix}.read1.barcoded.fastq.gz > ${dir}/${prefix}.ligation_efficiency.txt
