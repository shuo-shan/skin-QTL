#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -e "./step6.%J%I.err"
#BSUB -o "./step6.%J%I.out"
# script from https://github.com/GuttmanLab/sprite-pipeline/blob/master/scripts/python/get_hub_contacts.py
# get hub averages and print

module load condas/2018-05-11
source activate sprite

prefix=SPRITE-F37-KRT-PBS
python_dir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/python
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cd ${dir}


python ${python_dir}/get_hub_contacts.py \
	--heatmap ${dir}/${prefix}.raw.txt \
	--resolution 1000000 \
	--hub ${dir}/hub.bed \
	--assembly hg38 

