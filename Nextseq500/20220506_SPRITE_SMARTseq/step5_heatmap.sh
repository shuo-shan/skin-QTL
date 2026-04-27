#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=400000]
#BSUB -q long
#BSUB -W 24:00
#BSUB -e "./step5.5kb.%J%I.err"
#BSUB -o "./step5.5kb.%J%I.out"
# script from https://raw.githubusercontent.com/wiki/GuttmanLab/sprite-pipeline/5.-Heatmaps.md
# group all reads with the same barcodes into the same cluster. Also discard PCR duplicates.

module load condas/2018-05-11
source activate sprite
# python 3.9.12
# pysam 0.19.0

prefix=SPRITE-F37-KRT-PBS
python_dir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/python
HiCorrector_dir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/HiCorrector_1.2/bin/ic
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cd ${dir}

## same set-up from tutorial
#python ${python_dir}/get_sprite_contacts.py \
#	--clusters ${dir}/${prefix}.clusters \
#	--raw_contacts ${dir}/${prefix}.raw.txt \
#	--biases ${dir}/${prefix}.bias.txt \
#	--iced ${dir}/${prefix}.iced.txt \
#	--output ${dir}/${prefix}.final.txt \
#	--assembly hg38 \
#	--chromosome genome \
#	--min_cluster_size 2 \
#	--max_cluster_size 1000 \
#	--resolution 1000000 \
#	--downweighting none \
#	--hicorrector ${HiCorrector_dir} \
#	--iterations 100
	
# now let's try 5kb resolution.
python ${python_dir}/get_sprite_contacts.py \
       --clusters ${dir}/${prefix}.clusters \
       --raw_contacts ${dir}/${prefix}.raw.txt \
       --biases ${dir}/${prefix}.bias.txt \
       --iced ${dir}/${prefix}.iced.txt \
       --output ${dir}/${prefix}.final.txt \
       --assembly hg38 \
       --chromosome genome \
       --min_cluster_size 2 \
       --max_cluster_size 1000 \
       --resolution 10000 \
       --downweighting none \
       --hicorrector ${HiCorrector_dir} \
       --iterations 100
