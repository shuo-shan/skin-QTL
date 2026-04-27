#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=80400]
#BSUB -q long
#BSUB -W 8:00
#BSUB -e "./step1.%J%I.err"
#BSUB -o "./step1.%J%I.out"

module load java/1.8.0_77

prefix=SPRITE-F37-KRT-PBS
scriptDir=/nl/umw_manuel_garber/human/skin/eQTLs/SPRITE/pipeline/sprite-pipeline-master/scripts/java
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cd ${dir}

java -jar ${scriptDir}/BarcodeIdentification_v1.2.0.jar \
  --input1 ${dir}/SPRITE-F37-KRT-PBS_R1.fastq.gz \
  --input2 ${dir}/SPRITE-F37-KRT-PBS_R2.fastq.gz \
  --output1 ${dir}/${prefix}.read1.barcoded.fastq.gz \
  --output2 ${dir}/${prefix}.read2.barcoded.fastq.gz \
  --config ${scriptDir}/example_config.txt 

