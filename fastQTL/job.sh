#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=200000]
#BSUB -q long
#BSUB -W 121:00
### script to interface Gencove to retrieve, upload, and process genotyping data

### crystal shan 09/2021
### bsub -Is -q interactive -W 8:00 -n1 -R rusage[mem=450000] -R "span[hosts=1]" /bin/bash

############################
#### packages & env
date
module load condas/2018-05-11
source activate sshan_isoform
Rscript genotype_PCA.R
date
echo "done with PCA."
