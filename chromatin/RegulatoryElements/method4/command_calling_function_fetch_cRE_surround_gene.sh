#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20400]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./log/%J%I.fetch_cRE_surround_gene.err"
#BSUB -o "./log/%J%I.fetch_cRE_surround_gene.out"

dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method4
tss=/pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_transcription_start_sites.bed
cd ${dir}
cat ${tss} | cut -f5 | sort | uniq > genes.txt
while read gene; do 
        bash ${dir}/function_fetch_cRE_surrounding_gene.sh ${gene} KRT
        bash ${dir}/function_fetch_cRE_surrounding_gene.sh ${gene} FRB
        bash ${dir}/function_fetch_cRE_surrounding_gene.sh ${gene} MEL
        echo ${gene} >> completed_genes.txt
done < genes.txt

