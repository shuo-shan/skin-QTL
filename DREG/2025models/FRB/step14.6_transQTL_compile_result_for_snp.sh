#!/bin/bash

snp=$1
ct=$2
condition=$3
QTLtype=$4

#snp=rs1010167
#ct=FRB
#condition=TNF
#QTLtype=eQTL

dir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/${ct}
workingdir=${dir}/transQTL/temp_output/${snp}

# overlap with ATACseq peaks
bedF=${workingdir}/SNP_genotype_${ct}_${condition}_${QTLtype}.bed
if [[ $(wc -l < "${bedF}") -gt 1 ]]; then
	bash ${dir}/step14.6_overlap_SNPgenotypeBED_with_ATACseq_peak.sh ${workingdir}/SNP_genotype_${ct}_${condition}_${QTLtype}.bed
fi

# if cisGene has GWAS colocalization
colocF=${dir}/transQTL/coloc_table_${ct}_${condition}_${QTLtype}.txt
if grep -qw "${snp}" "${colocF}"; then
    grep -w "${snp}" "${colocF}" > "${workingdir}/coloc_${ct}_${condition}_${QTLtype}.txt"

    trait=$(grep -w ${snp} ${colocF} | cut -f3)
    g=$(grep -w ${snp} ${colocF} | cut -f2)
    
    cp ${dir}/coloc/${trait}/plots/${ct}_${g}.locus_tracks.pdf ${workingdir}/
    cp ${dir}/coloc/${trait}/plots_table/${ct}_${g}.table.pdf ${workingdir}/
fi

#while read line;do
#	snp=$(echo ${line} | cut -d' ' -f1)
#	g=$(echo ${line} | cut -d' ' -f2)
#	trait=$(echo ${line} | cut -d' ' -f3)
#
#	workingdir=${dir}/transQTL/temp_output/${snp}
#	cp ${dir}/coloc/${trait}/plots/${ct}_${g}.locus_tracks.pdf ${workingdir}/
#	cp ${dir}/coloc/${trait}/plots_table/${ct}_${g}.table.pdf ${workingdir}/
#done < ${dir}/transQTL/coloc_table_${ct}_${condition}_${QTLtype}.txt


# cis-QTL SNP annotations
bigtable=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/QTL_1E-02/compiled_table_${ct}.txt
head -1 ${bigtable} | tr '\t' '\n' | awk '{print NR,$0}' > ${workingdir}/rownames.txt
grep -w ${snp} ${bigtable} | tr '\t' '\n' > ${workingdir}/body.txt
if [ $(wc -l < ${workingdir}/body.txt) -gt 1 ]; then
        echo "yaay!"
        paste ${workingdir}/rownames.txt ${workingdir}/body.txt > ${workingdir}/annotation_${snp}.txt
        rm ${workingdir}/body.txt
fi
rm ${workingdir}/rownames.txt

# make plots
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="cisGene") col=i} NR>1{print $col}' \
    ${workingdir}/compiled_table_${ct}_${condition}_${QTLtype}.txt | sort -u > ${workingdir}/genes.txt
awk 'NR==1{for(i=1;i<=NF;i++) if($i=="transGene") col=i} NR>1{print $col}' \
    ${workingdir}/compiled_table_${ct}_${condition}_${QTLtype}.txt | sort -u >> ${workingdir}/genes.txt

while read gene;do
	bash ${dir}/step6_plotGeneSnpPair_any_pair.sh ${gene} ${snp}
	echo "${gene} done"
	mv ${dir}/plots/temp_output/plot_*_${gene}_${snp}.pdf ${workingdir}
done < ${workingdir}/genes.txt
rm ${workingdir}/plot_PBS_*.pdf


