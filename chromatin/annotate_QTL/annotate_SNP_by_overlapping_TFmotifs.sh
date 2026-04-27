#!/bin/bash
module load bcftools
module load bedtools

#
workingdir=$1
SNPbed=$2
prefix=$3

## toy example
#workingdir=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/transQTL/temp_output/rs10416689
#SNPbed=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/2025models/FRB/transQTL/temp_output/rs10416689/QTL_temp.bed
#prefix=temp

##############
cd ${workingdir}
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
bedtools slop -b 50 -i ${SNPbed} -g ${genome} > ${workingdir}/QTL_100bp_${prefix}.bed

#############
mkdir log
jobname=fimo${prefix}
cut -f4 ${workingdir}/QTL.bed > ${workingdir}/QTL.txt
while read snp;do
	echo "bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/run_fimo_perSNP.sh ${workingdir} ${snp} ${prefix}" >> ${workingdir}/commands_runFimo_${prefix}.txt
done < ${workingdir}/QTL.txt
while read c;do
	echo ${c} | bsub -W 03:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=800] -q long -e "${workingdir}/log/fimo_%J%I.err" -o "${workingdir}/log/fimo_%J%I.out"
done < ${workingdir}/commands.joined.txt

while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; date;sleep 120; done
# this creates a fimo_output.txt file

############# filter output file to motifs that contain the SNP
cat ${dir}/fimo_output.txt | awk '{OFS="\t"}{print $3,$4,$5,$1,$6,$2,$7,$8,$9,$10}' > fimo_output.bed
bedtools intersect -wb -a QTL.bed -b fimo_output.bed > fimo_output_SNP-containing.bed
# 69,009,013 entries --> 14,732,289 entries

############# filter by a lenient q-value 1E-04
cat fimo_output_SNP-containing.bed | awk '$15<0.0001' > fimo_output_SNP-containing_qval1E-4.bed # 686,501 entries

############# filter by a moderate q-value 1E-06
cat fimo_output_SNP-containing.bed | awk '$15<0.000001' > fimo_output_SNP-containing_qval1E-6.bed # 142,678 entries

############# filter by a strict q-value 1E-08. result seems to capture bad motifs like AAAAAAAA or ATATATAT
cat fimo_output_SNP-containing.bed | awk '$15<0.00000001' > fimo_output_SNP-containing_qval1E-8.bed # 38,557 entries

############# filter by a lenient q-value 1E-04 but also the TF from ENCODE TF ChIPseq dataset
cat  /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/TF_peak_overlapping/QTL_1E-02/QTL_overlapping_TF_peaks.bed | cut -f7 | sort | uniq > TF.txt
cat fimo_output_SNP-containing.bed | grep -w -f TF.txt | awk '$15<0.0001' > fimo_output_SNP-containing_qval1E-4_ENCODETF.bed # 15,754 entries

############# compile high-level summary in R
Rscript compile_summary.R
# final output is: QTL_1E-02_overlapping_TFmotif.bed













