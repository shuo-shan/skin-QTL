module load bcftools
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/TF_motif_overlapping/QTL_1E-02
cd ${dir}

########## QTL pval cutoff is 1E-02 from rankNormCPM featureSelected model
# fetch MEL QTL 
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/MEL_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > MEL_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > MEL_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > MEL_IFNeQTL.txt
cat MEL_*QTL.txt | sort | uniq > MEL_QTL.txt # (334,722)
# fetch KRT QTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/KRT_minimal/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > KRT_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > KRT_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > KRT_IFNeQTL.txt
cat KRT_*QTL.txt | sort | uniq > KRT_QTL.txt # (320,971)
# fetch FRB QTL
f=/pi/manuel.garber-umw/human/skin/eQTLs/DREG/FRB_new/masteroutput_all_with_colnames.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($7!="." && $7<0.01) print $1}' | sort | uniq > FRB_reQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($15!="." && $15<0.01) print $1}' | sort | uniq > FRB_PBSeQTL.txt
cat ${f} | awk '{OFS=FS="\t"}{if ($23!="." && $23<0.01) print $1}' | sort | uniq > FRB_IFNeQTL.txt
cat FRB_*QTL.txt | sort | uniq > FRB_QTL.txt # (323,391)
# get unique QTLs and the vcf file and bed file
cat MEL_QTL.txt KRT_QTL.txt FRB_QTL.txt | sort | uniq > QTL.txt #(691,034)
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/vcf/merged/bd841628-fcc2-487a-8460-f5428237f0c9.merged.with_missing_allele.minAC3.AFtagged.vcf.gz
bcftools view --include ID==@QTL.txt ${vcf} -Oz -o ${dir}/QTL.vcf.gz
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\t%ID\n' ${dir}/QTL.vcf.gz | awk '{OFS=FS="\t"}{print $1,$2-1,$2,$6,$4,$5}' | bedtools sort -i stdin > ${dir}/QTL.bed


##############
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
bedtools slop -b 50 -i ${dir}/QTL.bed -g ${genome} > QTL_100bp.bed

#############
mkdir log
jobname=fimo
cut -f4 ${dir}/QTL.bed > ${dir}/QTL.txt
while read snp;do
	echo "bash ${dir}/run_fimo_perSNP.sh ${snp}" >> ${dir}/commands.txt
done < ${dir}/QTL.txt
bash /pi/manuel.garber-umw/sshan/scripts/function_collapse_commands.sh ${dir} ${dir}/commands.txt ${dir}/commands.joined.txt 1000
while read c;do
	echo ${c} | bsub -W 03:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=800] -q long -e "${dir}/log/fimo_%J%I.err" -o "${dir}/log/fimo_%J%I.out"
done < ${dir}/commands.joined.txt

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













