module load bcftools
module load bedtools
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/TF_peak_overlapping/QTL_1E-02
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


########## fetch QTLs that overlap ChIPseq peaks
QTL=${dir}/QTL.bed
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
rm QTL_overlapping_TF_peaks.bed
for f in *.bed.gz;do
        TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
        ID=$(echo ${f} | cut -d'_' -f1)
	bedtools intersect -a ${QTL} -b ${f} | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> QTL_overlapping_TF_peaks.bed
        echo ${f}
done
mv QTL_overlapping_TF_peaks.bed ${dir}
cd ${dir}
Rscript ${dir}/intersect_reQTL_with_TFChIPseqPeaks.R # -----> this create the sorted bed file of high-level summary
