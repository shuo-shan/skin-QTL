#!/bin/bash
module load bedtools

workingdir=$1
snp_id=$2
snp_bed="${3:-}"   # empty string if not provided

########## create SNP bed if needed
if [[ ! -f "${snp_bed}" ]]; then
        echo "Warning, SNP bed not provided, generating now"
        bash /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/01_compile_snp_bed.sh ${workingdir} ${snp_id}
        snp_bed=${workingdir}/${snp_id}.bed
fi

########## fetch QTLs that overlap ChIPseq peaks
genome=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cd /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
outfile=${workingdir}/SNP_overlapping_TF_peaks_${snp_id}.bed
rm -f ${outfile}
for f in *.bed.gz;do
        TF=$(echo ${f} | cut -d'_' -f2 | sed 's/-human//g')
        ID=$(echo ${f} | cut -d'_' -f1)
	bedtools intersect -a ${snp_bed} -b ${f} | cut -f1-6 | awk -v TF=${TF} -v ID=${ID} '{OFS=FS="\t"}{print $0,TF,ID}' >> ${outfile}
        echo ${f}
done
cd ${workingdir}


########## compile summary table
n_lines=$(wc -l < "${outfile}")

echo "Number of overlaps: ${n_lines}"

if [[ "${n_lines}" -eq 0 ]]; then
    echo "No overlap found for ${snp_id}, skipping downstream analysis."
else
    echo "Overlap found, running R script..."
    
    singularity exec /share/pkg/containers/rstudio_example/r_4.5.2.sif \
    Rscript /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/annotate_QTL/03_overlap_tfbs_compile.R \
    ${workingdir} "${outfile}" ${snp_id}
fi
