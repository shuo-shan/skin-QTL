#!/bin/bash
module load bedtools
this_snp=$1

dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/TF_motif_overlapping/QTL_1E-02
rowidx=$(grep -w -n ${this_snp} ${dir}/QTL.bed | cut -d':' -f1)
mkdir temp${rowidx}
cd temp${rowidx}

### fetch SNP bed and fa files
awk -v n=${rowidx} 'NR==n' ${dir}/QTL_100bp.bed > temp${rowidx}.bed
this_snp=$(cat temp${rowidx}.bed | cut -f4)
genome_fa=/share/GHPCC/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa
bedtools getfasta -fi ${genome_fa} -bed temp${rowidx}.bed -fo temp${rowidx}.fa

### run FIMO on SNP against all motifs
export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv12_core_HUMAN_mono_meme_format.meme
#fimo --oc results${rowidx} --thresh 1 --no-qvalue --verbosity 1 ${memeF} temp${rowidx}.fa
fimo --oc results${rowidx} --thresh 0.001 --verbosity 1 ${memeF} temp${rowidx}.fa
cat results${rowidx}/fimo.tsv | grep .H12CORE. | head -100 > results${rowidx}/top100.fimo.tsv

### acquire lock before writing
output_file="${dir}/fimo_output.txt"
lock_file="${dir}/fimo_output.lock"
exec 9>"$lock_file" # opens the lock file for writing and associates it with file descriptor 9. lock file is created if it doesn't already exist.
flock 9 # apply an advisory lock on this file
# write to the output file
cat results${rowidx}/top100.fimo.tsv | awk -v snp=${this_snp} '{OFS=FS="\t"}{print snp,$0}' >> ${output_file}
# release the lock
flock -u 9

### clean-up
cd ${dir}
rm -r temp${rowidx}
