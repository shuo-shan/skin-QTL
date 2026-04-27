#!/bin/bash
# this script in ENCODE ATACseq pipeline performs bowtie2 alignment and post-alignment filtering. It also generates mito fraction log.

# -------------------------------------------------------------------------
# set-up
prefix=$1
dir=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline
trimmed_dir=${dir}/fastqs_trimmed
script_dir=${dir}
alignment_dir=${dir}/alignment
cd ${dir}
mkdir -p ${alignment_dir}

# --------------------------------------------------------------------------
# Align to Human Genome
module load samtools
module load bowtie2

# collect and sort R1 and R2 files
trimmed_r1=${trimmed_dir}/${prefix}_trim.R1.fastq.gz
trimmed_r2=${trimmed_dir}/${prefix}_trim.R2.fastq.gz
bam=${alignment_dir}/${prefix}.bam
log_file=${alignment_dir}/${prefix}.align.log

bowtie2_idx=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/Bowtie2Index/genome
multimapping=4
n_threads=10
echo "bowtie2 alignment started for ${prefix}"; date
# Allow reads with up to $multimapping hits
bowtie2 -k 5 -X1000 --mm --threads ${n_threads} -x ${bowtie2_idx} -1 ${trimmed_r1} -2 ${trimmed_r2} | samtools view -Su - | samtools sort -@ ${n_threads} -o ${bam}
samtools index ${bam}
echo "bowtie2 alignment done"; date

# --------------------------------------------------------------------------
# generate alignment stats
conda activate ATACseq
flagstat_qc=${alignment_dir}/${prefix}.flagstat.qc
echo "generate alignment stats"; date
samtools sort -n --threads 10 ${bam} -O SAM | SAMstats --sorted_sam_file - --outf ${flagstat_qc}

# create a mito-free BAM to calculate fraction of mito reads
echo "create mito-free bam file"
non_mito_bam=${alignment_dir}/${prefix}_non_mito.bam
samtools idxstats ${bam} | cut -f 1 | grep -v -P "^chrM$" | xargs samtools view -@ 2 -b ${bam} > ${non_mito_bam}

# calculate fraction of mitochondria reads
echo "calc fraction of mito reads"
non_mito_samstat_qc=${alignment_dir}/${prefix}_non_mito.samstat.qc
mito_samstat_qc=${alignment_dir}/${prefix}.mito.samstat.qc
mito_frac_logfile=${alignment_dir}/${prefix}.mito_frac_summary.txt

samtools sort -n --threads 10 ${non_mito_bam} -O SAM | SAMstats --sorted_sam_file - --outf ${non_mito_samstat_qc}
R_non_mito=$(grep "mapped (" $non_mito_samstat_qc | awk '{print $1}')

samtools sort -n --threads 10 "$bam" -O SAM | SAMstats --sorted_sam_file - --outf "$mito_samstat_qc"
R_total=$(grep "mapped (" $mito_samstat_qc | awk '{print $1}')
R_mito=$(( R_total - R_non_mito ))
fraction_mito_reads=$(echo "scale=4; $R_mito / $R_total" | bc)

echo "R_total: $R_total" > $mito_frac_logfile
echo "R_non_mito: $R_non_mito" > $mito_frac_logfile
echo "R_mito: $R_mito" >> $mito_frac_logfile
echo "Fraction of mito reads R_mito / R_total: $fraction_mito_reads" >> $mito_frac_logfile

# --------------------------------------------------------------------------
# Post alignment filtering
#  - remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 524), reads mapping to chrM.
#  - retain properly paired reads -f 2
#  - If max. number of multimappers is activated then use explicit multimapping filtering code
#  - Remove PCR duplicates (using Picard’s MarkDuplicates or FixSeq)

filtered_alignment_dir=${dir}/filtered_alignment
mkdir -p "$filtered_alignment_dir"

raw_bam_file=${alignment_dir}/${prefix}.bam
raw_alignment_file=${alignment_dir}/${prefix}_non_mito.bam

echo "Post alignment filtering for: $raw_bam_file"; date
# =============================
# remove  unmapped, mate unmapped, too short
# not primary alignment, reads failing platform
# only keep properly paired reads
# obtain name sorted BAM file
# ==================
filt_bam_file="${filtered_alignment_dir}/${prefix}.filt.bam"
tmp_filt_bam_prefix="tmp.${prefix}.filt.nmsrt"
tmp_filt_bam_file="${filtered_alignment_dir}/${tmp_filt_bam_prefix}.bam"
tmp_filt_fixmate_bam_file="${filtered_alignment_dir}/${tmp_filt_bam_prefix}.fixmate.bam"
tmp_filt_multimapper_bam_file="${filtered_alignment_dir}/${tmp_filt_bam_prefix}.multimapper.bam"

# note: flag 524 = 512 (QC fail) + 8 (mate unmapped) + 4 (read unmapped), -f 2: include only properly paired reads, -u: output in uncompressed BAM for piping
samtools view -F 524 -f 2 -u ${raw_bam_file} | \
awk 'BEGIN{OFS="\t"}
     /^@/ {print; next}
     {
        cigar = $6;
        len = 0;
        while (match(cigar, /^[0-9]+[MDN=X]/)) {
            len += substr(cigar, RSTART, RLENGTH - 1)
            cigar = substr(cigar, RSTART + RLENGTH)
        }
        if (len > 4) print
     }' | samtools sort -n -@ 10 -o ${tmp_filt_bam_file}  
# note: custom script to 
${script_dir}/assign_multimapper.py -k ${multimapping} --paired-end ${tmp_filt_bam_file} ${tmp_filt_multimapper_bam_file}
samtools view -h $tmp_filt_multimapper_bam_file | samtools fixmate -r - ${tmp_filt_fixmate_bam_file}

# note: Remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
samtools view -F 1804 -f 2 -u ${tmp_filt_fixmate_bam_file} | samtools sort -@ 10 -o ${filt_bam_file} -


# --------------------------------------------------
# clean-up
rm ${tmp_filt_multimapper_bam_file}
rm ${tmp_filt_fixmate_bam_file}
rm ${tmp_filt_bam_file}




    





































