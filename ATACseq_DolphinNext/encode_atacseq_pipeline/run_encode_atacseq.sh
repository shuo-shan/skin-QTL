#!/bin/bash

eval "$(conda shell.bash hook)"  # Ensures Conda is properly initialized
conda activate bifx2 || { echo "Error: Conda environment 'bifx2' not found"; exit 1; }

picard_path=$(readlink -f $(which picard))
export PICARD=$(dirname $picard_path)/picard.jar

assign_multimapper="/home/azita.ghodssi-umw/eQTL_chrombpnet/assign_multimapper.py"


# Check if cutadapt is installed
if ! command -v cutadapt &> /dev/null; then
    echo "Error: cutadapt not found in environment 'bifx2'"
    exit 1
fi

#-------------------------------------------------
# Set up the Env
DEBUG=1
INFO=1
script_dir="/home/azita.ghodssi-umw/eQTL_chrombpnet/atac-seq-pipeline/src"
#fastq_dir="/pi/manuel.garber-umw/human/skin/eQTLs/ATAC-seq/fastqs/merged"
fastq_dir="./fastqs"
out_dir="./out"

if [ ! -d $out_dir ]; then
    mkdir -p "$trimmed_dir"
fi


# references
source ./refs/hg38_refs.tsv


#-------------------------------------------------
# Input FastQ files

r1_files=()
r2_files=()

for file in "$fastq_dir"/*.fastq.gz; do
    if [[ $file =~ .R1 ]]; then
	r1_files+=("$file")
    elif [[ $file =~ .R2 ]]; then
	r2_files+=("$file")
    fi
done

IFS=$'\n' r1_files=($(sort <<<"${r1_files[*]}"))
IFS=$'\n' r2_files=($(sort <<<"${r2_files[*]}"))

# ensure both arrays have the same length
if [ ${#r1_files[@]} -ne ${#r2_files[@]} ]; then
    echo "Error: Mismatched number of R1 and R2 files."
    exit 1
fi

if [[ "$DEBUG" -eq 1 ]]; then
    r1_files=("${r1_files[0]}")
    r2_files=("${r2_files[0]}")
fi


if [[ "$INFO" -eq 1 ]]; then
    echo "R1 files:"
    echo "========="
    echo "${r1_files[*]}" | sed 's/ /\n/g'

    echo "R2 files:"
    echo "========="
    echo "${r2_files[*]}" | sed 's/ /\n/g'
fi

#-------------------------------------------------
# Run FASTQC

# ...


#-------------------------------------------------
# Detect and remove FASTQ adapters
adapter_err_rate=0.2
trimmed_dir="$out_dir/trimmed_fastqs"
n_threads=10

if [ ! -d $trimmed_dir ]; then
    mkdir -p "$trimmed_dir"
fi


for ((i=0; i<${#r1_files[@]}; i++)); do
    r1="${r1_files[i]}"
    r2="${r2_files[i]}"

    # Extract prefix from the filename (remove path and extension)
    prefix=$(basename "$r1" .fastq.gz)
    prefix="${prefix//.R1/_R1}"
    
    trimmed_r1="$trimmed_dir/${prefix}_trim_R1.fastq.gz"
    trimmed_r2="$trimmed_dir/${prefix/_R1/_R2}_trim_R2.fastq.gz"

    echo "Processing: $r1 and $r2"

    if [ -e "$trimmed_r1" ] && [ -e "$trimmed_r2" ]; then
        echo "$trimmed_r1 and $trimmed_r2 exist"
    else
        # Log files for adapter detection
        log_r1="$trimmed_dir/${prefix}_adapter_R1.txt"
        log_r2="$trimmed_dir/${prefix/_R1/_R2}_adapter_R2.txt"

        # Detect adapter sequences separately for R1 and R2
        python3 "$script_dir/detect_adapter.py" "$r1" > "$log_r1"
        python3 "$script_dir/detect_adapter.py" "$r2" > "$log_r2"

        # Parse logs to extract adapter sequences
        adapter_r1=$(awk 'NR==1 {print $1}' "$log_r1")
        adapter_r2=$(awk 'NR==1 {print $1}' "$log_r2")

        echo "Detected adapter for R1: $adapter_r1"
        echo "Detected adapter for R2: $adapter_r2"

        # FASTQ read adapter trimming with separate adapters
        log_file="$trimmed_dir/${prefix}_cutadapt.log"

	#bsub -o $log_file  -q long -W 10:00 -n 20 -R span[hosts=1] -R rusage[mem=1000] "cutadapt -j $n_threads -m 5 -e $adapter_err_rate  -a $adapter_r1 -A $adapter_r2 --pair-filter=both -o $trimmed_r1 -p $trimmed_r2 $r1 $r2"

        cutadapt -j "$n_threads" -m 10 -e "$adapter_err_rate"  -a "$adapter_r1" -A "$adapter_r2" --pair-filter=both  -o "$trimmed_r1" -p "$trimmed_r2" "$r1" "$r2" 2> "$log_file"

    fi
done


# --------------------------------------------------------------------------
# Align to Human Genome

alignment_dir="$out_dir/alignment"
multimapping=4

if [ ! -d $alignment_dir ]; then
    mkdir -p "$alignment_dir"
fi

# collect and sort R1 and R2 files
trimmed_r1_files=($(ls "${trimmed_dir}"/*R1.fastq.gz | sort))
trimmed_r2_files=($(ls "${trimmed_dir}"/*R2.fastq.gz | sort))

# ensure both arrays have the same length
if [ ${#trimmed_r1_files[@]} -eq 0 ]; then
    echo "Error: No trimmed fastq files found"
    exit 1
fi

if [ ${#trimmed_r1_files[@]} -ne ${#trimmed_r2_files[@]} ]; then
    echo "Error: Mismatched number of R1 and R2 files."
    exit 1
fi

# Loop through paired files
for ((i=0; i<${#trimmed_r1_files[@]}; i++)); do
    fq1="${trimmed_r1_files[i]}"
    fq2="${trimmed_r2_files[i]}"
    echo "Processing: $fq1, $fq2"

    prefix=$(basename "$fq1" .fastq.gz)
    prefix=$(echo "$prefix" | sed 's/_R1//g; s/_trim//g')
    
    bam="$alignment_dir/$prefix.bam"
    log_file="$alignment_dir/$prefix.align.log"

    n_threads=10
    if [ -e $bam ]; then
	echo "$bam exits"
    else
	if [ "$multimapping" -eq 0 ]; then
            # unique mapping reads only
	    echo "bowtie2 unique mapper"
	    #bsub -o $log_file -q long -W 10:00 -n 20 -R span[hosts=1] -R rusage[mem=1000] "bowtie2 -X2000 --mm --threads $n_threads -x $bowtie2_idx -1 $fq1 -2 $fq2  | samtools view -Su - | samtools sort -o $bam"
            bowtie2 -X2000 --mm --threads $n_threads -x $bowtie2_idx -1 $fq1 -2 $fq2 2>$log | samtools view -Su - | samtools sort -o "$bam"

	else
	    echo "bowtie2 multi-mapper"
            # Allow reads with up to $multimapping hits
	    #bsub -o $log_file -q long -W 10:00 -n 20 -R span[hosts=1] -R rusage[mem=1000] "bowtie2 -k $((multimapping+1)) -X2000 --mm --threads $n_threads -x $bowtie2_idx -1 $fq1 -2 $fq2 | samtools view -Su - | samtools sort -o $bam"
	    bowtie2 -k $((multimapping+1)) -X2000 --mm --threads $n_threads -x $bowtie2_idx -1 $fq1 -2 $fq2 | samtools view -Su - | samtools sort -o $bam
	    samtools index $bam
	fi
    fi

    
    # generate alignment stats
    flagstat_qc="$alignment_dir/$prefix.flagstat.qc"
    echo "generate alignment stats"
    if [ -e $flagstat_qc ]; then
	echo "$flagstat_qc exists"
    else
	samtools sort -n --threads 10 "$bam" -O SAM | SAMstats --sorted_sam_file - --outf "$flagstat_qc"
    fi

    # create a mito-free BAM to calculate fraction of mito reads
    echo "create mito-free bam file"
    non_mito_bam="$alignment_dir/${prefix}_non_mito.bam"
    if [ -e $non_mito_bam ]; then
	echo "$non_mito_bam exists"
    else
	samtools idxstats "$bam" | cut -f 1 | grep -v -P "^chrM$" | xargs samtools view -@ 2 -b "$bam" > "$non_mito_bam"
    fi
    
    # calculate fraction of mitochondria reads
    echo "calc fraction of mito reads"
    non_mito_samstat_qc="$alignment_dir/$prefix_non_mito.samstat.qc"
    mito_samstat_qc="$alignment_dir/$prefix.mito.samstat.qc"
    mito_frac_logfile="$alignment_dir/$prefix.mito_frac_summary.txt"

    if [ -e $mito_frac_logfile ]; then
	echo "$mito_frac_logfile exists"
    else
	samtools sort -n --threads 10 "$non_mito_bam" -O SAM | SAMstats --sorted_sam_file - --outf "$non_mito_samstat_qc"
	R_non_mito=$(grep "mapped (" $non_mito_samstat_qc | awk '{print $1}')

	samtools sort -n --threads 10 "$bam" -O SAM | SAMstats --sorted_sam_file - --outf "$mito_samstat_qc"
	R_total=$(grep "mapped (" $mito_samstat_qc | awk '{print $1}')
	R_mito=$(( R_total - R_non_mito ))

	fraction_mito_reads=$(echo "scale=4; $R_mito / $R_total" | bc)

	echo "R_total: $R_total" > $mito_frac_logfile
	echo "R_non_mito: $R_non_mito" > $mito_frac_logfile
	echo "R_mito: $R_mito" >> $mito_frac_logfile
	echo "Fraction of mito reads R_mito / R_total: $fraction_mito_reads" >> $mito_frac_logfile
     fi
done

# --------------------------------------------------------------------------
# Post alignment filtering
#  - remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 524), reads mapping to chrM.
#  - retain properly paired reads -f 2
#  - If max. number of multimappers is activated then use explicit multimapping filtering code
#  - If unique mapping is activated then remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)
#  - emove PCR duplicates (using Picard’s MarkDuplicates or FixSeq)


filtered_alignment_dir="${out_dir}/filtered_alignment"
if [ ! -d $filtered_alignment_dir ]; then
    mkdir -p "$filtered_alignment_dir"
fi

raw_alignment_files=($(find "${alignment_dir}/" -maxdepth 1 -type f -name "*_non_mito.bam" | sort))

# ensure there are alignment files
if [ ${#raw_alignment_files[@]} -eq 0 ]; then
    echo "Error: No alignment file found."
    exit 1
fi


# Loop through alignment files
for ((i=0; i<${#raw_alignment_files[@]}; i++)); do
    raw_bam_file="${raw_alignment_files[i]}"
    echo "Processing: $raw_bam_file"

    prefix=$(basename "$raw_bam_file" .bam)
    prefix=$(echo "$prefix" | sed 's/_non_mito//g;')

    # =============================
    # remove  unmapped, mate unmapped
    # not primary alignment, reads failing platform
    # only keep properly paired reads
    # obtain name sorted BAM file
    # ==================
    qname_sort_bam_file=”${prefix}.qnmsrt.bam”
    filt_bam_prefix="${prefix}.filt"
    filt_bam_file="${filtered_alignment_dir}/${filt_bam_prefix}.bam"

    tmp_filt_bam_prefix="tmp.${filt_bam_prefix}.nmsrt"
    tmp_filt_bam_file="${filtered_alignment_dir}/${tmp_filt_bam_prefix}.bam"
    tmp_filt_fixmate_bam_file="${filtered_alignment_dir}/${tmp_filt_bam_prefix}.fixmate.bam"

    n_threads=10
    if [ -e $filt_bam_file ]; then
	echo "$tmp_filt_fixmate_bam_file exists"
    else
	if [ "$multimapping" -eq 0 ]; then
	    mapq_threshold=30
	    samtools view -F 1804 -f 2 -q ${mapq_threshold} -u ${raw_bam_file} | samtools sort -n -o ${tmp_filt_bam_file} -
	    samtools fixmate -r ${tmp_filt_bam_file} ${tmp_filt_fixmate_bam_file}
	else
	    # Remove  unmapped, mate unmapped
	    # not primary alignment, reads failing platform
	    # Only keep properly paired reads
	    # Obtain name sorted BAM file
	    samtools view -F 524 -f 2 -u ${raw_bam_file} | samtools sort -n -o ${tmp_filt_bam_file}  -
	    
	    tmp_filt_multimapper_bam_file="${filtered_alignment_dir}/${tmp_filt_bam_prefix}.multimapper.bam"
	    $assign_multimapper -k $multimapping --paired-end ${tmp_filt_bam_file} $tmp_filt_multimapper_bam_file
	    samtools view -h $tmp_filt_multimapper_bam_file | samtools fixmate -r - ${tmp_filt_fixmate_bam_file}
	    #rm ${tmp_filt_multimapper_bam_file}
        fi


	# Remove orphan reads (pair was removed)
	# and read pairs mapping to different chromosomes
	# Obtain position sorted BAM
	samtools view -F 1804 -f 2 -u ${tmp_filt_fixmate_bam_file} | samtools sort -o ${filt_bam_file} -
	rm ${tmp_filt_fixmate_bam_file}
	rm ${tmp_filt_bam_file}
    fi


    # 
    # Mark duplicates
    #
    set -x
    tmp_filt_bam_file="${filtered_alignment_dir}/${filt_bam_prefix}.dupmark.bam"
    dup_file_qc="${filtered_alignment_dir}/${filt_bam_prefix}.dup.qc"

    java -Xmx4G -jar $PICARD MarkDuplicates INPUT=${filt_bam_file} OUTPUT=${tmp_filt_bam_file} METRICS_FILE=${dup_file_qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
    mv ${tmp_filt_bam_file} ${filt_bam_file}

    exit

    # 
    # Remove duplicates
    # Index final position sorted BAM
    # Create final name sorted BAM
    #
    FINAL_BAM_PREFIX="${OFPREFIX}.nodup"
    FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
    FINAL_BAM_INDEX_FILE="${FINAL_BAM_FILE}.bai"
    FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

    samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}
    # Index Final BAM file
    samtools index ${FINAL_BAM_FILE}
    samtools sort -n --threads 10 ${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${FINAL_BAM_FILE_MAPSTATS}

done

exit


    #
    # Compute library complexity
    #
    # - sort by name
    # - convert to bedPE and obtain fragment coordinates
    # - sort by position and strand
    # - obtain unique count statistics
    PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

    # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab]
    #       TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

    samtools sort -n ${FILT_BAM_FILE} -o ${OFPREFIX}.srt.tmp.bam
    bedtools bamtobed -bedpe -i ${OFPREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}

    rm ${OFPREFIX}.srt.tmp.bam
    rm ${FILT_BAM_FILE}



    





































