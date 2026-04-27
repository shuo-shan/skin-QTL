# this script runs adaptor detection and trimming
prefix=$1

dir=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/encode_atacseq_pipeline
trimmed_dir=${dir}/fastqs_trimmed
script_dir=${dir}
cd ${dir}
mkdir -p ${trimmed_dir}

r1=${dir}/fastq/${prefix}.R1.fastq.gz
r2=${dir}/fastq/${prefix}.R2.fastq.gz
trimmed_r1=${trimmed_dir}/${prefix}_trim.R1.fastq.gz
trimmed_r2=${trimmed_dir}/${prefix}_trim.R2.fastq.gz

# Detect adapter sequences separately for R1 and R2
echo "Processing: $r1 and $r2"; date
python3 ${script_dir}/detect_adapter.py ${r1} > ${trimmed_dir}/${prefix}_adapter_R1.txt
python3 ${script_dir}/detect_adapter.py ${r2} > ${trimmed_dir}/${prefix}_adapter_R2.txt
adapter_r1=$(awk 'NR==1 {print $1}' ${trimmed_dir}/${prefix}_adapter_R1.txt)
adapter_r2=$(awk 'NR==1 {print $1}' ${trimmed_dir}/${prefix}_adapter_R2.txt)

# Check if either adapter is empty
if [[ -z ${adapter_r1} || -z ${adapter_r2} ]]; then
	echo "No adapter detected for ${prefix}, using original FASTQ files for alignment"

	# Use original files
	cp ${r1} ${trimmed_r1}
	cp ${r2} ${trimmed_r2}

else

	# Parse logs to extract adapter sequences
	echo "Detected adapter for R1: $adapter_r1"
	echo "Detected adapter for R2: $adapter_r2"
	
	# FASTQ read adapter trimming with separate adapters
	module load cutadapt
	log_file=${trimmed_dir}/${prefix}_cutadapt.log
	adapter_err_rate=0.2
	n_threads=10

	cutadapt -j ${n_threads} -m 10 --cut 0 -e ${adapter_err_rate} -a ${adapter_r1} -A ${adapter_r2} --pair-filter=both -o ${trimmed_r1} -p ${trimmed_r2} ${r1} ${r2} 2> ${log_file}
	echo "Done with cutadapt for ${prefix}"; date
fi
