#!/bin/bash
# modularized pipeline to demultiplex and barcode modify 24-sample pooled fastq files for celseq2 libraries
# Shuo Shan 05/30/2025

########################################
# set-up
Dir=/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/03312025_3ctk_celseq2_atac_batch1/fastq_celseq2
mkdir -p ${Dir}/log
s=$1 # format like skineQTL-1 to skineQTL-26
cd ${Dir}

######################################### ~20min
## unzip files for demultiplexing
#cd ${Dir}
#echo "unzipping now"; date
#bsub -J unzip_${s}_r1 -W 1:00 -n 16 -R "span[hosts=1]" -R rusage[mem=1000] -q long -e "./log/log_unzip_%J_${s}_R1.err" -o "./log/log_unzip_%J_${s}_R1.out" \
#	bash -c "cd ${Dir}; pigz -p 16 -dc ${s}_R1_001.fastq.gz > ${s}_temp_R1.fastq"
#
#bsub -J unzip_${s}_r2 -W 1:00 -n 16 -R "span[hosts=1]" -R rusage[mem=1000] -q long -e "./log/log_unzip_%J_${s}_R2.err" -o "./log/log_unzip_%J_${s}_R2.out" \
#        bash -c "cd ${Dir}; pigz -p 16 -dc ${s}_R2_001.fastq.gz > ${s}_temp_R2.fastq"
#sleep 30
#while [ "$(bjobs -w | grep unzip_${s} | wc -l)" -gt 0 ];do
#	echo "Waiting for unzipping to finish..."
#	sleep 30
#done
#
######################################### ~60min
## split library into samples in Celseq2 by CS2 barcode
#echo 'starting to demultiplex'${s}; date
#bsub -J demux_${s} -W 3:00 -n 1 -R "span[hosts=1]" -R rusage[mem=8000] -q long -e "${Dir}/log/log_demux_%J_${s}.err" -o "${Dir}/log/log_demux_%J_${s}.out" \
#	bash -c "
#module load openjdk # this loads java
#cd ${Dir}
#
## Demultiplex
#java -Xmx8g -jar /pi/manuel.garber-umw/sshan/scripts/splitter_09.21.15_13.57.jar \
#	F1=${s}_temp_R1.fastq \
#	F2=${s}_temp_R2.fastq \
#	B=NNNNNNSSSSSS \
#	M=${Dir}/barcodeMap/barcode_${s}.txt \
#	HD=1 \
#	O=${s}
#
## Delete temporary files
#rm -f ${s}_temp_R1.fastq ${s}_temp_R2.fastq
#"
#echo 'done demultiplexing for '${s}; date
#
#sleep 30
## don't proceed unless demultiplexing job is done
#while [ "$(bjobs -w | grep -w demux_${s} | wc -l)" -gt 0 ];do
#	echo "Waiting for demultiplexing ${s} to finish..."
#	sleep 60
#done
#
#
######################################### ~5min
## Zip all demultiplexed files
## Use pigz for faster compression
#echo "zipping all demultiplexed fq..."; date
#cd ${Dir}/${s}
#
#find ${Dir}/${s} -name "*.fq" | grep -v unassigned | awk -v Dir=${Dir} -v s=${s} '{print "cd ${Dir}/${s}; xargs -P 8 pigz -f "$0}' > commands_sample.txt
#find ${Dir}/${s} -name "*.fq" | grep unassigned | awk -v Dir=${Dir} -v s=${s} '{print "cd ${Dir}/${s}; xargs -P 16 pigz -f "$0}' > commands_unassigned.txt
#
#while read c; do
#	echo ${c} | bsub -J pigz_${s} -W 1:00 -n 8 -R "span[hosts=1]" -R rusage[mem=100] -q long -e "./%J%I.err" -o "./%J%I.out"
#done < commands_sample.txt
#
#while read c; do
#	echo ${c} | bsub -J pigz_${s} -W 1:00 -n 16 -R "span[hosts=1]" -R rusage[mem=100] -q long -e "./%J%I.err" -o "./%J%I.out"
#done < commands_unassigned.txt
#
## don't proceed unless zipping is done
#while [ "$(bjobs -w | grep -w pigz_${s} | wc -l)" -gt 0 ];do
#	echo "Waiting for zipping to finish..."
#	sleep 10
#done

################
# trim read2 end 50 bases and then modify the barcode to UMI:barcode
echo "trimming read2...for ${s}"; date
cd ${Dir}/${s}

# Generate commands and save to commands.txt
find . -maxdepth 1 -type f -name "*.p2.fq.gz" | awk -v d="$Dir" '{print "bash " d "/trim_and_modify_barcode.sh " substr($0, 3)}' > commands_trimBCmod_${s}.txt
while read c; do
  echo ${c} | bsub -J trimBCmod_${s} -W 8:00 -n 8 -R "span[hosts=1]" -R rusage[mem=1000] -q long -e "${Dir}/log/log_trimBCmod_%J.err" -o "${Dir}/log/log_trimBCmod_%J.out"
done < commands_trimBCmod_${s}.txt

sleep 10
while [ "$(bjobs -w | grep -w trimBCmod_${s} | wc -l)" -gt 0 ];do
	echo "Waiting for trimBCmod to finish..."
	sleep 10
done

# rename to replicate number (S1mod means 1st technical replicate, modified by trimming end 50bp)
repnumber="S1"
for old_name in *.p2_trimmed_bcmodified.fq.gz; do
	sample_id=$(echo "$old_name" | sed 's/\..*//')
	new_name="${sample_id}_3ctk_${repnumber}mod.p2.fq.gz"
	mv ${old_name} ${new_name}
	echo "Renamed: $old_name -> $new_name"
done


