#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=1020]
#BSUB -q long
#BSUB -W 48:00
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/01152025_shallowtest/%J%I.demultiplex_trim_modify.out"
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/01152025_shallowtest/%J%I.demultiplex_trim_modify.err"

########################################
# split library into samples in Celseq2 by CS2 barcode
module load openjdk # this loads java
Dir=/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/01152025_shallowtest
barcodeMap=${Dir}/barcodeMap/barcode.txt
s=skineQTL01152025shallowseq
cd ${Dir}/fastq

## Set trap for cleanup in case of error
#trap 'rm -f ${s}_temp_R1.fastq ${s}_temp_R2.fastq' EXIT
#
## Decompress fastq to new temporary file
#echo "decompressing...";date
#gunzip -c ${s}_R1_001.fastq.gz > ${s}_temp_R1.fastq
#gunzip -c ${s}_R2_001.fastq.gz > ${s}_temp_R2.fastq
## Demultiplex
#echo "demultiplexing...";date
#java -Xmx100g -jar /pi/manuel.garber-umw/sshan/scripts/splitter_09.21.15_13.57.jar \
#	F1=${s}_temp_R1.fastq \
#	F2=${s}_temp_R2.fastq \
#	B=NNNNNNSSSSSS \
#	M=$barcodeMap \
#	HD=1 \
#	O=${s}
#
## Delete temporary files
#echo "removing temp files..."; date
#rm -f ${s}_temp_R1.fastq ${s}_temp_R2.fastq
#
## Zip all demultiplexed files
## Use pigz for faster compression
#echo "zipping all demultiplexed fq..."; date
#cd ${Dir}/fastq/${s}
#find . -name "*.fq" | xargs -P 4 pigz -f
## find . -name "*.fq" | xargs -P 4 gzip -f (use this if pigz isn't available.)


#################
## trim read2 end 50 bases and then modify the barcode to UMI:barcode
#cd ${Dir}/fastq/${s}
## Generate commands and save to commands.txt
#find . -maxdepth 1 -type f -name "*.p2.fq.gz" | awk -v d="$Dir" '{print "bash " d "/trim_and_modify_barcode.sh " substr($0, 3)}' > ${Dir}/commands.txt
#while read c; do
#  echo ${c} | bsub -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=10200] -q long -e "./%J%I.err" -o "./%J%I.out" -J trimNmod
#done < ${Dir}/commands.txt
#
## clean-up
#mkdir trimmed_barcode_modified_fastq
#mv *trimmed_bcmodified.fq.gz trimmed_barcode_modified_fastq
#
# rename to replicate number (S1mod means 1st technical replicate, modified by trimming end 50bp)
cd ${Dir}/fastq/${s}/trimmed_barcode_modified_fastq
repnumber="S1"
for old_name in *.p2_trimmed_bcmodified.fq.gz; do
	sample_id=$(echo "$old_name" | sed 's/\..*//')
	new_name="${sample_id}_3ctk_${repnumber}mod.p2.fq.gz"
	mv ${old_name} ${new_name}
	echo "Renamed: $old_name -> $new_name"
done

