while read s
do 
	gunzip ${s}_R1_001.fastq.gz && gunzip ${s}_R2_001.fastq.gz && java -Xmx5g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_R1_001.fastq F2=${s}_R2_001.fastq B=NNNNNNSSSSSS M=barcodeMap HD=1 O=${s} && gzip ${s}_R1_001.fastq && gzip ${s}_R2_001.fastq && echo ${s}; 
done < samples

for i in $(ls */*fq)
do
	gzip $i && echo $i
	echo $i
done
