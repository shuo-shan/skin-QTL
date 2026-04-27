#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20400]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "/pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/%J%I.err"
#BSUB -o "/pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/%J%I.out"
# specify downloaded output destination

outDir=/pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/humanTFChIP
cd ${outDir}

## ENCODE search URL --> metadata plus download URL
#python /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/download_ENCODE_narrowPeak_bed_files.py

# remove genetically modified samples
cat matadata.txt | grep -v modified > filtered_metadata.txt

# download
while read line;do
	ID=$(echo ${line} | cut -d';' -f3)
	target=$(echo ${line} | cut -d';' -f4 | cut -d'/' -f3)
	assembly=$(echo ${line} | cut -d';' -f6)
	url=$(echo ${line} | cut -d';' -f7)

	newName=${ID}_${target}_${assembly}.bed.gz
	wget -O ${newName} ${url}
	echo ${line}"...done"
done < filtered_metadata.txt

# copy the code for later bookkeeping
cp /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/download_ENCODE_narrowPeak_bed_files.sh ${outDir}
cp /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/download_ENCODE_narrowPeak_bed_files.py ${outDir}
mv /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/*.err ${outDir}
mv /pi/manuel.garber-umw/human/skin/eQTLs/literature/ENCODE/*.out ${outDir}
