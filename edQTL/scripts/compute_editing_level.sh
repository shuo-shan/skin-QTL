#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=2040]
#BSUB -q long
#BSUB -W 08:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

donorName=$1 # e.g. F25_gencove_hg38

# set-up working directory
module load samtools/1.16.1
dir=/pi/manuel.garber-umw/human/skin/eQTLs/edQTL
cd ${dir}/output
mkdir ${donorName}
cd ${donorName}

# merge all RNA-seq bam files for this donor
#bamdir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/bam
bamdir1=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/melanocytes/pipeline_12022022/DolphinNext/report7615/star
bamdir2=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fibroblasts/pipeline_11192022/DolphinNext/report7561/star
bamdir3=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/keratinocytes/pipeline_11192022/DolphinNext/report7559/star
ls ${bamdir1} | grep ${donorName}_ | grep .bam | grep -v .bam.bai | awk -v dir=${bamdir1} '{print dir"/"$0}' > bamfiles.txt
ls ${bamdir2} | grep ${donorName}_ | grep .bam | grep -v .bam.bai | awk -v dir=${bamdir2} '{print dir"/"$0}' >> bamfiles.txt
ls ${bamdir3} | grep ${donorName}_ | grep .bam | grep -v .bam.bai | awk -v dir=${bamdir3} '{print dir"/"$0}' >> bamfiles.txt
samtools merge -o ${donorName}.bam -b bamfiles.txt
samtools index ${donorName}.bam

# run the script
editing_sites=${dir}/data/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed
bamFile=${donorName}.bam
outName=${donorName}.out
perl ${dir}/scripts/parse_pileup_query.pl ${editing_sites} ${bamFile} ${outName} # need to provide 3 input:Edit Site list, INDEXED BAM alignment file and output file name

# clean-up
cp ${outName} ${dir}/output
cd ${dir}/output
#rm -r ${donorName}
