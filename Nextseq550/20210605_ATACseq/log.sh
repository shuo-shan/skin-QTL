#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=45000]
#BSUB -q long
#BSUB -W 8:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210605_ATACseq/%J%I.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210605_ATACseq/%J%I.err"


#module load bcl2fastq/2.17.1.14
#
## 06_05_21 sequencing data with Jake, Ken, Wei
#indir=/nl/umw_manuel_garber/kgellatly/10X/cellranger_6/6_5_21/Files
#samplesheet=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210605_ATACseq/SampleSheet.csv
#outdir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210605_ATACseq/fastqs
#
#bcl2fastq \
#--runfolder-dir $indir \
#--sample-sheet $samplesheet \
#--output-dir $outdir \
#--interop-dir $outdir \
#--use-bases-mask Y*,I8n*,n*,Y* \
#--mask-short-adapter-reads 0 \
#--minimum-trimmed-read-length 0 \
#--barcode-mismatches 1
#
#
#### check md5sum then backup on Amazon
#cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210605_ATACseq/fastqs
#for i in `find . -type f`; do md5sum $i > $i.md5sum; done
#cd /nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210605_ATACseq/
#/project/umw_biocore/bin/amazonBackup.bash fastqs s3://biocorebackup/garberlab/human/skin/eQTL/Nextseq550/20210605_ATACseq

### merge
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210605_ATACseq/fastqs
for sample in F25M_PBS_S1 F25M_IFN_S2 F49M_PBS_S3 F49M_IFN_S4 F55M_PBS_S5 F55M_IFN_S6;do
  lib=$(echo $sample | cut -d"_" -f1,2)
  echo "lib is " $lib
  zcat $dir/$lib/${sample}_L001_R1_001.fastq.gz $dir/$lib/${sample}_L002_R1_001.fastq.gz $dir/$lib/${sample}_L003_R1_001.fastq.gz $dir/$lib/${sample}_L004_R1_001.fastq.gz > $dir/merged/${lib}_R1.fastq.gz
  zcat $dir/$lib/${sample}_L001_R2_001.fastq.gz $dir/$lib/${sample}_L002_R2_001.fastq.gz $dir/$lib/${sample}_L003_R2_001.fastq.gz $dir/$lib/${sample}_L004_R2_001.fastq.gz > $dir/merged/${lib}_R2.fastq.gz
  echo "done with "$lib
done
