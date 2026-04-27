#!/bin/bash
#BSUB -n 6
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./%J%I.err"
#BSUB -o "./%J%I.out"

# https://github.com/shanshan950/Hi-C-data-preprocess/blob/master/documents/Fastq-to-FragmentContact.Tissue_example.md
module load condas/2018-05-11
source activate sshan_isoform

module load bowtie/1.3.0
module load samtools/1.9
module load perl/5.28.1
module load python3/3.5.0
module load python3/3.5.0_packages/numpy/1.18.5

dir=/nl/umw_manuel_garber/human/skin/eQTLs/literature/fibroblast/dixon_lung_fibroblast_imr90_hic/hg19
cd ${dir}

#### step 1: download data
##for file in SRR400268 SRR400267 SRR400266 SRR400265 SRR400264; do fastq-dump --split-files $SRR &;done wait
## check for the shortest read length in files
#for file in $(ls *_1.fastq);do 
#  echo $file $(cat $file | head -2 | tail -1 | wc -c)
#done # Different read length, therefore using 36bp for mapping for a fair processing
## We chose 36bp as read length because the shortest read length for this example is 36bp.

### step 2: bowtie mapping
#hg19=Your_hg19_BowtieIndexPath/YourIndexPrefix
#hg19Ind=Your_PathTo_hg19.fa.fai
#HiCorrPath=<where you put HiCorr>
#lib=Path_to_lib
#bed=Path_to_fragbed # <chr> <start> <end> <frag_id>
#outputname=Adrenal
hg19=/share/data/umw_biocore/genome_data/human/hg19/hg19
hg19Ind=/share/data/umw_biocore/genome_data/human/hg19/hg19.fa.fai
HiCorrPath=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/HiCorr/HiCorr
lib=/nl/umw_manuel_garber/sshan/scripts/hic_preprocess
bed=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/HiCorr/HiCorr/ref/HindIII/hg19.HindIII.frag.bed # <chr> <start> <end> <frag_id>
outputname=IMR90
file=SRR400264 # SRR400268 SRR400267 SRR400266 SRR400265 SRR400264

for expt in ${file};do
  fq1=${expt}_1.fastq
  fq2=${expt}_2.fastq
  length=$(head ${fq1} | tail -1 | wc -m)
  let length=${length}-1
  let trlen=${length}-36
  bowtie -v 3 -m 1 --trim3 ${trlen} --best --strata --time -p 5 --sam ${hg19} ${fq1} ${expt}.R1.sam &
  bowtie -v 3 -m 1 --trim3 ${trlen} --best --strata --time -p 5 --sam ${hg19} ${fq2} ${expt}.R2.sam &
  wait
  echo Total reads count for $expt is $(samtools view ${expt}.R1 | grep -vE ^@ | wc -l) >> ${expt}.summary.total.read_count &
  samtools view -u ${expt}.R1.sam | samtools sort -@ 5 -n -T ${expt}.R1 -o ${expt}.R1.sorted.bam &
  samtools view -u ${expt}.R2.sam | samtools sort -@ 5 -n -T ${expt}.R2 -o ${expt}.R2.sorted.bam &
  wait
  $lib/pairing_two_SAM_reads.pl <(samtools view ${expt}.R1.sorted.bam) <(samtools view ${expt}.R2.sorted.bam) | samtools view -bS -t $hg19Ind -o - - > ${expt}.bam
  echo "done pairing two SAM reads for "${expt}; date;
done
