#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=100000]
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

for expt in SRR400268 SRR400267 SRR400266 SRR400265 SRR400264; do
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
done &
wait
##### step 3: process bam files
#bamdir=/nl/umw_manuel_garber/human/skin/eQTLs/literature/fibroblast/dixon_lung_fibroblast_imr90_hic/report6435/merged_bams
#samtools merge ${outputname}.bam ${bamdir}/SRR400264_sorted.bam  ${bamdir}/SRR400265_sorted.bam  ${bamdir}/SRR400266_sorted.bam  ${bamdir}/SRR400267_sorted.bam  ${bamdir}/SRR400268_sorted.bam
#samtools sort -@ 5 ${outputname}.bam | samtools view -@ 5 - > temp1
#echo "done sorting merged bam";date
#cat temp1 | ${lib}/remove_dup_PE_SAM_sorted.pl > ${outputname}_deduped.sam
#wait
#samtools view -@ 5 -bS -t ${hg38Ind} -o ${outputname}.sorted.nodup.bam ${outputname}_deduped.sam
#echo "done removing duplicates for merged bam";date
#samtools view ${outputname}.sorted.nodup.bam | cut -f2-8 | $lib/bam_to_temp_HiC.pl > ${outputname}.temp
#wait
#echo "merged and de-duplicated bam files";date
#echo "done with step 3";date

#### step 4: map reads pair to fragment pairs, 36 is the read length for mapping
#echo "working now...";date
#${lib}/reads_2_cis_frag_loop.pl ${bed} 36 ${outputname}.loop.inward ${outputname}.loop.outward ${outputname}.loop.samestrand summary.frag_loop.read_count ${outputname} ${outputname}.temp
#echo "done with reads_2_cis_frag_loop";date
#${lib}/reads_2_trans_frag_loop.pl ${bed} 36 ${outputname}.loop.trans ${outputname}.temp
#echo "done with reads_2_trans_frag_loop";date
#wait
#echo "done with step 4.1";date
#perl ${lib}/summary_sorted_frag_loop.pl ${bed} ${dir}/${outputname}.loop.inward >> ${dir}/temp.${outputname}.loop.inward
#echo "done with step 4.2";date
#wait
#perl ${lib}/summary_sorted_frag_loop.pl ${bed} ${dir}/${outputname}.loop.outward >> ${dir}/temp.${outputname}.loop.outward
#wait
#perl ${lib}/summary_sorted_frag_loop.pl ${bed} ${dir}/${outputname}.loop.samestrand >> ${dir}/temp.${outputname}.loop.samestrand
#wait
#perl ${lib}/summary_sorted_trans_frag_loop.pl ${dir}/${outputname}.loop.trans >> ${dir}/temp.${outputname}.loop.trans
#wait
#mv temp.${outputname}.loop.inward ${outputname}.loop.inward
#mv temp.${outputname}.loop.outward ${outputname}.loop.outward
#mv temp.${outputname}.loop.samestrand ${outputname}.loop.samestrand
#mv temp.${outputname}.loop.trans ${outputname}.loop.trans
#echo "done with step 4.2";date

#perl ${lib}/resort_by_frag_id.pl ${bed} ${outputname}.loop.inward
#perl ${lib}/resort_by_frag_id.pl ${bed} ${outputname}.loop.outward
#wait
#perl ${lib}/resort_by_frag_id.pl ${bed} ${outputname}.loop.samestrand
#wait
#echo "done with step 4.3";date
#for file in ${outputname}.loop.inward ${outputname}.loop.outward ${outputname}.loop.samestrand;do
#        ${lib}/resort_by_frag_id.pl ${bed} temp.${file} &
#done
#cat ${outputname}.loop.trans | ${lib}/summary_sorted_trans_frag_loop.pl - > temp.${outputname}.loop.trans
#wait
#echo "done with step 4.3";date
#${lib}/merge_sorted_frag_loop.pl ${outputname}.loop.samestrand > frag_loop.${outputname}.samestrand
#wait
#${lib}/merge_sorted_frag_loop.pl <(cat ${outputname}.loop.inward | awk '{if($4>1000)print $0}') > frag_loop.${outputname}.inward 
#wait
#${lib}/merge_sorted_frag_loop.pl <(cat ${outputname}.loop.outward | awk '{if($4>25000)print $0}') > frag_loop.${outputname}.outward 
#wait
#echo "done with step 4.4";date
#${lib}/merge_sorted_frag_loop.pl frag_loop.${outputname}.samestrand frag_loop.${outputname}.inward frag_loop.${outputname}.outward > frag_loop.${outputname}.cis 
#wait
#${lib}/merge_sorted_frag_loop.pl ${outputname}.loop.trans > frag_loop.${outputname}.trans 
#wait
#echo "done with step 4.5";date
#frag_loop.$outputname.cis and frag_loop.$outputname.trans will be used to run HiCorr
#
#### step 5: run HiCorr
#echo "running HiCorr...";date
#${HiCorrPath}/HiCorr HindIII frag_loop.${outputname}.cis frag_loop.${outputname}.trans ${outputname} hg38
#echo "done running HiCorr, yaay!!";date
#
#
#### clean-up
#mv ./%J%I.err ${dir}
#mv ./%J%I.out ${dir}
















#
#
#
#
