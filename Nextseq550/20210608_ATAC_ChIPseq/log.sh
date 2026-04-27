#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=45000]
#BSUB -q long
#BSUB -W 8:00
#BSUB -o "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210608_ATAC_ChIPseq/%J%I.out"
#BSUB -e "/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210608_ATAC_ChIPseq/%J%I.err"


### merge
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq550/20210608_ATAC_ChIPseq/fastqs

#cd $dir/skin_eQTL_ATAC
#for x in ATAC_F25K_IFN_S2 ATAC_F25K_PBS_S1 ATAC_F49K_IFN_S4 ATAC_F49K_PBS_S3 ATAC_F55K_IFN_S6 ATAC_F55K_PBS_S5;do
#  lib=$(echo ${x} | cut -d"_" -f1,2,3)
#  sample=$(echo $x | cut -d"_" -f2,3,4)
#  echo "lib is " $lib
#  zcat $dir/skin_eQTL_ATAC/${lib}/${sample}_L001_R1_001.fastq.gz $dir/$lib/${sample}_L002_R1_001.fastq.gz $dir/$lib/${sample}_L003_R1_001.fastq.gz $dir/$lib/${sample}_L004_R1_001.fastq.gz > $dir/merged/${lib}_R1.fastq.gz
#  zcat $dir/skin_eQTL_ATAC/${lib}/${sample}_L001_R2_001.fastq.gz $dir/$lib/${sample}_L002_R2_001.fastq.gz $dir/$lib/${sample}_L003_R2_001.fastq.gz $dir/$lib/${sample}_L004_R2_001.fastq.gz > $dir/merged/${lib}_R2.fastq.gz
#  echo "done with "$lib
#done
#
#cd $dir/skin_eQTL_ChIP
#for x in ChIP_F25K_IFN_S8 ChIP_F25K_PBS_S7 ChIP_F49K_IFN_S10 ChIP_F49K_PBS_S9 ChIP_F55K_IFN_S12 ChIP_F55K_PBS_S11 ChIP_Input_K_S13;do
#  lib=$(echo ${x} | cut -d"_" -f1,2,3)
#  sample=$(echo $x | cut -d"_" -f2,3,4)
#  echo "lib is " $lib
#  zcat $dir/skin_eQTL_ChIP/${lib}/${sample}_L001_R1_001.fastq.gz $dir/$lib/${sample}_L002_R1_001.fastq.gz $dir/$lib/${sample}_L003_R1_001.fastq.gz $dir/$lib/${sample}_L004_R1_001.fastq.gz > $dir/merged/${lib}_R1.fastq.gz
#  zcat $dir/skin_eQTL_ChIP/${lib}/${sample}_L001_R2_001.fastq.gz $dir/$lib/${sample}_L002_R2_001.fastq.gz $dir/$lib/${sample}_L003_R2_001.fastq.gz $dir/$lib/${sample}_L004_R2_001.fastq.gz > $dir/merged/${lib}_R2.fastq.gz
#  echo "done with "$lib
#done

#### check md5sum then backup on Amazon
cd $dir/merged
for i in `find . -type f`; do md5sum $i > $i.md5sum; done
cd ..
/project/umw_biocore/bin/amazonBackup.bash merged s3://biocorebackup/garberlab/human/skin/eQTL/Nextseq550/20210608_ATAC_ChIPseq

