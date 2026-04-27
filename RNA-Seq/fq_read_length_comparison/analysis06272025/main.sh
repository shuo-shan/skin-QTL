# log of commands used to trim reads. reads were previously trimmed 3' end from 150bp to 100bp to prevent mapping to polyA. reads were also barcode modified to cater to ESAT requirements.

Dir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fq_read_length_comparison
cd ${Dir}
module load cutadapt
module load seqtk

### make trimmed 40bp reads from 5'end and also another version from the 3' end. 
### for each sample, downsample reads to match the other version.
outDir=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fq_read_length_comparison/analysis06272025/fastq
# F47_MEL_IFNG
f42=/pi/manuel.garber-umw/human/skin/eQTLs/Nextseq550/20210812_ChIPseq_Celseq2/basespace/fastqs/laneMerged/F47_MEL_IFN_S2.p2.fq.gz
zcat ${f42} | seqtk sample -s100 - 21224030 |  gzip > ${outDir}/42bp_F47_MEL_IFNG_downsampled.p2.fq.gz
f100=/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/03312025_3ctk_celseq2_atac_batch1/fastq_celseq2/skineQTL-9/TGCAGA_F47_MEL_IFNG_3ctk_S1mod.p2.fq.gz
cp ${f100} ${outDir}/100bp_F47_MEL_IFNG.p2.fq.gz
cutadapt --cores=8 -u -58 -o ${outDir}/trimmed5prime42bp_F47_MEL_IFNG.p2.fq.gz ${outDir}/100bp_F47_MEL_IFNG.p2.fq.gz
cutadapt --cores=8 -u 58 -o ${outDir}/trimmed3prime42bp_F47_MEL_IFNG.p2.fq.gz ${outDir}/100bp_F47_MEL_IFNG.p2.fq.gz

# F56_FRB_PBS
f40=/pi/manuel.garber-umw/human/skin/eQTLs/Nextseq2000/20211208_Celseq2_ATACseq/fastq/Celseq2/F56_FRB_PBS_S1.p2.fq.gz
cp ${f40} ${outDir}/orig40bp_F56_FRB_PBS_S1.p2.fq.gz
f100=/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/03312025_3ctk_celseq2_atac_batch1/fastq_celseq2/skineQTL-3_origBC/AGTGTC_F56_FRB_PBS_3ctk_S1mod.p2.fq.gz
zcat ${f100} | seqtk sample -s100 - 18795876 | gzip > ${outDir}/100bp_F56_FRB_PBS_downsampled.p2.fq.gz
cutadapt --cores=8 -u -60 -o ${outDir}/trimmed5prime40bp_F56_FRB_PBS.p2.fq.gz ${outDir}/100bp_F56_FRB_PBS_downsampled.p2.fq.gz
cutadapt --cores=8 -u 60 -o ${outDir}/trimmed3prime40bp_F56_FRB_PBS.p2.fq.gz ${outDir}/100bp_F56_FRB_PBS_downsampled.p2.fq.gz
