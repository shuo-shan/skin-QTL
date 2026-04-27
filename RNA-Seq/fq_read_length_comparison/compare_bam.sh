# compare the bam files of 100bp vs 40bp
module load samtools
dir40=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fq_read_length_comparison/DolphinNext/40bp/report10990/star
f40=${dir40}/trimmed40bp_CTTCGA_F114_FRB_PBS_3ctk_S1mod_sorted.bam

dir100=/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/fq_read_length_comparison/DolphinNext/100bp/report10991/star
f100=${dir100}/CTTCGA_F114_FRB_PBS_3ctk_S1mod_sorted.bam

# check read alignment rate
samtools flagstat ${f100}
samtools flagstat ${f40}
