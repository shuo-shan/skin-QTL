
module load RSEM/1.2.29
rsem-calculate-expression -p 10 --forward-prob 0 --paired-end --no-bam-output --star --star-path /share/pkg/star/2.5.3a/ --star-gzipped-read-file F60_1.fastq.gz F60_2.fastq.gz /share/data/umw_biocore/dnext_data/genome_data/human/hg38/gencode_v34/RSEM_ref_STAR/genome rsem/F60.rsem

