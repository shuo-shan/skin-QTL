module load bedtools/2.30.0
genome_sizes=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes

cat *F25_F*.narrowPeak *F49_F*.narrowPeak *F55_F*.narrowPeak | cut -f -4 | awk '{OFS="\t"}{print $1,$2,$3,$4"_"NR}' | bedtools sort -i stdin > temp.bed
bedtools slop -i temp.bed -g ${genome_sizes} -b 100 | bedtools merge -i stdin | awk '{print $0"\t"$1"_"$2"_"$3}' | bedtools sort -i stdin > ATAC_FRB_merged.bed

cat *F25_K*.narrowPeak *F49_K*.narrowPeak *F55_K*.narrowPeak | cut -f -4 | awk '{OFS="\t"}{print $1,$2,$3,$4"_"NR}' | bedtools sort -i stdin > temp.bed
bedtools slop -i temp.bed -g ${genome_sizes} -b 100 | bedtools merge -i stdin | awk '{print $0"\t"$1"_"$2"_"$3}' | bedtools sort -i stdin > ATAC_KRT_merged.bed

cat *F25_M*.narrowPeak *F49_M*.narrowPeak *F55_M*.narrowPeak | cut -f -4 | awk '{OFS="\t"}{print $1,$2,$3,$4"_"NR}' | bedtools sort -i stdin > temp.bed
bedtools slop -i temp.bed -g ${genome_sizes} -b 100 | bedtools merge -i stdin | awk '{print $0"\t"$1"_"$2"_"$3}' | bedtools sort -i stdin > ATAC_MEL_merged.bed

rm temp.bed
