module load bedtools/2.30.0
genome_sizes=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cat *F25*.narrowPeak *F49*.narrowPeak *F55*.narrowPeak | cut -f -4 | awk '{OFS="\t"}{print $1,$2,$3,$4"_"NR}' | bedtools sort -i stdin > temp.bed
bedtools slop -i temp.bed -g ${genome_sizes} -b 100 | bedtools merge -i stdin | awk '{print $0"\t"$1"_"$2"_"$3}' | bedtools sort -i stdin > merged.bed
rm temp.bed

module load bedtools/2.30.0
genome_sizes=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cat *S1*.narrowPeak *S3*.narrowPeak *S5*.narrowPeak | cut -f -4 | awk '{OFS="\t"}{print $1,$2,$3,$4"_"NR}' | bedtools sort -i stdin > temp.bed
bedtools slop -i temp.bed -g ${genome_sizes} -b 100 | bedtools merge -i stdin | awk '{print $0"\t"$1"_"$2"_"$3}' | bedtools sort -i stdin > PBS_merged.bed
rm temp.bed

module load bedtools/2.30.0
genome_sizes=/share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes
cat *S2*.narrowPeak *S4*.narrowPeak *S6*.narrowPeak | cut -f -4 | awk '{OFS="\t"}{print $1,$2,$3,$4"_"NR}' | bedtools sort -i stdin > temp.bed
bedtools slop -i temp.bed -g ${genome_sizes} -b 100 | bedtools merge -i stdin | awk '{print $0"\t"$1"_"$2"_"$3}' | bedtools sort -i stdin > IFN_merged.bed
rm temp.bed

