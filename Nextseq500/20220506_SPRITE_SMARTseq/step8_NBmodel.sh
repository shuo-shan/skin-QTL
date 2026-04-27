# haha
# step 1. prepare input files needed to run NB.R
# cluster coverage hg38.5kb.coverage_2-100
cat /share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes | cut -f1 > chroms
while read s; do
  grep -w ${s} ${dir}/hg38.5kb.windows.bed > tmp2
  bedtools coverage -sorted -counts -a tmp2 -b ${dir}/barcoded_elements_sorted.bed >> hg38.5kb.coverage_2-100
  echo ${s}
done < chroms

cat hg38.5kb.coverage_2-100 | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$4}' > hg38.5kb.coverage_2-100.bed

bedtools nuc -fi /share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.fa -bed hg38.5kb.coverage_2-100.bed > tmp
cat tmp | awk '{OFS="\t"}{print $1,$2,$3,$5,$7}' > hg38.5kb.coverage_2-100 # take sprite coverage and GC% for each 5kb tile
rm tmp chroms hg38.5kb.coverage_2-100.bed

# step 2. call NB.R script
mkdir all_cis
