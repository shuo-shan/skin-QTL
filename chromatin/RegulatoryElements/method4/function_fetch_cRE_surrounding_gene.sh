g=$1
celltype=$2


# set-up
genomeFile=/share/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.chrom.sizes
module load bedtools/2.30.0
dir=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method4
cd ${dir}

# DEFINE TSS from Ensembl GTF file, by /pi/manuel.garber-umw/human/skin/eQTLs/chromatin/log.sh
# columns are: TSS position, strand, Ensembl gene ID + transcript ID, Gene name + transcript ID, gene name, gene biotype
# TSS bed file is created from the start of all transcripts. script: /pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/ensembl2tss.sh
tss=/pi/manuel.garber-umw/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_transcription_start_sites.bed

# find all cREs within 300kbp distance from a gene's TSS
# if a gene has mutliple TSS, the 300kbp flanking distance starts from the start of first TSS and end of last TSS.
awk -v gene=${g} '{OFS="\t"}{if ($5==gene) print $0}' ${tss} | bedtools merge -i stdin | bedtools sort -i stdin > temp_${g}_tss.bed
chr=$(cat temp_${g}_tss.bed | head -1 | cut -f1)
tempstart=$(cat temp_${g}_tss.bed | head -1 | cut -f2)
tempend=$(cat temp_${g}_tss.bed | tail -1 | cut -f3)
echo -e "$chr\t$tempstart\t$tempend" | bedtools slop -b 300000 -i stdin -g ${genomeFile} | sort -k 1,1 -k2,2n > temp_${g}_tss_600kbp_window.bed
bedtools intersect -u -a temp_cRE_${celltype}.bed -b temp_${g}_tss_600kbp_window.bed | awk -v gene=${g} -v ct=${celltype} '{OFS="\t"}{print $1,$2,$3,$4,gene,ct}' > temp_cRE_surrounding_${g}_${celltype}_part0.bed
cat temp_cRE_surrounding_${g}_${celltype}_part0.bed | cut -f1-3 > temp_cRE_surrounding_${g}_${celltype}_part1.bed
cat temp_cRE_surrounding_${g}_${celltype}_part0.bed | cut -f4-6 | sed 's/merged_all_skin-eQTL_ATACseq_files_//g' | sed 's/\t/_/g' | sed 's/peak_/peak/g' > temp_cRE_surrounding_${g}_${celltype}_part2.bed
cat temp_cRE_surrounding_${g}_${celltype}_part0.bed | cut -f4-6 > temp_cRE_surrounding_${g}_${celltype}_part3.bed
paste temp_cRE_surrounding_${g}_${celltype}_part1.bed temp_cRE_surrounding_${g}_${celltype}_part2.bed temp_cRE_surrounding_${g}_${celltype}_part3.bed > temp_cRE_surrounding_${g}_${celltype}.bed

# label cREs as promoters or enhancers. if a cRE overlaps 500bp flanking region of a TSS of a gene, label it as 'promoter'. Otherwise, label as 'enhancer'
# make 500bp flanking region of a gene's TSSs
bedtools slop -b 500 -i temp_${g}_tss.bed -g ${genomeFile} | sort -k 1,1 -k2,2n > temp_${g}_tss_500bp_window.bed
bedtools intersect -u -a temp_cRE_surrounding_${g}_${celltype}.bed -b temp_${g}_tss_500bp_window.bed | awk '{OFS="\t"}{print $0,"promoter"}' >> temp_cRE_surrounding_${g}_PandE_${celltype}.bed
bedtools intersect -v -a temp_cRE_surrounding_${g}_${celltype}.bed -b temp_${g}_tss_500bp_window.bed | awk '{OFS="\t"}{print $0,"enhancer"}' >> temp_cRE_surrounding_${g}_PandE_${celltype}.bed
cat temp_cRE_surrounding_${g}_PandE_${celltype}.bed | sort -k 1,1 -k2,2n >> promoters_and_enhancers_surrounding_genes_${celltype}.bed

# clean-up
rm temp_${g}_tss.bed temp_${g}_tss_600kbp_window.bed temp_cRE_surrounding_${g}_${celltype}_part0.bed temp_cRE_surrounding_${g}_${celltype}_part1.bed temp_cRE_surrounding_${g}_${celltype}_part2.bed temp_cRE_surrounding_${g}_${celltype}_part3.bed temp_cRE_surrounding_${g}_${celltype}.bed temp_${g}_tss_500bp_window.bed temp_cRE_surrounding_${g}_PandE_${celltype}.bed

