#!/bin/bash
cd /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO
export PATH=/home/shuo.shan-umw/meme-5.5.1/meme/bin:/home/shuo.shan-umw/meme-5.5.1/meme/libexec/meme-5.5.1:$PATH
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO
mkdir fasta 

#### step 1. prepare fasta sequences for regions of interest as FIMO input
cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_300bp_flanking_window.bed
dict_gene_promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_allcts.bed
dict_gene_enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_allcts.bed
module load bedtools/2.30.0
cat ${dict_gene_promoter} | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' | sort -k4 | uniq --skip-fields=3 > promoter.bed
cat ${dict_gene_enhancer} | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' | sort -k4 | uniq --skip-fields=3 > enhancer.bed

split -l 100 -d promoter.bed promoter.split.
for f in promoter.split.*;do
	bedtools getfasta -fi /share/GHPCC/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa -bed ${f} > ${dir}/fasta/${f}.fa
	echo done with ${f}
done
rm *split*

split -l 100 -d enhancer.bed enhancer.split.
for f in enhancer.split.*;do
	bedtools getfasta -fi /share/GHPCC/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa -bed ${f} > ${f}.fa
	echo done with ${f}
done

#### step 2. prepare motif file into MEME format as FIMO input
cd ${dir}
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
TF=RUNX1
a=$(cat ${memeF} | grep -n ${TF}_ | head -1 | cut -d':' -f1)
b=$(cat ${memeF} | grep -n ${TF}_ | tail -1 | cut -d':' -f1)
cat ${memeF} | head -9 > header
cat ${memeF} | head -n $(expr ${b}) | tail -n $(expr ${b} - ${a} + 1) > body
cat header body > ${TF}.meme
rm header body

#### step 3. run FIMO and fetch top ranked and bottom ranked regions
cd ${dir}
rm -r ${dir}/${TF}
for f in fasta/promoter.split.*.fa;do
	outF=$(echo ${f} | sed 's/fasta\///g')
	fimo --oc ${dir}/${TF} --thresh 1 --no-qvalue ${dir}/${TF}.meme ${dir}/${f}
	cd ${dir}/${TF}
	mv fimo.tsv fimo.${outF}.tsv
done
cd ${dir}/${TF}
cat fimo.*.tsv | grep chr | cut -f3,4,5,6,7,8 | sort -r -n -k5 > fimo.txt
cat fimo.txt | cut -f1 > temp1
cat fimo.txt | cut -f2,3,4,5,6 > temp2
cat temp1 | sed 's/:/\t/g' | sed 's/-/\t/g' > tempa
paste tempa temp2 > fimo.bed
cat fimo.bed | awk '{if ($8 < 0.0001) print $0}' > fimo_topRanked_promoters.bed
cat fimo.bed | awk '{if ($8 > 0.0001) print $0}' | tail -n 1000 > fimo_bottomRanked_promoters.bed
rm *split*.tsv *xml fimo.gff fimo.html temp1 temp2 tempa

#### step 4. generalized function to find motif center and 1kb flanking region, make bed files for each
function get_motif_center_and_flanking_region {
	inF=$1
	inFname=$(echo ${inF} | sed 's/.*\///g' | sed 's/\.bed//g')
	cat ${inF} | awk '{if ($6=="+") print $0}' > temp_${inFname}_plus.bed
	cat temp_${inFname}_plus.bed | awk '{OFS="\t"}{print $1,$2+$4+(($5-$4)/2-1),$2+$4+(($5-$4)/2),$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_plus.bed
        cat temp_${inFname}_plus.bed | awk '{OFS="\t"}{print $1,$2+$4+(($5-$4)/2-1)-500,$2+$4+(($5-$4)/2)+500,$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_flanking_1kb_plus.bed	
	cat ${inF} | awk '{if ($6=="-") print $0}' > temp_${inFname}_minus.bed
	cat temp_${inFname}_minus.bed | awk '{OFS="\t"}{print $1,$2+$4+(($5-$4)/2)-1,$2+$4+(($5-$4)/2),$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_minus.bed
	cat temp_${inFname}_minus.bed | awk '{OFS="\t"}{print $1,$2+$4+(($5-$4)/2)-1-500,$2+$4+(($5-$4)/2)+500,$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_flanking_1kb_minus.bed
	cat temp_${inFname}_motif_center_plus.bed temp_${inFname}_motif_center_minus.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$6,$5}' | sort -n -k5 > ${inFname}_motif_center.bed
	cat temp_${inFname}_motif_center_flanking_1kb_plus.bed temp_${inFname}_motif_center_flanking_1kb_minus.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$6,$5}' | sort -n -k5 > ${inFname}_motif_center_flanking_1kb.bed
	rm temp_*.bed; echo "done with processing "${inF}
}
get_motif_center_and_flanking_region fimo_topRanked_promoters.bed
get_motif_center_and_flanking_region fimo_bottomRanked_promoters.bed

#### step 5. run FIMO to test RUNX1 motif top/bottom matching regions on GATA2 motif
cd ${dir}/${TF}
memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
TF2=GATA2
#regionF=fimo_topRanked_promoters_motif_center_flanking_1kb.bed
regionF=fimo_bottomRanked_promoters_motif_center_flanking_1kb.bed
regionTag=$(echo ${regionF} | sed 's/\.bed//g')

# fetch TF.meme motif file
a=$(cat ${memeF} | grep -n ${TF2}_ | grep H11MO.0 | head -1 | cut -d':' -f1)
b=$(cat ${memeF} | grep -n ${TF2}_ | grep H11MO.0 | tail -1 | cut -d':' -f1)
cat ${memeF} | head -9 > header
cat ${memeF} | head -n $(expr ${b}) | tail -n $(expr ${b} - ${a} + 1) > body
cat header body > ${TF2}.meme
rm header body

# make fasta file of regions of interest for fimo
rm -r fasta
mkdir fasta
bed_to_split=${regionF}
split -l 50 -d ${bed_to_split} split.${bed_to_split}
for f in split.*;do
        bedtools getfasta -fi /share/GHPCC/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa -bed ${f} > fasta/${f}.fa
        echo done with ${f}
done
rm split*

# run fimo
cd ${dir}/${TF}
rm -r ${dir}/${TF}/${TF2}
for f in ${dir}/${TF}/fasta/split*.fa;do
	outF=$(echo ${f} | sed 's/.*fasta\///g')
	fimo --oc ${dir}/${TF}/${TF2} --thresh 1 --no-qvalue ${dir}/${TF}/${TF2}.meme ${f}
	cd ${dir}/${TF}/${TF2}
	mv fimo.tsv fimo.${outF}.tsv
done

# compile fimo output and get motif center regions.
cd ${dir}/${TF}/${TF2}
cat fimo.*.tsv | grep chr | cut -f3,4,5,6,7,8 | sort -r -n -k5 > fimo.txt
cat fimo.txt | cut -f1 > temp1
cat fimo.txt | cut -f2,3,4,5,6 > temp2
cat temp1 | sed 's/:/\t/g' | sed 's/-/\t/g' > tempa
paste tempa temp2 > fimo.bed
cat fimo.bed | awk '{if ($8 < 0.0001) print $0}' > fimo_topRanked_regions.bed
rm *split*.tsv *xml fimo.gff fimo.html temp1 temp2 tempa
get_motif_center_and_flanking_region fimo_topRanked_regions.bed
mv fimo_topRanked_regions_motif_center.bed ${dir}/${TF}/${TF2}_motif_centers_topRanked_in_${TF}_${regionTag}.bed
rm -r ${dir}/${TF}/${TF2}
rm -r ${dir}/${TF}/fasta

# use R to calculate distance between GATA2 motif center and RUNX1 motif center







