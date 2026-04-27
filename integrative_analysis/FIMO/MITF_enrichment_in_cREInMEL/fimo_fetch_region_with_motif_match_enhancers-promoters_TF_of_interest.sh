#!/bin/bash
# NOTE TO SELF: don't forget to run the functions at the bottom of this file!
export PATH=/home/shuo.shan-umw/meme-5.5.1/meme/bin:/home/shuo.shan-umw/meme-5.5.1/meme/libexec/meme-5.5.1:$PATH
module load bedtools/2.30.0

#### resources
celltype=MEL
cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_1000bp_flanking_window.bed
list_gene_promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_${celltype}.txt
list_gene_enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_${celltype}.txt

# set-up
TF=MITF
regionTag=cREIn${celltype}
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/${TF}_enrichment_in_${regionTag}
mkdir ${dir}
cd ${dir}
cat ${list_gene_promoter} ${list_gene_enhancer} > list_gene_cREs.txt
regionF=${dir}/list_gene_cREs.txt

# step 1. generate bed file of regions of interest
cut -f2 ${regionF} > ${regionTag}.txt
cat ${cRE_bed} | grep -v chrX | grep -v chrY > autosome_cRE.bed
bash /pi/manuel.garber-umw/sshan/scripts/function_rapid_fgrep.sh ${dir}/${regionTag}.txt ${dir}/autosome_cRE.bed ${regionTag}.bed rapidGrep 5000

# step 2. run FIMO on the TF-of-interest for its motif enrichment against regions-of-interest
run_FIMO ${dir} ${TF} ${regionTag}.bed 100
mv ${dir}/${TF}/fimo_${TF}.bed ./

# step 3. find enhancers/promoters that intersect with the fimo matching region
bedtools intersect -wa -wb -a ${regionTag}.bed -b fimo_${TF}.bed | sort -f -k4 > ${regionTag}_with_${TF}_fimo_score.txt 
cat ${regionTag}_with_${TF}_fimo_score.txt | cut -f4 | sort | uniq > ${regionTag}_all.txt
grep -w -f ${regionTag}_all.txt ${regionF} | sort -f -k2 > genes_regulated_by_${regionTag}_with_${TF}_all.txt 
join -i -1 4 -2 2 -o 1.1,1.2,1.3,1.4,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.1,2.3 ${regionTag}_with_${TF}_fimo_score.txt genes_regulated_by_${regionTag}_with_${TF}_all.txt > temp

cat ${list_gene_promoter} | cut -f2 | sort -f | uniq > tempf1
grep -w -f tempf1 temp | awk '{OFS="\t"}{print $0,"promoter"}'> temp.promoters
cat ${list_gene_enhancer} | cut -f2 | sort -f | uniq > tempf2
grep -w -f tempf2 temp | awk '{OFS="\t"}{print $0,"enhancer"}'> temp.enhancers
cat temp.promoters temp.enhancers | sort -g -k12 | tr ' ' '\t' > ${regionTag}_with_${TF}_fimo_score_and_linked_gene.txt
cat temp.promoters temp.enhancers | sort -g -k12 | awk '!seen[$4]++' | tr ' ' '\t' > ${regionTag}_with_best_${TF}_fimo_score_and_linked_gene.txt

# step 4. compare with known MITF binding sites from primary human melanocyte MITF ChIPseq
cat ${regionTag}_with_${TF}_fimo_score_and_linked_gene.txt | cut -f1-4 | bedtools sort -i stdin | awk '!seen[$4]++'> temp.sorted.bed
bedtools intersect -loj -a /pi/manuel.garber-umw/human/skin/eQTLs/literature/MITF_CHIPseq/GSE50681_MITF_chipseq_peaks_hg38.bed \
                       -b temp.sorted.bed > MITF_ChIPseq_peaks_intersecting_with_fimo_predicted_${TF}_binding_sites_in_${regionTag}.bed	
bedtools intersect -loj -a temp.sorted.bed -b /pi/manuel.garber-umw/human/skin/eQTLs/literature/MITF_CHIPseq/GSE50681_MITF_chipseq_peaks_hg38.bed >\
                       fimo_predicted_${TF}_binding_sites_in_${regionTag}_intersecting_with_MITF_ChIPseq_peaks.bed	


# clean-up
rm -r ${TF}
rm -r rapid_fgrep_temp
rm ${regionTag}_top.txt ${regionTag}.txt ${regionTag}.bed 
rm ${regionTag}_with_${TF}_fimo_score.txt ${regionTag}_all.txt genes_regulated_by_${regionTag}_with_${TF}_all.txt
rm temp* *.bed


#### function 1. prepare motif file into MEME format as FIMO input
function get_meme_for_TF {
	echo "preparing HOCOMOCO motif file into MEME format as FIMO input"
	thisTF=$1
	memeF=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
	a=$(cat ${memeF} | grep -n ${thisTF}_ | grep H11MO.0 | head -1 | cut -d':' -f1)
	b=$(cat ${memeF} | grep -n ${thisTF}_ | grep H11MO.0 | tail -1 | cut -d':' -f1)
	cat ${memeF} | head -9 > header
	cat ${memeF} | head -n $(expr ${b}) | tail -n $(expr ${b} - ${a} + 1) > body
	cat header body > ${thisTF}.meme
	rm header body
}

#### function 2. prepare fasta sequences for regions of interest as FIMO input
function bed2fasta {
	# this function takes a bed file, split it, then convert to hg38 fasta file for each split bed in a ./fasta folder.
	echo "preparing fasta sequences for regions of interest as FIMO input"
	rm -r ./fasta
	mkdir ./fasta
	inF=$1 # a bed file to split and make fasta files
	splitn=${2:-100} # number of lines to split to. default is 100.
	split -l ${splitn} -d ${inF} split.regions.
	for f in split.regions.*;do
		bedtools getfasta -fi /share/GHPCC/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa -bed ${f} > ./fasta/${f}.fa
		echo done with ${f}
	done
	rm *split*
}

#### function 3
function get_motif_center_and_flanking_region {
        echo "generalized function to find motif center and 1kb flanking region, make bed files for each"
        inF=$1
        inFname=$(echo ${inF} | sed 's/.*\///g' | sed 's/\.bed//g')
        cat ${inF} | awk '{if ($6=="+") print $0}' > temp_${inFname}_plus.bed
        cat temp_${inFname}_plus.bed | awk '{OFS="\t"}{OFMT="%.0f";print $1,$2+$4+(($5-$4)/2-1),$2+$4+(($5-$4)/2),$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_plus.bed
        cat temp_${inFname}_plus.bed | awk '{OFS="\t"}{OFMT="%.0f";print $1,$2+$4+(($5-$4)/2-1)-500,$2+$4+(($5-$4)/2)+500,$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_flanking_1kb_plus.bed
        cat ${inF} | awk '{if ($6=="-") print $0}' > temp_${inFname}_minus.bed
        cat temp_${inFname}_minus.bed | awk '{OFS="\t"}{OFMT="%.0f";print $1,$2+$4+(($5-$4)/2)-1,$2+$4+(($5-$4)/2),$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_minus.bed
        cat temp_${inFname}_minus.bed | awk '{OFS="\t"}{OFMT="%.0f";print $1,$2+$4+(($5-$4)/2)-1-500,$2+$4+(($5-$4)/2)+500,$1"_"$2"_"$3"_"$4"_"$5,$6,$7,$8}' > temp_${inFname}_motif_center_flanking_1kb_minus.bed
        cat temp_${inFname}_motif_center_plus.bed temp_${inFname}_motif_center_minus.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$6,$5}' | sort -n -k5 > ${inFname}_motif_center.bed
        cat temp_${inFname}_motif_center_flanking_1kb_plus.bed temp_${inFname}_motif_center_flanking_1kb_minus.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$6,$5}' | sort -n -k5 > ${inFname}_motif_center_flanking_1kb.bed
        rm temp_*.bed; echo "done with processing "${inF}
}

#### function 4. run FIMO and fetch top ranked regions
function run_FIMO {
	echo "takes in TF name and bed file to query, outputs a fimo_${TF}.bed with all bed file entries and matching scores"
	thisDir=$1
	thisTF=$2
	thisRegions=$3
	thisRegionsSplitN=$4
	thisRegionTag=$5
	export PATH=/home/shuo.shan-umw/meme-5.5.1/meme/bin:/home/shuo.shan-umw/meme-5.5.1/meme/libexec/meme-5.5.1:$PATH
	module load bedtools/2.30.0
	get_meme_for_TF ${thisTF}                       # run custom function
	bed2fasta ${thisRegions} ${thisRegionsSplitN}   # run custom function
	cd ${thisDir}
	rm -r ${thisDir}/${thisTF}
	for f in fasta/split.*.fa;do
		outF=$(echo ${f} | sed 's/fasta\///g')
		fimo --oc ${thisDir}/${thisTF} --thresh 0.1 --no-qvalue ${thisDir}/${thisTF}.meme ${thisDir}/${f}
		cd ${thisDir}/${thisTF}
		mv fimo.tsv fimo.${outF}.tsv
	done
	cd ${thisDir}/${thisTF}
	cat fimo.*.tsv | grep chr | cut -f3,4,5,6,7,8,10 | sort -r -n -k5 > fimo.txt
	cat fimo.txt | cut -f1 > temp1
	cat fimo.txt | cut -f2,3,4,5,6,7 > temp2
	cat temp1 | sed 's/:/\t/g' | sed 's/-/\t/g' > tempa
	paste tempa temp2 | sort -r -n -k7 | awk '!seen[$1"_"$2"_"$3]++' > fimo_${thisTF}.bed
	mv fimo.txt fimo_output_${thisTF}.txt
	#rm *split*.tsv *xml fimo.gff fimo.html temp1 temp2 tempa fimo.bed 
	rm -r ${thisDir}/fasta
	cd ${thisDir}
}


#### function 5. calculate motif center distance
function calculate_motif_center_distance {
	TF2center=fimo_topRanked_${TF}_motif_center.bed 
	TF1flank=regions.bed

}




