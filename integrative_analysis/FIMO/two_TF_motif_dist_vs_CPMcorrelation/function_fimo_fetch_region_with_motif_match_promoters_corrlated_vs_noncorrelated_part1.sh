#!/bin/bash
export PATH=/home/shuo.shan-umw/meme-5.5.1/meme/bin:/home/shuo.shan-umw/meme-5.5.1/meme/libexec/meme-5.5.1:$PATH
module load bedtools/2.30.0
TF1=$1

#### resources
cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_300bp_flanking_window.bed
dict_gene_promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_allcts.bed
dict_gene_enhancer=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_closest_gene_enhancer_links_allcts.bed
source activate fastQTL

#### set-up
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/${TF1}
mkdir ${dir}
cd ${dir}

#### functions
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
		fimo --oc ${thisDir}/${thisTF} --thresh 1 --no-qvalue ${thisDir}/${thisTF}.meme ${thisDir}/${f}
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


#### function 5. run FIMO for TF2 in TF1 regions.
function run_FIMO_for_TF2 {
        echo "takes in TF name and bed file to query, outputs a fimo_${TF}.bed with all bed file entries and matching scores"
        thisDir=$1
        thisTF=$2
        thisRegions=$3
        thisRegionsSplitN=$4
        export PATH=/home/shuo.shan-umw/meme-5.5.1/meme/bin:/home/shuo.shan-umw/meme-5.5.1/meme/libexec/meme-5.5.1:$PATH
        module load bedtools/2.30.0
        get_meme_for_TF ${thisTF}                       # run custom function
        bed2fasta ${thisRegions} ${thisRegionsSplitN}   # run custom function
        for f in fasta/split.*.fa;do
                outF=$(echo ${f} | sed 's/fasta\///g')
                fimo --oc ${thisDir} --thresh 1 --no-qvalue ${thisDir}/${thisTF}.meme ${thisDir}/${f}
                mv fimo.tsv fimo.${outF}.tsv
        done
        cat fimo.*.tsv | grep chr | cut -f3,4,5,6,7,8,10 | sort -r -n -k5 > fimo.txt
        cat fimo.txt | cut -f1 > temp1
        cat fimo.txt | cut -f2,3,4,5,6,7 > temp2
        cat temp1 | sed 's/:/\t/g' | sed 's/-/\t/g' > tempa
        paste tempa temp2 | sort -r -n -k7 | awk '!seen[$1"_"$2"_"$3]++' > fimo_${thisTF}.bed
        mv fimo.txt fimo_output_${thisTF}.txt
        rm *split*.tsv *xml fimo.gff fimo.html temp1 temp2 tempa ${thisTF}.meme fimo_output_${thisTF}.txt
        rm -r ${thisDir}/fasta
}
#### __MAIN__
###  For the target gene, find highly and lowly correlated DE genes and store as .txt files
scriptDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_dist_vs_CPMcorrelation
Rscript ${scriptDir}/utility_find_high_and_low_CPM_correlated_genes_for_target_gene.R ${TF1}

###############################################################################################################################################################################

################# TF1 containing Enhancers of genes that are highly correlated with TF1 ##################
### For highly correlated query genes, fetch their enhancer regions, make bed file
# step 1. For TF1 highly correlating DE genes, get their enhancer regions, and find the regions with high confidence that TF1 occurs.
regionF=${dict_gene_enhancer}
geneF=${dir}/${TF1}_highly_correlated_genes.txt # generated by /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/utility_correlate_two_genes_CPM.R
regionTag=${TF1}HighCorrelationGenesEnhancers # syntax: camel
grep -f ${geneF} ${regionF} | awk '{OFS="\t"}{print $1"_"$2"_"$3}' | sort | uniq | tr '_' '\t' | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > regions.bed
run_FIMO ${dir} ${TF1} regions.bed 100
mv ${dir}/${TF1} ${dir}/${TF1}_${regionTag}
mv ${dir}/regions.bed ${dir}/${TF1}_${regionTag}/regions.bed
cd ${dir}/${TF1}_${regionTag}
Rscript ${scriptDir}/utility_calculate_bh_padj.R ${dir}/${TF1}_${regionTag}/fimo_${TF1}.bed # this generates padj_fimo_${TF1}.bed
cat padj_fimo_${TF1}.bed | awk '{if ($10 <= 0.05) print $0}' > fimo_topRanked_${TF1}.bed
get_motif_center_and_flanking_region fimo_topRanked_${TF1}.bed
rm fimo.split.regions*.tsv temp*
cd ${dir}
rm ${TF1}.meme


################# TF1 containing Enhancers of genes that are poorly correlated with TF1 ##################
### For poorly correlated query genes, fetch their enhancer regions, make bed file
# step 1. For TF1 poorly correlating DE genes, get their enhancer regions, and find the regions with high confidence that TF1 occurs.
regionF=${dict_gene_enhancer}
geneF=${dir}/${TF1}_poorly_correlated_genes.txt # generated by /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/utility_correlate_two_genes_CPM.R
regionTag=${TF1}LowCorrelationGenesEnhancers # syntax: camel
grep -f ${geneF} ${regionF} | awk '{OFS="\t"}{print $1"_"$2"_"$3}' | sort | uniq | tr '_' '\t' | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > regions.bed
run_FIMO ${dir} ${TF1} regions.bed 100
mv ${dir}/${TF1} ${dir}/${TF1}_${regionTag}
mv ${dir}/regions.bed ${dir}/${TF1}_${regionTag}/regions.bed
cd ${dir}/${TF1}_${regionTag}
Rscript ${scriptDir}/utility_calculate_bh_padj.R ${dir}/${TF1}_${regionTag}/fimo_${TF1}.bed # this generates padj_fimo_${TF1}.bed
cat padj_fimo_${TF1}.bed | awk '{if ($10 <= 0.05) print $0}' > fimo_topRanked_${TF1}.bed
get_motif_center_and_flanking_region fimo_topRanked_${TF1}.bed
rm fimo.split.regions*.tsv temp*
cd ${dir}
rm ${TF1}.meme

#################
conda deactivate
