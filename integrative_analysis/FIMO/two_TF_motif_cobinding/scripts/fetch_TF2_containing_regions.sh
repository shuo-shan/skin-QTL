#!/bin/bash
export PATH=/home/shuo.shan-umw/meme-5.5.1/meme/bin:/home/shuo.shan-umw/meme-5.5.1/meme/libexec/meme-5.5.1:$PATH
module load bedtools/2.30.0
TF1=$1
TF2=$2

#### resources
cRE_bed=/pi/manuel.garber-umw/human/skin/eQTLs/ATACseq_DolphinNext/peaks/merged_all_skin-eQTL_ATACseq_files_summits_300bp_flanking_window.bed
dict_gene_promoter=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method3/dictionary_gene_promoter_links_allcts.bed
scriptDir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_cobinding/scripts/
source activate fastQTL

#### set-up
dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/FIMO/two_TF_motif_cobinding/results/${TF1}/${TF2}
mkdir ${dir}
cd ${dir}

#### FUNCTIONS
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

#### function 4
function run_FIMO {
        echo "runs FIMO for the TF and regions of interest"
        thisTF=$1
        thisRegions=$2
        export PATH=/home/shuo.shan-umw/meme/bin:/home/shuo.shan-umw/meme/libexec/meme-5.5.5:$PATH
        module load bedtools/2.30.0
        get_meme_for_TF ${thisTF}                       # run custom function
        bedtools getfasta -fi /share/GHPCC/data/umw_biocore/genome_data/human/hg38_gencode_v34/hg38_gencode_v34.fa -bed ${thisRegions} > ${thisRegions}.fa
        # new FIMO version (5.5.5) reads .fa sequence name wrong. >chr10:129963242-129963543 doesn't result in sequence name.
        # replace >chr10:129963242-129963543 to >chr10_129963242_129963543
        cat ${thisRegions}.fa | tr ':' '-' | tr '-' '_' > renamed_${thisRegions}.fa
        fimo -o ./ --thresh 1 --text ${thisTF}.meme renamed_${thisRegions}.fa > fimo_all_output.txt
        #motif_id       motif_alt_id    sequence_name   start   stop    strand  score   p-value q-value matched_sequence
        rm ${thisTF}.meme ${thisRegions}.fa renamed_${thisRegions}.fa
}

######### REFERENCE

# folder: ${regionTag}
# input file: regions.bed
# output files: fimo_topRanked_${TF1}.bed, padj_fimo_${TF1}.bed, fimo_topRanked_${TF1}_motif_center.bed, fimo_topRanked_${TF1}_motif_center_flanking_1kb.bed
mkdir ${regionTag}
cd ${regionTag}
grep -w -f ${geneF} ${regionF} | awk '{OFS="\t"}{print $1"_"$2"_"$3}' | sort | uniq | tr '_' '\t' | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > regions.bed
### From the promoters, fetch the ones that TF1 motif occurs
run_FIMO ${TF1} regions.bed
Rscript ${scriptDir}/utility_calculate_bh_padj.R ${dir}/${regionTag}/fimo_all_output.txt # this fills the q-value column in the same file.
cat fimo_all_output.txt | awk '{OFS=FS="\t"}NR==1{print $0}NR>1{if ($8 <= 0.00001) print $0}' > fimo_sig_output.txt
cat fimo_sig_output.txt | awk 'NR>1' | cut -f3 > temp1
cat fimo_sig_output.txt | awk 'NR>1' | cut -f4,5,6,7,8,9 > temp2
cat temp1 | sed 's/_/\t/g' > tempa
paste tempa temp2 | sort -r -n -k7 | awk '!seen[$1"_"$2"_"$3]++' > fimo_sig_${TF1}.bed
get_motif_center_and_flanking_region fimo_sig_${TF1}.bed
rm temp1 temp2 tempa
cd ${dir}

##########








################# HIGH-Cor: TF2 in TF1-containing-Promoters of genes that are highly correlated with TF1 ##################
regionF=${dir}/../${TF1}_${TF1}HighCorrelationGenesPromoters/fimo_topRanked_${TF1}_motif_center_flanking_1kb.bed
# ^ first 3 columns: 1kb flanking region of the motif center matched in the region. column4: region location and motif matching start to finish. column5: score. column6: strand
regionTag=${TF1}TopMatching${TF1}HighCorrelationGenesPromoters
mkdir ${regionTag}
cd ${regionTag}
cat ${regionF} | awk '{OFS="\t"}{print $1"_"$2"_"$3"_"$6}' | sort | uniq | tr '_' '\t' | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$4}' > regions.bed
# ^ first 3 columns: 1kb flanking region of the TF1 motif center matched in region of interest, column4: unique ID, column5: strand
### From the TF1 centered 1kb window, fetch the ones that TF2 motif occurs
run_FIMO ${TF2} regions.bed 50
Rscript ${scriptDir}/utility_calculate_bh_padj.R ${dir}/${regionTag}/fimo_${TF2}.bed # this generates padj_fimo_${TF2}.bed
cat padj_fimo_${TF2}.bed | awk '{if ($8 <= 0.0001) print $0}' > fimo_topRanked_${TF2}.bed
get_motif_center_and_flanking_region fimo_topRanked_${TF2}.bed
cd ${dir}

################# LOW-Cor: TF2 in TF1-containing-Promoters of genes that are poorly correlated with TF1 ##################
regionF=${dir}/../${TF1}_${TF1}LowCorrelationGenesPromoters/fimo_topRanked_${TF1}_motif_center_flanking_1kb.bed
regionTag=${TF1}TopMatching${TF1}LowCorrelationGenesPromoters
mkdir ${regionTag}
cd ${regionTag}
cat ${regionF} | awk '{OFS="\t"}{print $1"_"$2"_"$3"_"$6}' | sort | uniq | tr '_' '\t' | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$4}' > regions.bed
### From the TF1 centered 1kb window, fetch the ones that TF2 motif occurs
run_FIMO ${TF2} regions.bed 50
Rscript ${scriptDir}/utility_calculate_bh_padj.R ${dir}/${regionTag}/fimo_${TF2}.bed # this generates padj_fimo_${TF2}.bed
cat padj_fimo_${TF2}.bed | awk '{if ($8 <= 0.0001) print $0}' > fimo_topRanked_${TF2}.bed
get_motif_center_and_flanking_region fimo_topRanked_${TF2}.bed
cd ${dir}


####### calculate motif center distance between TF1 and TF2 for both HIGH-COR and LOW-COR
#######
################# check for the occurence of TF1 and TF2 cobidning 
Rscript ${scriptDir}/utility_check_two_TF_cobinding.R ${TF1} ${TF2}

####### clean-up
cd ${dir}
#rm -r ${TF2}
conda deactivate


