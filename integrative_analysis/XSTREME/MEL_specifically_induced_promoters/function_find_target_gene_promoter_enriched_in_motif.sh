#!/bin/bash

TF=$1 # transcription factor name
dir=$2 # where I want the result to be
regions=$3 # full path of BED file # promoter bed file to query from
fimothreshold=0.00001
regionSplitN=1000
dict_linked_gene_cRE=/pi/manuel.garber-umw/human/skin/eQTLs/chromatin/RegulatoryElements/method4/promoters_and_enhancers_surrounding_genes_all3cts.bed # gene-promoter-surrounding-enhancer link


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
                fimo --oc ${thisDir}/${thisTF} --thresh ${fimothreshold} ${thisDir}/${thisTF}.meme ${thisDir}/${f}
                mv ${thisDir}/${thisTF}/fimo.tsv ${thisDir}/${thisTF}/fimo.${outF}.tsv
        done
        cd ${thisDir}/${thisTF}
        cat fimo.*.tsv | grep chr | cut -f3,4,5,6,7,8,10 | sort -r -n -k5 > fimo.txt
        cat fimo.txt | cut -f1 > temp1
        cat fimo.txt | cut -f2,3,4,5,6,7 > temp2
        cat temp1 | sed 's/:/\t/g' | sed 's/-/\t/g' > tempa
        paste tempa temp2 | sort -r -n -k7 | awk '!seen[$1"_"$2"_"$3]++' > fimo_${thisTF}.bed
        mv fimo.txt fimo_output_${thisTF}.txt
	rm ${thisDir}/${thisTF}.meme
	rm *split*.tsv *xml temp1 temp2 tempa
        rm -r ${thisDir}/fasta
        cd ${thisDir}
}


#### resources

#### MAIN
# generates a fimo_${TF}.bed with all regions containing motif occurence with qval <= threshold
run_FIMO ${dir} ${TF} ${regions} ${regionSplitN}
# overlap with promoter bed
bedtools intersect -wb -a ${dir}/${TF}/fimo_${TF}.bed -b ${dict_linked_gene_cRE} | awk '$17=="promoter"' | awk '$16=="MEL"' | awk '{print $15}' | sort | uniq > ${dir}/FIMO_predicted_genes_with_${TF}_in_promoter.txt
# clean-up
cd ${dir}
rm -r ${dir}/${TF}













