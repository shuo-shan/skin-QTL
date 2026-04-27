#!/bin/bash
# adapted from Pranitha Vangala's script: /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/perElement_contacts_5kb_2-100.sh
# turns SPRITE pipeline output: a cluster file, into a folder containing all the pairwise interactions. Each file is a 5kb tile across genome, and contains all the interactions with this 5kb tile.  

############## my adapted script
######## Section 1. Prepare the 4 input files needed for this script
# step 1. make subset of cluster file with clusters of size 2 to 100
dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
cat ${dir}/SPRITE-F37-KRT-PBS.clusters | awk 'NF>2' | awk 'NF<100' | grep -v 'NOT_FOUND' > ${dir}/SPRITE-F37-KRT-PBS.clusters_2-100 # results in 16,882 clusters

# step 2. make 5kb genome tiles
module load bedtools/2.29.2
bedtools makewindows -g /share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes -w 5000 > ${dir}/hg38.5kb.windows.bed

# step 3. make perLocus_5kb_clusters file
rm barcoded_elements.bed
while read line;do 
  tag=$(echo ${line} | cut -d' ' -f1); 
  echo ${line} | cut -d' ' -f2- | tr ' ' '\n' | sed 's/.*_//g' | sed 's/:/\t/g' | sed 's/-/\t/g' > tmp.this; 
  cat tmp.this | awk -v t=${tag} '{OFS="\t"}{print $0,t}' - >> barcoded_elements.bed
  rm tmp.this; 
done < ${dir}/SPRITE-F37-KRT-PBS.clusters_2-100 # for each fragment with barcode ID in cluster file, list it as BED format with the barcode ID being its name.
bedtools sort -i ${dir}/barcoded_elements.bed > ${dir}/barcoded_elements_sorted.bed
bedtools intersect -a ${dir}/hg38.5kb.windows.bed -b ${dir}/barcoded_elements_sorted.bed -wo -F 0.50 > tmp1
cat tmp1 | awk '{OFS="\t"}{print $0,$4"_"$5"_"$6"_"$7}' > tmp2
cat tmp2 | awk '!seen[$9]++' > tmp3
cat tmp3 | awk '{OFS="\t"}{print $1"_"$2"_"$3,$7}' > ${dir}/perLocus_5kb_clusters # remove duplicates (in case fragment overlaps 50:50 in two tiles, keep the first tile seen)
rm tmp1 tmp2 tmp3
rm ${dir}/barcoded_elements.bed
# this code chunk list 5kb genome tiles that intersect sprite-barcoded fragments, one tile : one fragment per line.

# step 4. annotate 5kb tiles that overlap putative regulatory regions defined by k27ac peaks
# for k27ac peak definition, see /nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/define_RE.sh script
k27acPeaks=/nl/umw_manuel_garber/human/skin/eQTLs/chromatin/RegulatoryElements/KRT_k27acPeaks_sorted.bed
bedtools intersect -a ${dir}/hg38.5kb.windows.bed -b ${k27acPeaks} -wo -F 0.50 > tmp1
cat tmp1 | awk '!seen[$7]++' > tmp2 # if element equally spans two 5kb tiles, keep the first one that's seen.
cat tmp2 | awk '{OFS="\t"}{print $1,$2,$3,$7}' | sed 's/;/\t/g' | awk '{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3,$7,$8}' > ${dir}/k27ac_peaks_annotations_5k.txt 
rm tmp1 tmp2

# step 5. focus on the 5kb tiles that overlap a gene body
cat /nl/umw_manuel_garber/human/skin/eQTLs/literature/UCSC_tracks/Ensembl_GRCh38.105_biotype_annotation.txt | awk '{OFS="\t"}{print $1,$2,$3}' > tmp.genes
cat ${dir}/perLocus_5kb_clusters | cut -f1 | tr '_' '\t' | awk '{OFS="\t"}{print $0,$1"_"$2"_"$3}' > tmp1
bedtools intersect -a tmp1 -b tmp.genes -u > ${dir}/perLocus_5kb_clusters_overlapping_gene # 36,341 tiles out of 69,250 

######## Section 2. Main body of the script
# I will turn this into a function in the same folder.
while read s;do
  gene=$(echo ${s} | cut -d' ' -f4)
  echo "bash function_get_sprite_pairs.sh ${gene}" >> commands.txt
done < ${dir}/perLocus_5kb_clusters_overlapping_gene

jobname=sprite
split -l 300 commands.txt temp.commands.
for f in temp.commands.*;do
  cat $f | tr '\n' '; ' | sed 's/$/\n/g' > ${f}.joined
done 
cat *.joined > commands.joined.txt
rm temp.commands.* commands.txt

while read c; do
  echo ${c} | bsub -W 8:00 -J ${jobname} -n 1 -R "span[hosts=1]" -R rusage[mem=2040] -R "select[rh=6]" -q long -e "./log/step7.%J%I.err" -o "./log/step7.%J%I.out"
done < commands.joined.txt
sleep 5
while [[ $(bjobs | grep ${jobname} | wc -l) != 0 ]] ; do echo $(bjobs | grep ${jobname} | wc -l) "jobs remaining"; sleep 5; done

# append all files into one file all_cis
rm ${dir}/all_cis
cd ${dir}/perLocus_5k
for f in *.bed;do
  id=$(echo $f | sed 's/anno//g' | sed 's/\.bed//g')
  awk -v name=${id} '{OFS="\t"}{print $0,name}' ${f} >> ${dir}/all_cis
  echo "done with ${id}"
done
cd ${dir}

# example script for one gene:
#gene=$1
gene=chr16_27310000_27315000
outdir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
mkdir -p ${outdir}/perLocus_5k
cd ${outdir}/perLocus_5k
chr=$(grep -w ${gene} ${dir}/k27ac_peaks_annotations_5k.txt | awk '{print $1}' | sort -u)
cat ${outdir}/perLocus_5kb_clusters | grep -w ${gene} | awk '{print $2}' | sort -u > tmp_${gene} # find all sprite-barcode IDs that reside in the 5kb tile(s) containing the gene's promoter
fgrep -f tmp_${gene} ${dir}/SPRITE-F37-KRT-PBS.clusters_2-100 | awk 'NF>2' > ${gene}_clusters # find clusters for the gene's 5kb tile's sprite barcode ID
rm ${gene}_clusters.bed
while read s; do
  tag=$(echo ${s} | awk '{print $1}')
  echo ${s} | sed 's/ /\n/g' | grep -v "DPM" | sed 's/.*]_//g' | sed 's/:/\t/g' | sed 's/-/\t/g' | awk -v k=${tag} '{OFS="\t"}{print $0,k}' | awk '{OFS="\t"}{print $1,$2,$2+200,$4}' | grep -w ${chr} >> ${gene}_clusters.bed
done < ${gene}_clusters # for all fragments interacting with this gene's promoter-residing 5kb tile, fetch the cis-fragments in same chromosome, and extend the fragments to be 200bp, and list as a bed file for the gene's cis-interacting partners in the same cluster. if this promoter is in multiple clusters, append them all in the same BED file
bedtools sort -i ${gene}_clusters.bed | bedtools groupby -g 1 -c 2,3 -o min,max -i stdin > ${gene}_range.bed # if interacting fragments span multiple regions in the chromosome, take the start and end as the range and store as one line in BED format
bedtools intersect -wa -a ${dir}/hg38.5kb.windows.bed -b ${gene}_range.bed > ${gene}.promoter.bed # all 5kb tiles that overlap the gene's range of interacting fragments are considered as "promoters"??? <- it's just a poor choice for the name. it really is a way to narrow down the 5kb tiles that overlap the gene or gene's interaction partners.

# 1. get 5kb tiles that intersect the 200bp fragments in clusters containing the gene's promoter
# 2. get 5kb tiles, count number of distinct clusters occur, list cluster ID
# 3. get 5kb tiles in this cluster that intersect regulatory regions marked by k27ac, list reg region annotation
# 4. format it to bedfile, remove duplicated 5kb tiles
bedtools intersect -wa -wb -a ${gene}.promoter.bed -b ${gene}_clusters.bed |\
bedtools groupby -g 1,2,3 -i stdin -c 7 -o count_distinct,distinct |\
bedtools intersect -wa -wb -a stdin -b ${dir}/k27ac_peaks_annotations_5k.txt |\
awk -v OFS="\t" '{print $1,$2,$3,$NF,$4,$5}' | awk '!seen[$1,$2,$3]++' > anno${gene}.bed

# 1. get 5kb tiles that intersect the 200bp fragments in clusters containing the gene's promoter
# 2. subset the 5kb tiles to find ones that are not annotated with k27ac mark, annotate these tiles as "none"
# 3. format it to bedfile
bedtools intersect -wa -wb -a ${gene}.promoter.bed -b ${gene}_clusters.bed |\
bedtools groupby -g 1,2,3 -i stdin -c 7 -o count_distinct,distinct |\
bedtools intersect -wa -v -a stdin -b ${dir}/k27ac_peaks_annotations_5k.txt |\
awk -v OFS="\t" '{print $1,$2,$3,"none",$4,$5}' | awk '!seen[$1,$2,$3]++' >> anno${gene}.bed
cat ${dir}/k27ac_peaks_annotations_5k.txt | grep ${gene} | bedtools intersect -v -a anno${gene}.bed -b stdin > tmp_${gene} # discard this gene's promoter
bedtools sort -faidx /share/data/umw_biocore/dnext_data/genome_data/human/hg38/main/genome.chrom.sizes -i tmp_${gene} | awk '!seen[$1,$2,$3]++' > anno${gene}.bed # sort all 5kbp tiles that interact with the gene's promoter 

# 4. remove file if it's empty
nline=$(wc -l ./perLocus_5k/anno${gene}.bed | cut -d' ' -f1)
if [ ${nline} == 0 ] ; then
  echo anno${gene}.bed has zero lines.
  rm anno${gene}.bed
fi

rm tmp_${gene}
rm ${gene}_clusters
rm ${gene}_clusters.bed
rm ${gene}_range.bed
rm ${gene}.promoter.bed

# this file outputs a bed file that contains all the unique interacting 5kb tiles to a 5kb tile of interest. Each interacting 5kb tile is annotated by whether it overlaps a putative regulation region defined by k27ac peak (none, promoter, or enhancer). it also counts how many times the interaction partner is seen with this tile of interest (i.e. how many unique clusters do this pair appear).


############### REF SCRIPT
# adapted from Pranitha Vangala's script: /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/perElement_contacts_5kb_2-100.sh
#gene=$1
#outdir=/nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/
#mkdir -p ${outdir}/perLocus_5k
#cd ${outdir}/perLocus_5k
#chr=`cat /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/k27ac_peaks_annotations_5k |grep "\<$gene\>" | awk '{print $1}' | sort -u`
#grep "\<$gene\>" ${outdir}/perElement_5k_clusters | awk '{print $2}' |sort -u > tmp_${gene}
# v rate limiting step
#fgrep -f tmp_${gene} ${outdir}/clusters_hMDM_2-100 | awk 'NF>2' > ${gene}_clusters
#rm ${gene}_clusters.bed
#while read s; do tag=`echo ${s} | awk '{print $1}'`; echo ${s} | sed 's/ /\n/g'| grep -v "DPM"| awk -v k="$tag" '{print $1"\t"k}' | sed 's/:/\t/g' | awk -v OFS="\t" '{print $1,$2,$2+200,$3}' |grep -w $chr >> ${gene}_clusters.bed; done < ${gene}_clusters
#bedtools sort -i ${gene}_clusters.bed | bedtools groupby -g 1 -c 2,3 -o min,max -i stdin > ${gene}_range.bed
#bedtools intersect -wa -a /project/umw_garberlab/vangala/bin/hg19.5kb.win.bed -b ${gene}_range.bed > ${gene}.promoter.bed
#bedtools intersect -wa -wb -a ${gene}.promoter.bed -b ${gene}_clusters.bed | bedtools groupby -g 1,2,3 -i stdin -c 7 -o count_distinct,distinct |  bedtools intersect -wa -wb -a stdin -b /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/k27ac_peaks_annotations_5k | awk -v OFS="\t" '{print $1,$2,$3,$NF,$4,$5}' > anno${gene}.bed
#bedtools intersect -wa -wb -a ${gene}.promoter.bed -b ${gene}_clusters.bed | bedtools groupby -g 1,2,3 -i stdin -c 7 -o count_distinct,distinct |  bedtools intersect -wa -v -a stdin -b /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/k27ac_peaks_annotations_5k | awk -v OFS="\t" '{print $1,$2,$3,"none",$4,$5}' >> anno${gene}.bed
#cat /nl/umw_manuel_garber/human/hMDM/SPRITE/addDepth_09_07_19/split_fastqs/bams/k27ac_peaks_annotations_5k | grep "\<$gene\>"| bedtools intersect -v -a anno${gene}.bed -b stdin > tmp_${gene}
#bedtools sort -faidx /share/data/umw_biocore/genome_data/human/hg19/hg19.chrom.sizes -i tmp_${gene} | awk '!seen[$1,$2,$3]++' > anno${gene}.bed

### v this command was commented out in the original script, so it was probably not used
##bedtools intersect -sorted -wa -wb -a anno${gene}.bed -b clust_2-1000-final_SIP0523-*_merged.bam -g /share/data/umw_biocore/genome_data/human/hg19/hg19.chrom.sizes | sed 's/::/\t/g' | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$12}' | sed 's/]\[/\./g;s/]//g;s/\[//g' | bedtools groupby -sorted -g 1,2,3,4,5,6 -i stdin -c 7 -o count_distinct > perElement_interactions/anno${gene}.bed

#rm tmp_${gene}
#rm ${gene}_clusters
#rm ${gene}_clusters.bed
#rm ${gene}_range.bed
#rm ${gene}.promoter.bed
