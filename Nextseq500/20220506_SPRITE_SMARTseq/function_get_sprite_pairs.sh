#!/bin/bash
# this function outputs a BED file with all interaction partners of a 5kb tile of interest. number of how many clusters this pair is seen together is logged as integers.
# adapted from Pranitha Vangala's script listed below, by shuo.shan@umassmed.edu, in June 2022

dir=/nl/umw_manuel_garber/human/skin/eQTLs/Nextseq500/20220506_SPRITE_SMARTseq
gene=$1
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
