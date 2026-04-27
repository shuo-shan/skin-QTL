module load bedtools/2.29.2
f1=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_FRB_merged_peaks.bed
f2=/nl/umw_manuel_garber/human/skin/eQTLs/ChIP-seq/DolphinNext/FRB/report5357/chip/ChIP_FRB_merged.narrowPeak
#ln -s ${f1} ./
#ln -s ${f2} ./
bedtools intersect -a ${f1} -b ${f2} | awk '{OFS="\t"}{print $1,$2,$3,"ATAC_H3K27ac_intersect_FRB_"NR}' > temp.ATAC_H3K27acChIP_intersect_FRB.bed
bedtools sort -i temp.ATAC_H3K27acChIP_intersect_FRB.bed > ATAC_H3K27acChIP_intersect_FRB.bed

f1=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_KRT_merged_peaks.bed
f2=/nl/umw_manuel_garber/human/skin/eQTLs/ChIP-seq/DolphinNext/KRT_merged/report5430/chip/ChIP_KRT_merged.narrowPeak
#ln -s ${f1} ./
#ln -s ${f2} ./
bedtools intersect -a ${f1} -b ${f2} | awk '{OFS="\t"}{print $1,$2,$3,"ATAC_H3K27ac_intersect_KRT_"NR}' > temp.ATAC_H3K27acChIP_intersect_KRT.bed
bedtools sort -i temp.ATAC_H3K27acChIP_intersect_KRT.bed > ATAC_H3K27acChIP_intersect_KRT.bed

f1=/nl/umw_manuel_garber/human/skin/eQTLs/ATACseq_DolphinNext/softlinks/narrowPeak/ATAC_MEL_merged_peaks.bed
f2=/nl/umw_manuel_garber/human/skin/eQTLs/ChIP-seq/DolphinNext/MEL/report5356/chip/ChIP_MEL_merged.narrowPeak
#ln -s ${f1} ./
#ln -s ${f2} ./
bedtools intersect -a ${f1} -b ${f2} | awk '{OFS="\t"}{print $1,$2,$3,"ATAC_H3K27ac_intersect_MEL_"NR}' > temp.ATAC_H3K27acChIP_intersect_MEL.bed
bedtools sort -i temp.ATAC_H3K27acChIP_intersect_MEL.bed > ATAC_H3K27acChIP_intersect_MEL.bed

### DEFINE TSS from Ensembl GTF file
gene=/nl/umw_manuel_garber/human/skin/eQTLs/literature/UCSC_tracks/Homo_sapiens.GRCh38.105.gtf.gz
zcat $gene | awk '{if ($3=="transcript") print $0}' | awk '{OFS="\t"}{if ($7=="+") print $4,$0; else print $5,$0}' > temp.output
cat temp.output | cut -f2,1,8,10 | awk '{OFS=FS="\t"}{print $2,$1,$3,$4}' | awk '{if ($0~/gene_name/) print "chr"$0}' > temp.output2
cat temp.output2 | cut -f1,2,3 | awk '{OFS="\t"}{print $1,$2,$2+1,$3}' > temp.output3 # end position added one to match bed format. Gtf format end position is included. Bed format end position is NOT included.
cat temp.output2 | sed 's/.*gene_id \"//g' | sed 's/\"; .*//g' > temp.id
cat temp.output2 | sed 's/.*gene_name \"//g' | sed 's/\"; .*//g' > temp.name
cat temp.output2 | sed 's/.*transcript_id \"//g' | sed 's/\"; .*//g' > temp.transcriptid
cat temp.output2 | sed 's/.*gene_biotype \"//g' | sed 's/\"; .*//g' > temp.biotype
paste -d"_" temp.id temp.transcriptid > temp.id2
paste -d"_" temp.name temp.transcriptid > temp.name2

paste temp.output3 temp.id2 temp.name2 temp.name temp.biotype > Ensembl_GRCh38.105_TSS_annotation.txt
bedtools sort -i Ensembl_GRCh38.105_TSS_annotation.txt > Ensembl_GRCh38.105_TSS_annotation_sorted.txt

### annotate each RE with proximity to closest TSS
# columns: RE chromosome, RE start, RE end (half open, half closed), RE id, TSS chromosome, TSS start, TSS end, TSS strand, TSS_ensembl_gene&transcript_id, TSS_Ensembl_genesymbol&transcript id, TSS_Ensembl_geneName, TSS_Ensembl_gene_biotype, distance 
module load bedtools/2.29.2
bedtools closest -t first -a ATAC_H3K27acChIP_intersect_FRB.bed -b Ensembl_GRCh38.105_TSS_annotation_sorted.txt -d > ATAC_H3K27acChIP_intersect_FRB_with_TSS_distance.txt
bedtools closest -a ATAC_H3K27acChIP_intersect_KRT.bed -b Ensembl_GRCh38.105_TSS_annotation_sorted.txt -d > ATAC_H3K27acChIP_intersect_KRT_with_TSS_distance.txt
bedtools closest -a ATAC_H3K27acChIP_intersect_MEL.bed -b Ensembl_GRCh38.105_TSS_annotation_sorted.txt -d > ATAC_H3K27acChIP_intersect_MEL_with_TSS_distance.txt
