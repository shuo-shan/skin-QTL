celltype=$1
peak=$2

cat ${celltype}_with_TSS_distance_annotated.txt | grep -w ${peak} | sort -k 7n > temp.${celltype}.${peak}
head -1 temp.${celltype}.${peak} | awk -v ct=${celltype} '{OFS="\t"}{print $4,ct"_"$9"_"$8"_"$6"_"$7}' > temp1.${celltype}.${peak}
tail -1 temp.${celltype}.${peak} | awk -v ct=${celltype} '{OFS="\t"}{print $4,ct"_"$9"_"$8"_"$6"_"$7}' > temp2.${celltype}.${peak}
join temp1.${celltype}.${peak} temp2.${celltype}.${peak} > temp3.${celltype}.${peak}
cat temp3.${celltype}.${peak} >> ${celltype}_atac_peak_and_neighboring_genes.txt
rm temp.${celltype}.${peak} temp1.${celltype}.${peak} temp2.${celltype}.${peak} temp3.${celltype}.${peak}
