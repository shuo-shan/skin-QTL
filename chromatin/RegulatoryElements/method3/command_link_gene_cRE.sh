##### step 6a. create gene-promoter dictionary (negative distance means gene tss is upstream of promoter.)
# header: gene, promoter name, distance, celltype
celltype=MEL
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,2 | grep _promoter_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd | sort -k1,1 -k3,3n > dictionary_gene_promoter_links_${celltype}.txt
rm temp*

celltype=KRT
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,2 | grep _promoter_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd | sort -k1,1 -k3,3n > dictionary_gene_promoter_links_${celltype}.txt
rm temp*

celltype=FRB
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,2 | grep _promoter_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd | sort -k1,1 -k3,3n > dictionary_gene_promoter_links_${celltype}.txt
rm temp*

##### step 6b. create gene-enhancer dictionary
celltype=MEL
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,2 | grep _enhancer_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd > temp.dictionary_gene_enhancer_links_${celltype}.txt
rm temp*
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,3 | grep _enhancer_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd >> temp.dictionary_gene_enhancer_links_${celltype}.txt
cat temp.dictionary_gene_enhancer_links_${celltype}.txt | sort -k1,1 -k3,3n > dictionary_gene_enhancer_links_${celltype}.txt
rm temp*

celltype=KRT
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,2 | grep _enhancer_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd > temp.dictionary_gene_enhancer_links_${celltype}.txt
rm temp*
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,3 | grep _enhancer_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd >> temp.dictionary_gene_enhancer_links_${celltype}.txt
cat temp.dictionary_gene_enhancer_links_${celltype}.txt | sort -k1,1 -k3,3n > dictionary_gene_enhancer_links_${celltype}.txt
rm temp*

celltype=FRB
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,2 | grep _enhancer_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd > temp.dictionary_gene_enhancer_links_${celltype}.txt
rm temp*
cat ${celltype}_atac_peak_and_neighboring_genes.txt | cut -d' ' -f1,3 | grep _enhancer_  > temp1
cat temp1 | cut -d' ' -f2 | cut -d'_' -f4 > tempa
cat temp1 | cut -d' ' -f1 > tempb
cat temp1 | cut -d' ' -f2 | cut -d'_' -f3,5 | tr '_' '\t' | awk '{if ($1=="upstream") print -$2; else print $2}' > tempc
cat temp1 | cut -d' ' -f2 | cut -d'_' -f1 > tempd
paste tempa tempb tempc tempd >> temp.dictionary_gene_enhancer_links_${celltype}.txt
cat temp.dictionary_gene_enhancer_links_${celltype}.txt | sort -k1,1 -k3,3n > dictionary_gene_enhancer_links_${celltype}.txt
rm temp*

