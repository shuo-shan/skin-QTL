dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot/eQTL_vs_GWAS
cd $dir
# GWAS results:
gwas=/pi/manuel.garber-umw/human/skin/eQTLs/GWAS_SNPs/vitiligo/JinY2016/all_GWAS123cmh.txt
# eQTL results:
# MEL PBSeQTL: /pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot/best_associated_PBSeQTL_pairs_and_pval_and_position.txt
eqtl=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/manhattan_plot/best_associated_PBSeQTL_pairs_and_pval_and_position.txt

# the follow code makes these plots:
# correlation plot of common SNPs. y-axis = GWAS -log10Pval. x-axis = eQTL -log10Pval
# get a list of common SNPs in GWAS and eQTL results. Plot GWAS manhattan plot, highlighting eQTLs and labeling eQTL pval
# get a list of common SNPs in GWAS and eQTL results. Plot eQTL manhattan plot, highlighting significant GWAS SNPs and labeling GWAS pval


cat $gwas | awk 'NR>1{print $2}' | sed 's/RS/rs/g' > gwas_tested_snps.txt
cat $eqtl | cut -f1 > eqtl_tested_snps.txt
comm -12 <(sort gwas_tested_snps.txt) <(sort eqtl_tested_snps.txt) > common_tested_snps.txt
cat $gwas | awk 'NR>1{if ($8 < 0.00000007) print $2}' | sed 's/RS/rs/g' > gwas_sig_snps.txt
# plotted in script01.R




