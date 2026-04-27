# list of commands to run jaspar enrichment
# useful links: https://bitbucket.org/CBGR/jaspar_enrichment/src/master/, https://jaspar.genereg.net/enrichment/

# Enrichment within a given universe of genomic regions
jaspar_dir=/pi/manuel.garber-umw/human/skin/eQTLs/integrative_analysis/JASPAR/jaspar_enrichment/bin
${jaspar_dir}/JASPAR_enrich.sh oneSetBg <loladb_dir> <S bed> <U bed> <output dir> <API_URL> <n_cores>
# <loladb_dir> is the directory containing all the LOLA database batches for an assembly
# <API_URL> is the URL to the matrix API in JASPAR, e.g. http://jaspar.genereg.net/api/v1/matrix/
# <n_cores> is the number of cores to use when parallelizing.
# <S bed> is the genomic regions we are interested in, .bed format
# <U bed> is the universie set of genomic regions. every region in S should overlap with one region in U.

# example 1:
# this gives us a list of enriched TFs in the regions. You can make a beeswarm plot with the -log10(P-value) and color code the TF class
${jaspar_dir}/JASPAR_enrich.sh oneSetBg ${hg38_dir} promoters.bed all_regions.bed ${output_dir}



# Differential Enrichment
${jaspar_dir}/JASPAR_enrich.sh twoSets <loladb_dir> <S1 bed> <S2 bed> <output dir> <API URL> <n_cores>
# the <U bed> is not needed here because JASPAR treats the U as the merge between S1 and S2 bed files.
