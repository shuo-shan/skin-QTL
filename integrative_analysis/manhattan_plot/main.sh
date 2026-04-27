
# step 1: get SNPs tested for PBS-eQTL
script_step1.sh  

# step 2: pick best association per SNP, and annotate SNP with genome location
python script_step2.py

# step 3: annotate SNP with genome location
join -t $'\t' -1 3 -2 1 -o 1.3,1.1,1.2,2.2,2.3 <(sort -k 3,3 temp_dictionary.txt) <(sort -k 1,1 best_associated_PBSeQTL_pairs_and_pval.txt) > best_associated_PBSeQTL_pairs_and_pval_and_position.txt

# step 4. plot Manhattan plot
echo "source activate fastQTL; Rscript script_step3.R" | bsub -J manhattan -W 8:00 -n 1 -R "span[hosts=1]" -R rusage[mem=408000] -q long -e "./%J%I.err" -o "./%J%I.out"


