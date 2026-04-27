#!/bin/bash

snp=$1
dir=$2

cd ${dir}

# acquire lock before writing
outputName=best_associated_PBSeQTL_pairs_and_pval
output_file="${dir}/${outputName}.txt"
lock_file="${dir}/${outputName}.lock"
exec 9>"$lock_file" # this line opens the lock file for writing and associates it with file descriptor 9. lock file is created if it does not already exist.
flock 9 # apply an advisory lock on this file

# write to the output file
# column 15 is featureSelected_PBSeQTL_pval
cat temp_result.txt | awk -v v=${snp} '{OFS="\t"}{if ($1==v) print $1,$2,$15}' | sort -nk3 | head -1 >> ${output_file}

# release the lock
flock -u 9

