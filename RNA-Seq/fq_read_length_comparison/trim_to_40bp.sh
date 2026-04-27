# trim the read 2 end 60 bases so the resulting 40bp fragment matches the length of previous libraries
# this script runs on 8 nodes and 1000MB memory per node
#!/bin/bash

# Load required module
module load cutadapt

# Input file
f_in=$1
f_out=$2

# Step 1: Trim the last 50 bases using cutadapt
date
cutadapt --cores=8 -u -60 -o ${f_out} ${f_in}
echo "done trimming ${f_in}"; date

