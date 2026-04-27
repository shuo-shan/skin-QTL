# trim the read 2 end 50 bases to prevent polyA mapping to mRNA during alignment; modify the barcode from UMI_barcode format to UMI:barcode format for ESAT input
#!/bin/bash

# Load required module
module load cutadapt

# Set working directory
Dir=/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/01152025_shallowtest/fastq/skineQTL01152025shallowseq
cd ${Dir}

# Input file
fname=$1
f=${fname%.p2.fq.gz}

# Step 1: Trim the last 50 bases using cutadapt
cutadapt -u -50 -o "${f}.p2_trimmed.fq.gz" "${f}.p2.fq.gz"

# Print filename and timestamp
echo "done trimming ${f}.p2.fq.gz"; date

# Step 2: Modify the barcode in the header
zcat ${f}.p2_trimmed.fq.gz | awk '{if(NR%4==1) gsub("_", ":"); print}' | gzip > ${f}.p2_trimmed_bcmodified.fq.gz
echo "Processing complete: ${f}"; date

# cleanup
rm ${f}.p2_trimmed.fq.gz
