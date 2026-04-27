# trim the read 2 end 50 bases to prevent polyA mapping to mRNA during alignment; modify the barcode from UMI_barcode format to UMI:barcode format for ESAT input
# this script runs on 8 nodes and 1000MB memory per node
#!/bin/bash

# Load required module
module load cutadapt

# Set working directory
Dir=/pi/manuel.garber-umw/human/skin/eQTLs/Azenta_sequencing/03312025_3ctk_celseq2_atac_batch1/fastq_celseq2

# Input file
fname=$1
f=${fname%.p2.fq.gz}

# Step 1: Trim the last 50 bases using cutadapt
date
cutadapt --cores=8 -u -50 -o "${f}.p2_trimmed.fq.gz" "${f}.p2.fq.gz"
echo "done trimming ${f}.p2.fq.gz"; date

# Step 2: Modify the barcode in the header
date
pigz -p 8 -dc ${f}.p2_trimmed.fq.gz > ${f}.p2_trimmed.fq
echo "done with unzipping!"; date
awk '{if(NR%4==1) gsub("_", ":"); print}' ${f}.p2_trimmed.fq | pigz -p 8 > ${f}.p2_trimmed_bcmodified.fq.gz
echo "Processing complete: ${f}"; date

# cleanup
rm ${f}.p2_trimmed.fq.gz ${f}.p2_trimmed.fq
