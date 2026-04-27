#!/bin/bash

f=$1
fname=$(basename "$f" .p2.fq.gz)

# Define the lock file
lockfile="adapter_counts.lock"

# Perform your processing on the `.fq.gz` file
zcat "$f" | \
    tee >(grep -c "TGGAATTCTCGG" > "a1_${fname}.txt") \
        >(grep -c "AGATCGGAAGAGC" > "a2_${fname}.txt") \
        >(grep -c "CTGTCTCTTATA" > "a3_${fname}.txt") \
        >(grep -c "@" > "a4_${fname}.txt") \
    > /dev/null  # Discard the output from zcat, we only care about the counts

  # Read the counts
  a1=$(cat a1_${fname}.txt)
  a2=$(cat a2_${fname}.txt)
  a3=$(cat a3_${fname}.txt)
  a4=$(cat a4_${fname}.txt)

# Attempt to acquire the lock and run the following commands if successful
{
  # Lock the file for exclusive write access
  flock 200 || { echo "Could not acquire lock, exiting."; exit 1; }

  # Append the results to the output file
  echo -e "${f}\t${a1}\t${a2}\t${a3}\t${a4}" >> adapter_counts.txt

  # Clean up temporary files
  rm a1_${fname}.txt a2_${fname}.txt a3_${fname}.txt a4_${fname}.txt

} 200>"$lockfile"  # File descriptor 200 used for the lock file

