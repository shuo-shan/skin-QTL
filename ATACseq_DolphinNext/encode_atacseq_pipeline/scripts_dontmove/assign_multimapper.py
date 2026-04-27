#!/usr/bin/env python3

###### assign_multimappers.py
##
## When mapping with Bowtie2, a read may have multiple valid alignments (multimapping reads).
## This script ensures that reads with too many alignments (≥ -k) are filtered out, leaving only
## ambiguously mapped reads below the threshold.
##
## Logic:
##     Reads are grouped by query name (qname).
##        If a read aligns to more locations than the specified -k threshold, it is discarded.
##        If it has fewer or equal alignments, all alignments are retained.
##     Handles paired-end reads (--paired-end flag)
##        If paired-end mode is enabled, -k is doubled (since each fragment has two reads).
##     Processes a QNAME-sorted SAM file
##         Reads are stored in current_reads while processing.
##         If a new QNAME is encountered:
##           If the previous read has alignments exceeding -k, it is discarded.
##           Otherwise, all previous alignments are written to stdout.


#!/usr/bin/env python2

import sys
import random
import argparse
import pysam

def parse_args():
    '''
    Gives options
    '''
    parser = argparse.ArgumentParser(description='Saves reads below a alignment threshold and discards all others')
    parser.add_argument('-k', help='Alignment number cutoff', required=True)
    parser.add_argument('--paired-end', dest='paired_ended', action='store_true', help='Data is paired-end')
    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('output_bam', help='Output BAM file')
    args = parser.parse_args()
    
    alignment_cutoff = int(args.k)
    paired_ended = args.paired_ended

    return alignment_cutoff, paired_ended, args.input_bam, args.output_bam


if __name__ == "__main__":
    '''
    Runs the filtering step of choosing multimapped reads
    '''

    alignment_cutoff, paired_ended, input_bam, output_bam = parse_args()

    if paired_ended:
        alignment_cutoff = alignment_cutoff * 2

    # Open input and output BAM files using pysam
    input_samfile = pysam.AlignmentFile(input_bam, "rb")
    output_samfile = pysam.AlignmentFile(output_bam, "wb", header=input_samfile.header)

    current_reads = []  # Store reads with the same qname
    current_qname = ''

    # Write header only once at the beginning
    for read in input_samfile:
        if read.qname.startswith('@'):
            continue  # Skip header lines in the loop

        # Process the read
        read_elems = read.qname

        # Handle reads with the same qname
        if read.qname == current_qname:
            current_reads.append(read)
        else:
            # Discard reads that have more than the alignment cutoff
            if len(current_reads) >= alignment_cutoff:
                current_reads = [read]
                current_qname = read.qname
            elif len(current_reads) > 0:
                # Output all reads for the current qname, then discard
                for r in current_reads:
                    output_samfile.write(r)

                current_reads = [read]
                current_qname = read.qname
            else:
                # First read in file
                current_reads.append(read)
                current_qname = read.qname

    # Write any remaining reads to output
    for read in current_reads:
        output_samfile.write(read)

    # Close files
    input_samfile.close()
    output_samfile.close()


