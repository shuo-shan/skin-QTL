import pysam
import collections
import argparse
import itertools


#this code is the optimized version of collapseDup_fast, where memory usage is minimzed 
#the v2 version:
#counting is done by aligning the 5' end of each reads, instead of leftmost coordinate in the aligned chromosome ( which is reference_start tag)
#also the PCR duplicates are written out for quality controls


def parseHeader(header):
    # extract the sequences that should be the UMI and BC:
    tmp1 = header.split(" ")
    tmp2 = tmp1[0].split(':')
    #use hard index to improve performance
    return tmp2[-2] + tmp2[-1]

def ToDict(key, value): 
    #convert a list to dict with value and key in seperate lists
    #the length of key and value must be the same
    res_dct = {}
    #update the value if the key exist
    for i in range(len(key)):
        if res_dct.has_key(key[i]):
            res_dct[key[i]].append(value[i])
        else:
            res_dct[key[i]]=[value[i]]
    return res_dct 

ap = argparse.ArgumentParser()
ap.add_argument('--inbam', help='Mapped reads in bam format')
ap.add_argument('--outbam', help='Output file to save reads after collapsing PCR duplicates')
#ap.add_argument("--u",help="a binary flag used to indicate that only uniquely mapped reads will be considered. By default uniquely mapped reads are defined as reads with MAPQ=60",action="store_true")
#ap.add_argument('--end',action='store_true',help='End position of read will not be considered (do NOT use False here! function not supported!')


args = ap.parse_args()
bam=args.inbam
out=args.outbam


chr_list=[]
position=[]
loaded_reads=[]
loaded_reads_flag=[]

samfile = pysam.Samfile(bam, "rb" )

#extract chr names frm BAm file
chr_list=[]
for i in samfile.header['SQ']:
    chr_list.append(i['SN'])


print ("List of chromosomes extracted from BAM")
print (chr_list)

mappedReads_counter=0
readSet_counter=0


bam_header = pysam.Samfile(bam, 'rb').header
outfile = pysam.Samfile(out, "wb", header=bam_header)
PCR_dup_out=pysam.Samfile('removed_PCR_dup.bam', "wb", header=bam_header)
print ("Open ",bam, "using pysam")


for chr in chr_list:
    #initial variables for each choromosome
    softClip=0
    position[:]=[]
    loaded_reads[:]=[]
    loaded_reads_flag[:]=[]
    numberReadsUniquePlusMultiMapped = 0
    
    for read in samfile.fetch(chr):
        mappedReads_counter+=1
        numberReadsUniquePlusMultiMapped+=1
        #locate the 5' end  a.get_reference_positions(full_length=True)
        if read.is_reverse:#5' end is the largest coordinate (original read is reverse complimented)
            if read.cigartuples[-1][0] == 4: #if the start of a read (5') is soft clipped
                myPosition=read.reference_end + read.cigartuples[-1][1] #always mark where the 5' end of the original read exactly mapped to 
                position.append(myPosition)
                loaded_reads.append(read)
                loaded_reads_flag.append('R')
                softClip += 1
            else:
                myPosition=read.reference_end
                position.append(myPosition)
                loaded_reads.append(read)
                loaded_reads_flag.append('R')
        else:#if the original reads is the same as forward genome (not reversed complimented)
            if read.cigartuples[0][0] == 4: #if the start of a read (5') is soft clipped
                myPosition=read.reference_start - read.cigartuples[0][1] #always mark where the 5' end of the original read exactly mapped to 
                position.append(myPosition)
                loaded_reads.append(read)
                loaded_reads_flag.append('F')
                softClip += 1
            else:
                myPosition=read.reference_start
                position.append(myPosition)
                loaded_reads.append(read)
                loaded_reads_flag.append('F')

    print ("numberReadsUniquePlusMultiMapped",numberReadsUniquePlusMultiMapped)
    print ("softClippedReads",softClip)
    counter_chr=collections.Counter(position) #count the frequency of all real position
    print ("Number of position with #reads sharing >=1", len(set(position)))

    count=0
    print  ("Processing", len(counter_chr.items()), "items")
    reads_dict = ToDict(position,loaded_reads) #save all reads to memory for speed!
    flag_dict = ToDict(position,loaded_reads_flag)#important to keep loaded reads and flag aligned!
    #free up memory 
    position[:]=[]
    loaded_reads[:]=[]
    loaded_reads_flag[:]=[]
    #start counting here
    for key,val in counter_chr.items():
        if count%10000==1:
            print (count)
        count+=1
        if val==1:#only one real reads, no need for checking UMI
            tmp = reads_dict[key]
            for tmp_read in tmp:
                 outfile.write(tmp_read)
                 readSet_counter+=1
                 
        elif val>1:#more than one real reads, check UMI
            notsetReads=set()
            notsetReads.clear()
            val_count=0
            readsCollection = reads_dict[key]
            flagCollection = flag_dict[key]
            for read, flag in itertools.izip(readsCollection,flagCollection): 
                #the clipped reads are only considered by their real position
                extended_read_name=parseHeader(read.query_name)+flag
                val_count+=1
                if extended_read_name not in notsetReads:
                    outfile.write(read)
                    notsetReads.add(extended_read_name)
                    readSet_counter+=1
                else:
                    PCR_dup_out.write(read)
            #error checking
            if val_count!=val:
                raise Exception('val_count != val')
outfile.close()
PCR_dup_out.close()



#-----------------------
#statistics


print ('Number of mapped reads',mappedReads_counter)
#print ('Number of reads mapped to unique location (UNIQUE map)',numberReadsUniqueGlobal)
print ('Number of reads after collapsing PCR dublicated (each read is present once)',readSet_counter)



samfile.close()
print ("DONE!")


