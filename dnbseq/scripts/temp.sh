#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=120000]
#BSUB -q long
#BSUB -W 121:00
#BSUB -o “/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/scripts/temp.%J%I.out” 
#BSUB -e “/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq/scripts/temp.%J%I.err”
### comment

module load java/1.8.0_171
module load condas/2018-05-11
source activate sshan_isoform
Dir=/nl/umw_manuel_garber/human/skin/eQTLs/dnbseq

##############################################################
### split fastq by CS2 barcode
lib=F30M_F31M
barcodeMap=$Dir/data/data_info/barcodeMap/$lib
cd $Dir/data/data_info
cat samples | grep $lib > this.sample
cd $Dir/data/fastq/$lib
while read s;do
    gunzip ${s}_1.fq.gz && gunzip ${s}_2.fq.gz
    java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_1.fq F2=${s}_2.fq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
    gzip ${s}_1.fq && gzip ${s}_2.fq
    cd $Dir/data/fastq/$lib/$s
    for i in $(ls */*fq);do
      gzip $i
      echo $i
    done
    echo ${s}
done < $Dir/data/data_info/this.sample

lib=F41M_F45M
barcodeMap=$Dir/data/data_info/barcodeMap/$lib
cd $Dir/data/data_info
cat samples | grep $lib > this.sample
cd $Dir/data/fastq/$lib
while read s;do
    gunzip ${s}_1.fq.gz && gunzip ${s}_2.fq.gz
    java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_1.fq F2=${s}_2.fq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
    gzip ${s}_1.fq && gzip ${s}_2.fq
    cd $Dir/data/fastq/$lib/$s
    for i in $(ls */*fq);do
      gzip $i
      echo $i
    done
    echo ${s}
done < $Dir/data/data_info/this.sample

lib=F46M_F47M
barcodeMap=$Dir/data/data_info/barcodeMap/$lib
cd $Dir/data/data_info
cat samples | grep $lib > this.sample
cd $Dir/data/fastq/$lib
while read s;do
    gunzip ${s}_1.fq.gz && gunzip ${s}_2.fq.gz
    java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_1.fq F2=${s}_2.fq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
    gzip ${s}_1.fq && gzip ${s}_2.fq
    cd $Dir/data/fastq/$lib/$s
    for i in $(ls */*fq);do
      gzip $i
      echo $i
    done
    echo ${s}
done < $Dir/data/data_info/this.sample

lib=F49M_F50M
barcodeMap=$Dir/data/data_info/barcodeMap/$lib
cd $Dir/data/data_info
cat samples | grep $lib > this.sample
cd $Dir/data/fastq/$lib
while read s;do
    gunzip ${s}_1.fq.gz && gunzip ${s}_2.fq.gz
    java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_1.fq F2=${s}_2.fq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
    gzip ${s}_1.fq && gzip ${s}_2.fq
    cd $Dir/data/fastq/$lib/$s
    for i in $(ls */*fq);do
      gzip $i
      echo $i
    done
    echo ${s}
done < $Dir/data/data_info/this.sample

lib=F51M_F52M
barcodeMap=$Dir/data/data_info/barcodeMap/$lib
cd $Dir/data/data_info
cat samples | grep $lib > this.sample
cd $Dir/data/fastq/$lib
while read s;do
    gunzip ${s}_1.fq.gz && gunzip ${s}_2.fq.gz
    java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_1.fq F2=${s}_2.fq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
    gzip ${s}_1.fq && gzip ${s}_2.fq
    cd $Dir/data/fastq/$lib/$s
    for i in $(ls */*fq);do
      gzip $i
      echo $i
    done
    echo ${s}
done < $Dir/data/data_info/this.sample

lib=F55M_F56M
barcodeMap=$Dir/data/data_info/barcodeMap/$lib
cd $Dir/data/data_info
cat samples | grep $lib > this.sample
cd $Dir/data/fastq/$lib
while read s;do
    gunzip ${s}_1.fq.gz && gunzip ${s}_2.fq.gz
    java -Xmx100g -jar /home/ed70w/bin/splitter_09.21.15_13.57.jar F1=${s}_1.fq F2=${s}_2.fq B=NNNNNNSSSSSS M=$barcodeMap HD=1 O=${s}
    gzip ${s}_1.fq && gzip ${s}_2.fq
    cd $Dir/data/fastq/$lib/$s
    for i in $(ls */*fq);do
      gzip $i
      echo $i
    done
    echo ${s}
done < $Dir/data/data_info/this.sample
