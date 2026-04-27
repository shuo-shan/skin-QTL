#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R rusage[mem=20400]
#BSUB -q long
#BSUB -W 121:00
#BSUB -e "./log/%J%I.join.err"
#BSUB -o "./log/%J%I.join.out"

rm commands.KRT.join.txt commands.MEL.join.txt commands.FRB.join.txt
celltype=KRT
cat temp_atac_active_in_krt_pbs.txt temp_atac_active_in_krt_ifn.txt | sort | uniq > temp.peaks
while read peak;do
  echo "bash function_join_table.sh ${celltype} ${peak}" >> commands.${celltype}.join.txt 
done < temp.peaks
echo "done with ${celltype}"
rm temp.peaks

celltype=MEL
cat temp_atac_active_in_mel_pbs.txt temp_atac_active_in_mel_ifn.txt | sort | uniq > temp.peaks
while read peak;do
  echo "bash function_join_table.sh ${celltype} ${peak}" >> commands.${celltype}.join.txt
done < temp.peaks
echo "done with ${celltype}"
rm temp.peaks

celltype=FRB
cat temp_atac_active_in_frb_pbs.txt temp_atac_active_in_frb_ifn.txt | sort | uniq > temp.peaks
while read peak;do
  echo "bash function_join_table.sh ${celltype} ${peak}" >> commands.${celltype}.join.txt
done < temp.peaks
echo "done with ${celltype}"
rm temp.peaks

# append the bjob header to each commands.txt file and submit them as bjob in the long queue
bsub -J joinKRT < commands.KRT.join.txt
bsub -J joinMEL < commands.MEL.join.txt
bsub -J joinFRB < commands.FRB.join.txt
