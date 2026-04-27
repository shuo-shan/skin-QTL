#!/bin/bash


####  build chrompbnet singularity image
singularity pull docker://kundajelab/chrombpnet:latest   
#singularity build <singularity-image>.sif docker://<docker-image-name>:<tag>
#it creates a .sif file in your current directory (e.g., chrombpnet_latest.sif). Since you're skipping singularity build (which is fine here), you can immediately proceed to the singularity inspect command without doing anything else

####  chrompbnet_lastest.sif: check image
singularity inspect -d chrombpnet_latest.sif
#     bootstrap: docker
#     from: kundajelab/chrombpnet:latest   

#### fix a bug in the chrombpnet by creating a new singularity image
cp chrombpnet_latest.sif my_container.sif
# creates a writable "sandbox" directory using .sif as source. regular .sif images are read-only. A sandbox allows me to enter and modify the image's filesystem like a regular directory.
singularity build --sandbox my_container_sandbox my_container.sif
# opens an interactive shell inside the container, allowing me to work within the modified image. --writable flag allows me to edit files and make changes.
singularity shell --writable my_container_sandbox
cd my_container_sandbox/
cd scratch/chrombpnet/chrombpnet/training/
# change the line in predict.py
sed -i.bak 's/from scipy import nanmean, nanstd/from numpy import nanmean, nanstd/' predict.py


#### fix a bug in the chrombpnet by creating a new singularity image
> singularity build --sandbox my_container_sandbox my_container.sif
> singularity shell --writable my_container_sandbox
> cd /scratch/chrombpnet/chrombpnet/training
> sed -i.bak 's/from scipy import nanmean, nanstd/from numpy import nanmean, nanstd/' predict.py
> exit


# Run chrombpnet: Option 1
> bsub -q gpu -W 8:00 -R rusage[mem=100G] -o chrom.log -e error.log singularity exec --nv /home/azita.ghodssi-umw/resources/images/chrombpnet_latest_rev.sif python /home/azita.ghodssi-umw/gpu_test.py

# scipy, which chrombpnet uses, does check the OMP_NUM_THREADS environment variable
> export OMP_NUM_THREADS=10
> echo $OMP_NUM_THREADS

# set the run time to max of 30 days, 100G memory divided by 10 cores
> bsub -J chrombpnet -q gpu -n 10 -R "rusage[mem=10240] span[hosts=1]" -W 720:00 -o chrom.log -e error.log singularity exec --nv /home/azita.ghodssi-umw/resources/images/chrombpnet_latest_rev.sif chrombpnet pipeline -ibam /home/azita.ghodssi-umw/chrombpnet_tutorial/data/downloads/merged.bam -d "ATAC" -g /home/azita.ghodssi-umw/chrombpnet_tutorial/data/downloads/hg38.fa -c /home/azita.ghodssi-umw/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes -p /home/azita.ghodssi-umw/chrombpnet_tutorial/data/downloads/peaks_no_blacklist.bed  -n /home/azita.ghodssi-umw/chrombpnet_tutorial/data/output_negatives.bed -fl /home/azita.ghodssi-umw/chrombpnet_tutorial/data/splits/fold_0.json -b /home/azita.ghodssi-umw/chrombpnet_tutorial/bias_model/ENCSR868FGK_bias_fold_0.h5 -o /home/azita.ghodssi-umw/chrombpnet_tutorial/chrombpnet_model/


# Run chrombpnet interactively: Option 2 
# start the gpu session
#alias isessiongpu="bsub -Is -q gpu -W 8:00 -R rusage[mem=100G] /bin/bash"
> isessiongpu

# start the singularity image 
> singularity shell --nv chrombpnet_latest.sif

# check whether tensorflow is installed
> python
> import tensorflow as tf
> print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
# should see "Num GPU Available: 1"
> exit()

> cd chrombpnet_tutorial/
> chrombpnet pipeline -ibam ./data/downloads/merged.bam -d "ATAC" -g ./data/downloads/hg38.fa -c ./data/downloads/hg38.chrom.sizes -p ./data/downloads/peaks_no_blacklist.bed  -n ./data/output_negatives.bed -fl ./data/splits/fold_0.json -b ./bias_model/ENCSR868FGK_bias_fold_0.h5 -o chrombpnet_model/
